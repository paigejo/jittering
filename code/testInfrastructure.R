# Shared infrastructure for sim-study test scripts.
# All sim files are model-tagged so SPDE and BYM2 worlds can coexist:
#   simPopsSurveys_<MODEL>.RData
#   intPtsDHS_simStudy1_<MODEL>_<i>_K<u>_<r>.RData
# Every entry point accepts `model` ("spde" | "bym2"), `nsim`, and `regenerate`.

# ---------- filename helpers ----------
.modelTag    <- function(model) toupper(match.arg(model, c("spde", "bym2")))
.simPopsFile <- function(model)
    sprintf("savedOutput/simStudy1/simPopsSurveys_%s.RData", .modelTag(model))
.intPtsFile  <- function(simIdx, model, KDHSu = 16, KDHSr = 21)
    sprintf("savedOutput/simStudy1/intPtsDHS_simStudy1_%s_%d_K%d_%d.RData",
            .modelTag(model), simIdx, KDHSu, KDHSr)

# Disable interactive browser() prompts inside simData1*.
.silenceBrowser <- function() {
    assign("browser", function(...) invisible(NULL), envir = globalenv())
}

# TRUE when outFile is missing OR is a checkpoint with seed=seed and fewer than
# `nsim` completed sims (i.e. simData1* has more work to do).
.checkpointNeedsMoreSims <- function(outFile, nsim, seed) {
    if(!file.exists(outFile)) return(TRUE)
    chk <- new.env()
    ok  <- tryCatch({ load(outFile, envir = chk); TRUE },
                    error = function(e) FALSE)
    if(!ok || !exists("surveysDHS", envir = chk)) return(TRUE)
    if(!exists(".seedUsed", envir = chk)) return(TRUE)         # pre-checkpoint file
    if(!isTRUE(chk$.seedUsed == seed))    return(TRUE)         # seed mismatch
    length(chk$surveysDHS) < nsim
}

# Canonicalize DHS covariate column names: "urb" -> "urban", "pop" -> "normPop".
canonDhsNames <- function(m) {
    cn <- colnames(m)
    if(!is.null(cn)) {
        cn <- gsub("^urb$", "urban",   cn)
        cn <- gsub("^pop$", "normPop", cn)
        colnames(m) <- cn
    }
    m
}

# ---------- 1. simulate populations + surveys ----------
# Returns the names loaded (surveysDHS, surveysMICS, areaPops, subareaPops, ...),
# injected directly into the caller's environment for downstream convenience.
simulateSurveys <- function(model, nsim = 1, seed = 123, regenerate = FALSE,
                            envir = parent.frame()) {
    model   <- match.arg(model, c("spde", "bym2"))
    outFile <- .simPopsFile(model)

    # simData1* now writes directly to the model-tagged outFile and checkpoints
    # after every iteration, so a partial file can be resumed transparently.
    # Decide whether we need to (re)invoke the simulator at all.
    needsRun <- regenerate || .checkpointNeedsMoreSims(outFile, nsim, seed)
    if(!needsRun) {
        cat("[sim] loading cached ", outFile, " (>=", nsim, " sims complete)\n", sep = "")
    } else {
        # On regenerate=TRUE we must wipe the existing checkpoint, otherwise
        # simData1* sees the cached "all nsim complete" file and returns without
        # doing anything — so the new fields (e.g. *_smoothRisk) never get
        # written. Rename to .bak so it can be recovered if needed.
        if(regenerate && file.exists(outFile)) {
            bak <- paste0(outFile, ".bak")
            cat("[sim] regenerate=TRUE: moving existing ", outFile, " -> ", bak, "\n", sep = "")
            file.rename(outFile, bak)
        }
        cat("[sim] generating/resuming ", outFile,
            " (nsim=", nsim, ", seed=", seed, ")\n", sep = "")
        .silenceBrowser()
        t0 <- proc.time()[3]
        if(model == "bym2") simData1BYM2(nsim = nsim, seed = seed)
        else                simData1    (nsim = nsim, seed = seed)
        cat(sprintf("[sim] done in %.1f min\n", (proc.time()[3] - t0)/60))
    }
    loaded <- load(outFile, envir = envir)
    invisible(loaded)
}

# ---------- 2. DHS integration points ----------
buildSimDhsIntPts <- function(simIdx, surveysDHS, model,
                              KDHSu = 16, KDHSr = 21,
                              JInnerUrban = 4, JInnerRural = 4, JOuterRural = 1,
                              regenerate = FALSE) {
    model   <- match.arg(model, c("spde", "bym2"))
    outFile <- .intPtsFile(simIdx, model, KDHSu, KDHSr)

    if(!regenerate && file.exists(outFile)) {
        cat("[intPts sim ", simIdx, "] loading cached ", outFile, "\n", sep = "")
        load(outFile)                       # -> intPtsDHS
        return(intPtsDHS)
    }
    cat("[intPts sim ", simIdx, "] building ", outFile, " ...\n", sep = "")
    datDHS <- surveysDHS[[simIdx]]
    if(!("n" %in% names(datDHS)) && "N" %in% names(datDHS)) datDHS$n <- datDHS$N
    if(!("y" %in% names(datDHS)) && "Z" %in% names(datDHS)) datDHS$y <- datDHS$Z

    t0 <- proc.time()[3]
    intPtsDHS <- makeAllIntegrationPointsDHS(
        cbind(datDHS$east, datDHS$north), datDHS$urban,
        areaNames = datDHS$subarea, popPrior = TRUE,
        numPointsUrban = KDHSu, numPointsRural = KDHSr,
        JInnerUrban = JInnerUrban, JInnerRural = JInnerRural,
        JOuterRural = JOuterRural,
        adminMap = adm2Full, outFile = outFile, saveOutput = FALSE)
    cat(sprintf("[intPts sim %d] done in %.1f min\n", simIdx, (proc.time()[3] - t0)/60))
    save(intPtsDHS, file = outFile)
    invisible(intPtsDHS)
}

# ---------- 3. build inputsMDM for one sim ----------
buildInputsForSim <- function(simIdx, surveysDHS, surveysMICS, intPtsMICS, model,
                              KMICS = 100, KDHSu = 16, KDHSr = 21) {
    model <- match.arg(model, c("spde", "bym2"))
    load(.intPtsFile(simIdx, model, KDHSu, KDHSr))   # -> intPtsDHS
    datDHS  <- surveysDHS [[simIdx]]
    datMICS <- surveysMICS[[simIdx]]
    for(p in list(c("N","n"), c("N","ns"), c("Z","y"), c("Z","ys"))) {
        if(!(p[2] %in% names(datDHS)))  datDHS [[p[2]]] <- datDHS [[p[1]]]
        if(!(p[2] %in% names(datMICS))) datMICS[[p[2]]] <- datMICS[[p[1]]]
    }
    if(!("Stratum" %in% names(datMICS)))
        datMICS$Stratum <- adm2ToStratumMICS(datMICS$subarea)

    inp <- makeInputsMDM(datDHS = datDHS, datMICS = datMICS,
                         intPtsDHS = intPtsDHS, intPtsMICS = intPtsMICS,
                         KMICS = KMICS, KDHSurb = KDHSu, KDHSrur = KDHSr,
                         saveNewIntPts = FALSE)
    inp$intPtsDHS$covsUrb <- canonDhsNames(inp$intPtsDHS$covsUrb)
    inp$intPtsDHS$covsRur <- canonDhsNames(inp$intPtsDHS$covsRur)
    list(datDHS = datDHS, datMICS = datMICS, inputsMDM = inp)
}

# ---------- 4. truth values per generative model ----------
# SPDE has no φ; we record sqrt(margVar) as the closest analogue for σ_spatial.
modelTruths <- function(model) {
    model <- match.arg(model, c("spde", "bym2"))
    if(model == "bym2")
        list(sigmaSpatial = sqrt(0.5), phi = 0.8,
             sigmaEps = sqrt(1.5), alpha = -1.25)
    else
        list(sigmaSpatial = sqrt(0.5), phi = NA_real_,
             sigmaEps = sqrt(1.5), alpha = -1.25)
}

# ---------- 5. pretty-print fit summary ----------
printFitSummary <- function(res, truths) {
    p   <- res$opt$par
    SD  <- res$TMBsd
    nll <- if(!is.null(res$opt$value)) res$opt$value else res$opt$objective
    cat(sprintf("convergence=%s  pdHess=%s  NLL=%.4f\n",
                as.character(res$opt$convergence),
                if(inherits(SD,"sdreport")) as.character(SD$pdHess) else "n/a",
                nll))
    cat(sprintf("  sigmaSpatial = %.4f  (truth %.4f)\n",
                exp(-0.5*p["log_tau"]), truths$sigmaSpatial))
    phiTruth <- if(is.na(truths$phi)) "NA (SPDE)" else sprintf("%.4f", truths$phi)
    cat(sprintf("  phi          = %.6f  (truth %s)   [logit_phi = %.3f]\n",
                plogis(p["logit_phi"]), phiTruth, p["logit_phi"]))
    cat(sprintf("  sigmaEps     = %.4f  (truth %.4f)\n",
                exp(-0.5*p["log_tauEps"]), truths$sigmaEps))
    if("alpha" %in% names(p))
        cat(sprintf("  alpha        = %.4f  (truth %.4f)\n",
                    p["alpha"], truths$alpha))
}

# NLL as a function of phi, holding other params at their MLE.
nllProfileInPhi <- function(res, phiVals = c(0.50, 0.70, 0.80, 0.90, 0.95, 0.99, 0.999)) {
    obj <- res$TMBobj
    mle <- res$opt$par
    nlls <- sapply(phiVals, function(p) {
        par2 <- mle; par2["logit_phi"] <- qlogis(p); obj$fn(par2)
    })
    data.frame(phi = phiVals,
               logit_phi = round(qlogis(phiVals), 3),
               NLL       = round(nlls, 4),
               diffNLL   = round(nlls - min(nlls), 4))
}

# ---------- 5b. simulate + int pts for multiple models ----------
# Runs each model's prepareSims sequentially by default. With parallel=TRUE,
# spawns one child R process per model via callr — useful only when the host
# machine has enough RAM/CPU headroom (running 2+ heavy sims concurrently is
# what previously crashed the workstation).
prepareSimsForModels <- function(models = c("bym2", "spde"), nsim = 1, seed = 123,
                                 KDHSu = 16, KDHSr = 21,
                                 regenerate = FALSE, parallel = FALSE) {
    models <- sapply(models, function(m) match.arg(m, c("spde", "bym2")))
    if(!parallel) {
        for(m in models)
            prepareSims(model = m, nsim = nsim, seed = seed,
                        KDHSu = KDHSu, KDHSr = KDHSr, regenerate = regenerate)
        return(invisible(NULL))
    }
    if(!requireNamespace("callr", quietly = TRUE))
        stop("parallel=TRUE requires the `callr` package")
    cat("[prepareSimsForModels] launching ", length(models),
        " parallel R children\n", sep = "")
    # The parent's getwd() after source("setup.R") is the jittering ROOT
    # (setup.R setwd's there). Each child Rscript launches in code/ — where
    # setup.R lives — so the source("setup.R") line works as in interactive
    # use, and setup.R itself does the setwd to the project root.
    parentDir <- getwd()
    codeDir   <- file.path(parentDir, "code")
    if(!file.exists(file.path(codeDir, "setup.R")))
        stop("Could not locate setup.R from getwd()=", parentDir,
             "; run this from the jittering project root (where code/ lives).")
    workers <- lapply(models, function(m) {
        logPath <- file.path(parentDir, sprintf("code/prepareSims_%s.log",
                                                toupper(m)))
        callr::r_bg(
            func = function(model, nsim, seed, KDHSu, KDHSr, regenerate, codeDir) {
                setwd(codeDir)
                source("setup.R")
                source("code/simData.R")
                source("code/makeIntegrationPoints.R")
                source("code/testInfrastructure.R")
                prepareSims(model = model, nsim = nsim, seed = seed,
                            KDHSu = KDHSu, KDHSr = KDHSr, regenerate = regenerate)
            },
            args = list(model = m, nsim = nsim, seed = seed,
                        KDHSu = KDHSu, KDHSr = KDHSr,
                        regenerate = regenerate,
                        codeDir = codeDir),
            stdout = logPath,
            stderr = "2>&1"
        )
    })
    for(w in workers) w$wait()
    invisible(NULL)
}

# ---------- 6. simulate + all DHS int pts, no fitting ----------
prepareSims <- function(model, nsim = 1, seed = 123,
                        KDHSu = 16, KDHSr = 21, regenerate = FALSE) {
    model <- match.arg(model, c("spde", "bym2"))
    env <- new.env(parent = environment())
    simulateSurveys(model, nsim = nsim, seed = seed,
                    regenerate = regenerate, envir = env)
    surveysDHS <- env$surveysDHS
    stopifnot(length(surveysDHS) >= nsim)
    for(i in seq_len(nsim))
        buildSimDhsIntPts(i, surveysDHS, model,
                          KDHSu = KDHSu, KDHSr = KDHSr,
                          regenerate = regenerate)
    invisible(NULL)
}

# ---------- 7. top-level fitMM driver ----------
# Returns a list of fit objects (one per element of simIdx).
testFitMM <- function(model, nsim = 1, simIdx = seq_len(nsim), seed = 123,
                      KMICS = 100, KDHSu = 16, KDHSr = 21, Qgh = 10,
                      regenerate = FALSE, profilePhi = TRUE) {
    model <- match.arg(model, c("spde", "bym2"))
    stopifnot(all(simIdx >= 1), all(simIdx <= nsim))

    env <- new.env(parent = environment())
    simulateSurveys(model, nsim = nsim, seed = seed,
                    regenerate = regenerate, envir = env)
    surveysDHS  <- env$surveysDHS
    surveysMICS <- env$surveysMICS
    stopifnot(length(surveysDHS) >= nsim, length(surveysMICS) >= nsim)

    for(i in simIdx) buildSimDhsIntPts(i, surveysDHS, model,
                                       KDHSu = KDHSu, KDHSr = KDHSr,
                                       regenerate = regenerate)

    micsEnv <- new.env()
    load("savedOutput/global/intPtsMICS_100.RData", envir = micsEnv)
    intPtsMICS <- micsEnv$intPtsMICS

    truths <- modelTruths(model)
    results <- vector("list", length(simIdx))
    names(results) <- paste0("sim", simIdx)

    for(j in seq_along(simIdx)) {
        i <- simIdx[j]
        cat(sprintf("\n=== fitMM(%s sim %d, PC prior on phi) ===\n",
                    toupper(model), i))
        ip <- buildInputsForSim(i, surveysDHS, surveysMICS, intPtsMICS, model,
                                KMICS = KMICS, KDHSu = KDHSu, KDHSr = KDHSr)
        t0 <- proc.time()[3]
        res <- fitMM(datDHS = ip$datDHS, datMICS = ip$datMICS,
                     inputsMDM = ip$inputsMDM,
                     KMICS = KMICS, Qgh = Qgh, getSDs = TRUE, verbose = FALSE)
        cat(sprintf("fit time: %.1f min\n", (proc.time()[3] - t0)/60))
        printFitSummary(res, truths)
        if(profilePhi) {
            cat(sprintf("\n=== NLL profile in phi (%s sim %d) ===\n",
                        toupper(model), i))
            print(nllProfileInPhi(res), row.names = FALSE)
        }
        results[[j]] <- res
    }
    invisible(results)
}
