# ============================================================================
# INLA-style vs Gaussian (with Pt-step bug fixed) score comparison for
# M_M_BYM2 and M_D_BYM2 on all 10 BYM2 + 10 SPDE sims. One fit per
# (sim, model, generative); both methods scored from that single fit so the
# only thing varying across the comparison is the posterior-draws routine.
#
# Outputs: savedOutput/simStudy1/scores_inla_compare/<TAG>/scores_<mod>_sim<i>.RData
# containing scoresFE_gauss, scoresHyper_gauss, scoresArea_gauss and
# the analogous _inla blocks, plus timing.
# ============================================================================

# ---- driver-side imports ----------------------------------------------------
# (Workers re-source the same files inside their callr child.)

INLA_COMPARE_MODELS <- c("M_M_BYM2", "M_D_BYM2")

# Score one (model, sim, generative) fit with both Gauss and INLA-style.
.scoreOneFitInlaCompare <- function(simIdx, modName, outFile,
                                    surveysDHS, surveysMICS, intPtsMICS,
                                    model, truths, areaPops,
                                    KMICS, KDHSu, KDHSr, Qgh, NDRAWS, COVS) {
    inp <- buildInputsForSim(simIdx, surveysDHS, surveysMICS, intPtsMICS, model,
                             KMICS = KMICS, KDHSu = KDHSu, KDHSr = KDHSr)

    # Fit once
    t0 <- proc.time()[3]
    res <- tryCatch(.fitOne(modName, inp, KMICS, KDHSu, KDHSr, Qgh, COVS),
                    error = function(e) {
                        cat("    [", modName, " sim ", simIdx, "] FIT ERROR: ",
                            conditionMessage(e), "\n", sep = "")
                        NULL
                    })
    if(is.null(res)) return(invisible(FALSE))
    fitTime <- (proc.time()[3] - t0)/60

    # ---- Gaussian (corrected Pt-step) scoring -------------------------------
    t1 <- proc.time()[3]
    scoresFE_gauss <- tryCatch(.scoreFE(res, NDRAWS),
        error = function(e) { cat("    scoreFE(g) err: ", conditionMessage(e), "\n"); NULL })
    scoresHyper_gauss <- tryCatch(.scoreHyper(res, truths, NDRAWS),
        error = function(e) { cat("    scoreHyper(g) err: ", conditionMessage(e), "\n"); NULL })
    scoresArea_gauss <- tryCatch({
        # predGrid useInla=FALSE
        grid <- predGrid(res$TMBsd, popMat = popMatNGAThresh, nsim = NDRAWS,
                         obj = res$TMBobj, admLevel = "stratMICS",
                         useInla = FALSE, res = res)
        agg <- predArea(grid, areaVarName = "area", orderedAreas = adm1@data$NAME_1)
        truth <- areaPops[, simIdx]
        estMat <- as.matrix(agg$aggregationResults$p)
        rownames(estMat) <- agg$aggregationResults$region
        estMat <- estMat[match(names(truth), rownames(estMat)), , drop = FALSE]
        if(ncol(estMat) == 1) {
            getScores(truth=truth, est=estMat[,1],
                      significance=c(.5,.8,.9,.95), doFuzzyReject=TRUE,
                      getAverage=TRUE, na.rm=TRUE)
        } else {
            getScores(truth=truth, estMat=estMat,
                      significance=c(.5,.8,.9,.95), doFuzzyReject=TRUE,
                      getAverage=TRUE, na.rm=TRUE)
        }
    }, error = function(e) { cat("    scoreArea(g) err: ", conditionMessage(e), "\n"); NULL })
    scoreTime_gauss <- (proc.time()[3] - t1)/60

    # ---- INLA-style scoring -------------------------------------------------
    t2 <- proc.time()[3]
    inla_ok <- TRUE
    inla_draws <- tryCatch(inlaStyleDraws(res, NDRAWS = NDRAWS,
                                          deltaZ = 1.0, deltaPi = 2.5,
                                          maxAxialSteps = 4),
                           error = function(e) {
                               cat("    inlaStyleDraws err: ", conditionMessage(e), "\n")
                               inla_ok <<- FALSE
                               NULL
                           })

    scoresFE_inla <- NULL; scoresHyper_inla <- NULL; scoresArea_inla <- NULL
    if(inla_ok && !is.null(inla_draws)) {
        # Hand-roll FE / Hyper scoring directly from the INLA draws matrix
        # (mirrors .scoreFE / .scoreHyper logic but uses the supplied draws).
        rn  <- rownames(inla_draws)
        aRow <- inla_draws[rn == "alpha", , drop = FALSE]
        bRows <- inla_draws[rn == "beta", , drop = FALSE]
        feDraws <- rbind(aRow, bRows)
        if(nrow(feDraws) == length(TRUE_FE)) {
            scoresFE_inla <- tryCatch(getScores(truth = TRUE_FE, estMat = feDraws,
                                                significance = c(.5,.8,.9,.95),
                                                doFuzzyReject = TRUE,
                                                getAverage = FALSE, na.rm = TRUE),
                                      error = function(e) { cat("    scoreFE(i) err: ", conditionMessage(e), "\n"); NULL })
        }
        # Hyper scoring on natural scale
        spec <- list(
            sigmaBYM2 = list(par="log_tau",    f=function(x) exp(-0.5*x), truth=truths$sigmaSpatial),
            phi       = list(par="logit_phi",  f=plogis,                  truth=truths$phi),
            sigmaEps  = list(par="log_tauEps", f=function(x) exp(-0.5*x), truth=truths$sigmaEps))
        rows <- list()
        for(name in names(spec)) {
            s <- spec[[name]]
            idx <- which(rn == s$par)
            if(length(idx) == 0) next
            nat <- s$f(inla_draws[idx, ])
            if(is.na(s$truth)) {
                sc <- getScores(truth = mean(nat), estMat = matrix(nat, nrow = 1),
                                significance = c(.5,.8,.9,.95), doFuzzyReject = TRUE,
                                getAverage = FALSE, na.rm = TRUE)
                sc[, c("Bias","Var","MSE","RMSE","CRPS")] <- NA
                sc[, grep("^(IntervalScore|Coverage)", colnames(sc))] <- NA
            } else {
                sc <- getScores(truth = s$truth, estMat = matrix(nat, nrow = 1),
                                significance = c(.5,.8,.9,.95), doFuzzyReject = TRUE,
                                getAverage = FALSE, na.rm = TRUE)
            }
            sc <- cbind(est = mean(nat), sd = sd(nat), truth = s$truth, sc)
            rownames(sc) <- name
            rows[[name]] <- sc
        }
        if(length(rows) > 0) scoresHyper_inla <- do.call(rbind, rows)

        # Area-level via predGrid(useInla=TRUE)
        scoresArea_inla <- tryCatch({
            grid <- predGrid(res$TMBsd, popMat = popMatNGAThresh, nsim = NDRAWS,
                             obj = res$TMBobj, admLevel = "stratMICS",
                             useInla = TRUE, res = res)
            agg <- predArea(grid, areaVarName = "area", orderedAreas = adm1@data$NAME_1)
            truth <- areaPops[, simIdx]
            estMat <- as.matrix(agg$aggregationResults$p)
            rownames(estMat) <- agg$aggregationResults$region
            estMat <- estMat[match(names(truth), rownames(estMat)), , drop = FALSE]
            if(ncol(estMat) == 1) {
                getScores(truth=truth, est=estMat[,1],
                          significance=c(.5,.8,.9,.95), doFuzzyReject=TRUE,
                          getAverage=TRUE, na.rm=TRUE)
            } else {
                getScores(truth=truth, estMat=estMat,
                          significance=c(.5,.8,.9,.95), doFuzzyReject=TRUE,
                          getAverage=TRUE, na.rm=TRUE)
            }
        }, error = function(e) { cat("    scoreArea(i) err: ", conditionMessage(e), "\n"); NULL })
    }
    scoreTime_inla <- (proc.time()[3] - t2)/60

    save(scoresFE_gauss, scoresHyper_gauss, scoresArea_gauss,
         scoresFE_inla,  scoresHyper_inla,  scoresArea_inla,
         fitTime, scoreTime_gauss, scoreTime_inla,
         file = outFile)
    cat(sprintf("    [%s sim %d] %s  fit=%.1f gauss=%.1f inla=%.1f min\n",
                model, simIdx, modName, fitTime, scoreTime_gauss, scoreTime_inla))
    invisible(TRUE)
}

scoreInlaCompare <- function(model, nsim = 10, nWorkers = 4,
                             KMICS = 100, KDHSu = 16, KDHSr = 21,
                             Qgh = 10, NDRAWS = 1000,
                             COVS = c("urban","access","elev","distRiversLakes","normPop"),
                             regenerate = FALSE) {
    if(!requireNamespace("callr", quietly = TRUE))
        stop("scoreInlaCompare needs callr")
    model <- match.arg(model, c("spde","bym2"))
    tag   <- toupper(model)
    outDir <- sprintf("savedOutput/simStudy1/scores_inla_compare/%s", tag)
    if(!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)

    # LPT scheduling (M_M_BYM2 longer than M_D_BYM2; same per-model cost)
    estMin <- c(M_M_BYM2 = 18, M_D_BYM2 = 7)
    tasks <- list()
    for(simIdx in seq_len(nsim)) {
        for(modName in INLA_COMPARE_MODELS) {
            outFile <- sprintf("%s/scores_%s_sim%d.RData", outDir, modName, simIdx)
            if(!regenerate && file.exists(outFile)) next
            tasks[[length(tasks)+1]] <- list(simIdx = simIdx, modName = modName,
                                             outFile = outFile,
                                             est = unname(estMin[modName]))
        }
    }
    if(length(tasks) == 0) {
        cat(sprintf("[%s] all files already present, skipping.\n", tag))
        return(invisible(NULL))
    }

    nW     <- min(nWorkers, length(tasks))
    tasks  <- tasks[order(-sapply(tasks, `[[`, "est"))]
    chunks <- vector("list", nW); loads <- numeric(nW)
    for(t in tasks) {
        w <- which.min(loads)
        chunks[[w]] <- c(chunks[[w]], list(t))
        loads[w] <- loads[w] + t$est
    }
    cat(sprintf("[%s] %d tasks across %d workers; est walltime %.0f min\n",
                tag, length(tasks), nW, max(loads)))

    parentDir <- getwd(); codeDir <- file.path(parentDir, "code")
    workers <- lapply(seq_len(nW), function(w) {
        logPath <- file.path(codeDir,
                             sprintf("inlaCompare_%s_w%d.log", tag, w))
        callr::r_bg(
            func = function(taskChunk, model, KMICS, KDHSu, KDHSr,
                            Qgh, NDRAWS, COVS, codeDir) {
                setwd(codeDir)
                source("setup.R")
                source("code/modFED.R");  source("code/modFEM.R")
                source("code/modFEMD.R"); source("code/modM_DSep.R")
                source("code/modM_MSep.R"); source("code/modM_DMSep.R")
                source("code/modMdSep.R"); source("code/makeInputsTMB.R")
                source("code/modBYM2.R")
                source("code/testInfrastructure.R")
                source("code/scoreSimStudy.R")
                source("code/inlaStyleDraws.R")
                source("code/scoreInlaCompare.R")

                simEnv <- new.env()
                simulateSurveys(model, nsim = max(sapply(taskChunk, `[[`, "simIdx")),
                                regenerate = FALSE, envir = simEnv)
                micsEnv <- new.env()
                load("savedOutput/global/intPtsMICS_100.RData", envir = micsEnv)
                truths <- modelTruths(model)
                for(t in taskChunk) {
                    .scoreOneFitInlaCompare(t$simIdx, t$modName, t$outFile,
                                            simEnv$surveysDHS, simEnv$surveysMICS,
                                            micsEnv$intPtsMICS, model, truths,
                                            simEnv$areaPops,
                                            KMICS, KDHSu, KDHSr, Qgh, NDRAWS, COVS)
                }
            },
            args = list(taskChunk = chunks[[w]], model = model,
                        KMICS = KMICS, KDHSu = KDHSu, KDHSr = KDHSr,
                        Qgh = Qgh, NDRAWS = NDRAWS, COVS = COVS,
                        codeDir = codeDir),
            stdout = logPath, stderr = "2>&1"
        )
    })
    for(w in workers) w$wait()
    cat(sprintf("[%s] all workers done.\n", tag))
}
