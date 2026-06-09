# ============================================================================
# scoreSimStudyFull — end-to-end scoring of all 8 sim-study models on N sims
# of either the BYM2 or SPDE generative process, with platform-aware
# parallelism and INLA-style auto-fallback for non-PSD jointPrecision.
#
# Pipeline:
#   Phase 1 — FE-only models (Md_FE, M_D_FE, M_M_FE, M_DM_FE) in parallel.
#   Phase 2 — BYM2 / Md models (Md, M_D_BYM2, M_M_BYM2, M_DM_BYM2) in parallel
#             with a stricter worker cap (BYM2 fits + INLA CCD walks peak at
#             >12 GB per process; over-parallelising crashes laptops).
#   Phase 3 — collate per-(sim, model) score files into a summary.
#
# Worker counts are chosen from sessionInfo()$running:
#   Mac / Windows  ->  8 workers for FE,  1 for BYM2 (serial; protects RAM)
#   Other (Linux)  -> 32 workers for FE,  8 for BYM2 (cluster nodes)
#
# Per-fit posterior draws default to useInla = "auto": Gaussian Cholesky of
# jointPrecision when it's PSD, falling back to INLA-style otherwise.
# ============================================================================

.detectWorkerCaps <- function() {
    inf <- sessionInfo()
    running <- if(is.null(inf$running)) "" else inf$running
    isLocal <- grepl("macOS|mac OS|Windows|darwin", running, ignore.case = TRUE)
    # nBYM2 is for the LIGHT BYM2 models (Md, M_D_BYM2, M_M_BYM2 — single data
    # source each, ~1-2 GB transient per fit). nMDM_BYM2 is for the HEAVY
    # combined M_DM_BYM2 fit which uses DHS + MICS + BYM2 + Laplace and peaks
    # at ~10 GB transient; we run it in a dedicated phase with far fewer
    # workers to avoid OOM-kills (every BYM2 worker that got SIGKILL'd in the
    # June-07 run died on an M_DM_BYM2 fit).
    if(isLocal) list(nFE = 8L,  nBYM2 = 1L, nMDM_BYM2 = 1L,
                     label = paste("local:",   running))
    else        list(nFE = 16L, nBYM2 = 12L, nMDM_BYM2 = 4L,
                     label = paste("cluster:", running))
}

# Compile any missing TMB templates serially up front. Prevents a race where
# multiple callr workers hit compile() on the same .cpp simultaneously and
# corrupt each other's .o/.so writes (a known cause of silent FE-worker death
# on the cluster).
.precompileTemplates <- function() {
    templates <- c("modM_BYM2_GH_v2", "modD_BYM2_GH_v2", "modMDM_BYM2_GH_v2",
                   "modM_FEnug_GH",   "modM_DSepRepar",  "modM_DSep")
    needed <- character(0)
    for(t in templates) {
        if(!any(file.exists(paste0("code/", t, c(".o", ".so", ".dll")))))
            needed <- c(needed, t)
    }
    if(length(needed) == 0) {
        cat("[precompile] all TMB templates already built\n")
        return(invisible(NULL))
    }
    cat(sprintf("[precompile] %d template(s) missing; building serially: %s\n",
                length(needed), paste(needed, collapse = ", ")))
    for(t in needed) {
        src <- paste0("code/", t, ".cpp")
        if(!file.exists(src)) {
            warning("[precompile] source missing: ", src)
            next
        }
        t0 <- proc.time()[3]
        TMB::compile(src, framework = "TMBad", safebounds = FALSE)
        cat(sprintf("[precompile] %s built in %.1f s\n", t, proc.time()[3]-t0))
    }
}

# Phase runner: schedule the given `models` across `nWorkers` callr children.
.runPhase <- function(phaseName, models, model, nsim, simIdxList,
                      outDir, nWorkers,
                      KMICS, KDHSu, KDHSr, Qgh, NDRAWS, COVS,
                      regenerate, useInla, truthQuantity) {
    if(!requireNamespace("callr", quietly = TRUE))
        stop("scoreSimStudyFull needs the `callr` package")

    # Per-model walltime estimates (minutes) — tuned from prior runs with
    # NDRAWS=1000. Used for LPT load balancing.
    estMin <- c(M_DM_BYM2 = 100, M_M_BYM2 = 20, M_D_BYM2 = 10, Md = 3,
                M_DM_FE = 6,  M_M_FE = 3,  M_D_FE = 3,  Md_FE = 3)

    tasks <- list()
    for(simIdx in simIdxList) {
        for(modName in models) {
            outFile <- sprintf("%s/scores_%s_sim%d.RData", outDir, modName, simIdx)
            if(!regenerate && file.exists(outFile)) next
            tasks[[length(tasks) + 1]] <- list(simIdx  = simIdx,
                                               modName = modName,
                                               outFile = outFile,
                                               est = unname(estMin[modName]))
        }
    }
    if(length(tasks) == 0) {
        cat(sprintf("[%s] all score files already present, skipping.\n", phaseName))
        return(invisible(NULL))
    }

    nW     <- min(nWorkers, length(tasks))
    tasks  <- tasks[order(-sapply(tasks, `[[`, "est"))]
    chunks <- vector("list", nW); loads <- numeric(nW)
    for(t in tasks) {
        w <- which.min(loads)
        chunks[[w]] <- c(chunks[[w]], list(t))
        loads[w]   <- loads[w] + t$est
    }
    cat(sprintf("[%s] %d tasks across %d workers; est walltime %.0f min\n",
                phaseName, length(tasks), nW, max(loads)))

    parentDir <- getwd(); codeDir <- file.path(parentDir, "code")
    # Logs go under scores_full/logs/ (shared across generative models, with
    # the generative tag in the filename) instead of dumping into code/.
    logsDir <- file.path(dirname(outDir), "logs")
    if(!dir.exists(logsDir)) dir.create(logsDir, recursive = TRUE)
    tag <- toupper(model)
    workers <- lapply(seq_len(nW), function(w) {
        logPath <- file.path(logsDir, sprintf("scoreFull_%s_%s_w%d.log",
                                              tag, phaseName, w))
        callr::r_bg(
            func = function(taskChunk, model, KMICS, KDHSu, KDHSr,
                            Qgh, NDRAWS, COVS, useInla, codeDir, truthQuantity) {
                setwd(codeDir)
                source("setup.R")
                source("code/modFED.R");   source("code/modFEM.R")
                source("code/modFEMD.R");  source("code/modM_DSep.R")
                source("code/modM_MSep.R");source("code/modM_DMSep.R")
                source("code/modMdSep.R"); source("code/makeInputsTMB.R")
                source("code/modBYM2.R")
                source("code/testInfrastructure.R")
                source("code/inlaStyleDraws.R")
                source("code/scoreSimStudy.R")

                simEnv <- new.env()
                simulateSurveys(model,
                                nsim = max(sapply(taskChunk, `[[`, "simIdx")),
                                regenerate = FALSE, envir = simEnv)
                micsEnv <- new.env()
                load("savedOutput/global/intPtsMICS_100.RData", envir = micsEnv)
                truths   <- modelTruths(model)
                areaPops <- .pickAreaTruth(simEnv, truthQuantity)
                for(t in taskChunk) {
                    .scoreOneFit(t$simIdx, t$modName, t$outFile,
                                 simEnv$surveysDHS, simEnv$surveysMICS,
                                 micsEnv$intPtsMICS, model, truths,
                                 areaPops,
                                 KMICS, KDHSu, KDHSr, Qgh, NDRAWS, COVS)
                }
            },
            args = list(taskChunk = chunks[[w]], model = model,
                        KMICS = KMICS, KDHSu = KDHSu, KDHSr = KDHSr,
                        Qgh = Qgh, NDRAWS = NDRAWS, COVS = COVS,
                        useInla = useInla, codeDir = codeDir,
                        truthQuantity = truthQuantity),
            stdout = logPath, stderr = "2>&1"
        )
    })
    for(w in workers) w$wait()

    # Audit worker exit statuses so silent crashes don't go unnoticed.
    badWorkers <- integer(0)
    for(wi in seq_along(workers)) {
        ec <- tryCatch(workers[[wi]]$get_exit_status(),
                       error = function(e) NA_integer_)
        if(is.na(ec) || ec != 0L) {
            badWorkers <- c(badWorkers, wi)
            cat(sprintf("[%s] worker %d exited with status %s\n",
                        phaseName, wi, as.character(ec)))
        }
    }
    # Audit on-disk completion vs scheduled count.
    nDone <- sum(sapply(tasks, function(t) file.exists(t$outFile)))
    if(nDone < length(tasks)) {
        cat(sprintf("[%s] WARNING: %d / %d scheduled tasks left no score file on disk\n",
                    phaseName, length(tasks) - nDone, length(tasks)))
        cat(sprintf("[%s] check %s/scoreFull_%s_%s_w*.log for errors.\n",
                    phaseName, logsDir, tag, phaseName))
    }
    cat(sprintf("[%s] phase done (%d / %d wrote files; %d worker(s) exited non-zero).\n",
                phaseName, nDone, length(tasks), length(badWorkers)))
}

# Collate per-(sim, model) score files into summary tables and save.
# Prints, in order:
#   * area-level prevalence scores (one row per model)
#   * fixed-effect / covariate scores (one block per model; rows = parameters)
#   * hyperparameter scores (one block per model that has them)
.collateScoresFull <- function(model, nsim, outDir, allModels) {
    if(!exists("MODELS")) {
        # scoreSimStudy.R defines helpers we re-use here
        source("code/scoreSimStudy.R")
    }
    tag <- toupper(model)
    modelData <- setNames(lapply(allModels, .loadModelScores,
                                 outDir = outDir, nsim = nsim), allModels)

    # ── Area-level (one row per model) ─────────────────────────────────────
    areaTab <- .buildAggTable(lapply(modelData,
                                     function(d) .avgScoreList(d$Area)))
    cat(sprintf("\n----- [%s] Area-level scores (mean across sims) -----\n", tag))
    if(!is.null(areaTab)) print(areaTab, row.names = FALSE, digits = 4)

    # ── Fixed-effect / covariate scores (one block per model) ──────────────
    cat(sprintf("\n----- [%s] Fixed-effect / covariate scores (mean across sims) -----\n", tag))
    feNames <- names(TRUE_FE)
    feTabs  <- list()                  # cache for the save below
    for(m in allModels) {
        a <- .avgScoreList(modelData[[m]]$FE)
        if(is.null(a)) next
        if(nrow(a) == length(feNames)) rownames(a) <- feNames
        cat(sprintf("\n[%s]  (n = %d sims)\n", m, length(modelData[[m]]$FE)))
        print(round(a, 4))
        feTabs[[m]] <- a
    }

    # ── Hyperparameter scores (only models that estimate them) ─────────────
    cat(sprintf("\n----- [%s] Hyperparameter scores (mean across sims) -----\n", tag))
    hypTabs <- list()
    anyHyper <- FALSE
    for(m in allModels) {
        a <- .avgScoreList(modelData[[m]]$Hyper)
        if(is.null(a)) next
        anyHyper <- TRUE
        cat(sprintf("\n[%s]\n", m))
        print(round(a, 4))
        hypTabs[[m]] <- a
    }
    if(!anyHyper) cat("(no models produced hyperparameter scores)\n")

    summaryFile <- file.path(outDir, sprintf("scoresSummary_%s.RData", tag))
    save(areaTab, feTabs, hypTabs, modelData, file = summaryFile)
    cat(sprintf("\n[%s] summary saved to %s\n", tag, summaryFile))
    invisible(list(areaTab = areaTab, feTabs = feTabs,
                   hypTabs = hypTabs, modelData = modelData))
}

# ============================================================================
# Main entry point. Call from a driver script.
# ============================================================================
scoreSimStudyFull <- function(model = c("bym2","spde"),
                              nsim = 100,
                              simIdxList = seq_len(nsim),
                              regenerate = FALSE,
                              useInla = "auto",
                              truthQuantity = c("smoothRisk","fineScalePrevalence"),
                              outBase = "savedOutput/simStudy1/scores_full",
                              KMICS = 100, KDHSu = 16, KDHSr = 21,
                              Qgh = 10, NDRAWS = 1000,
                              COVS = c("urban","access","elev","distRiversLakes","normPop")) {
    model         <- match.arg(model)
    truthQuantity <- match.arg(truthQuantity)
    tag           <- toupper(model)

    FE_MODELS        <- c("Md_FE", "M_D_FE", "M_M_FE", "M_DM_FE")
    # Split BYM2 into LIGHT (single data source — ~1-2 GB transient) and
    # HEAVY (combined DHS+MICS+BYM2 — ~10 GB transient). They go in separate
    # phases so the heavy one runs with a small nWorkers and doesn't OOM-kill
    # the lighter workers via shared-RAM pressure.
    BYM2_LIGHT       <- c("Md", "M_D_BYM2", "M_M_BYM2")
    BYM2_HEAVY       <- c("M_DM_BYM2")
    ALL_MODELS       <- c(FE_MODELS, BYM2_LIGHT, BYM2_HEAVY)

    outDir <- file.path(outBase, tag)
    if(!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)

    # If the caller asked for a full regeneration, wipe any pre-existing score
    # files for the (model, simIdx) cells we're about to score AND the per-
    # phase worker log files. Without this, OOM-killed workers from a previous
    # run would leave behind stale score files that .runPhase would skip
    # (regenerate is checked per-file, so a partial outFile from before could
    # survive even with regenerate=TRUE if the current fit fails) and ambiguous
    # log files mixing this run's output with a previous run's. Doing both
    # up-front guarantees a clean slate.
    if(isTRUE(regenerate)) {
        nRemoved <- 0
        for(m in c(FE_MODELS, BYM2_LIGHT, BYM2_HEAVY)) {
            for(simIdx in simIdxList) {
                f <- sprintf("%s/scores_%s_sim%d.RData", outDir, m, simIdx)
                if(file.exists(f)) {
                    file.remove(f); nRemoved <- nRemoved + 1
                }
            }
        }
        logsDir <- file.path(dirname(outDir), "logs")
        oldLogs <- character(0)
        if(dir.exists(logsDir)) {
            oldLogs <- list.files(
                logsDir,
                pattern = sprintf("^scoreFull_%s_(FE|BYM2|MDM_BYM2)_w[0-9]+\\.log$", tag),
                full.names = TRUE)
            if(length(oldLogs) > 0) file.remove(oldLogs)
        }
        cat(sprintf("[regenerate=TRUE] Removed %d score files from %s/ and %d worker logs from %s/\n",
                    nRemoved, outDir, length(oldLogs), logsDir))
    }

    caps <- .detectWorkerCaps()
    cat(sprintf("\n================ scoreSimStudyFull (%s) ================\n", tag))
    cat(sprintf("Platform: %s\n", caps$label))
    cat(sprintf("Worker caps: FE = %d, BYM2 light = %d, BYM2 heavy (M_DM) = %d\n",
                caps$nFE, caps$nBYM2, caps$nMDM_BYM2))

    # Precompile any missing TMB templates serially before any callr worker
    # launches, so parallel workers can dyn.load existing .so/.dll files
    # without racing on compile().
    .precompileTemplates()
    cat(sprintf("nsim = %d  (simIdx %d..%d)\n",
                length(simIdxList), min(simIdxList), max(simIdxList)))
    cat(sprintf("Output dir: %s\n", outDir))
    cat(sprintf("useInla = %s   regenerate = %s   truthQuantity = \"%s\"\n",
                deparse(useInla), as.character(regenerate), truthQuantity))
    t0 <- proc.time()[3]

    # Phase 1: FE-only models (light, parallelise widely)
    cat(sprintf("\n--- Phase 1: FE-only models (nWorkers = %d) ---\n", caps$nFE))
    .runPhase("FE", FE_MODELS, model, nsim, simIdxList, outDir, caps$nFE,
              KMICS, KDHSu, KDHSr, Qgh, NDRAWS, COVS, regenerate, useInla, truthQuantity)

    # Phase 2: LIGHT BYM2 models (single data source, modest memory)
    cat(sprintf("\n--- Phase 2: BYM2 light models %s (nWorkers = %d) ---\n",
                paste(BYM2_LIGHT, collapse=","), caps$nBYM2))
    .runPhase("BYM2", BYM2_LIGHT, model, nsim, simIdxList, outDir, caps$nBYM2,
              KMICS, KDHSu, KDHSr, Qgh, NDRAWS, COVS, regenerate, useInla, truthQuantity)

    # Phase 3: HEAVY BYM2 model (M_DM_BYM2 alone — ~10 GB per worker peak)
    cat(sprintf("\n--- Phase 3: BYM2 heavy M_DM_BYM2 (nWorkers = %d) ---\n",
                caps$nMDM_BYM2))
    .runPhase("MDM_BYM2", BYM2_HEAVY, model, nsim, simIdxList, outDir, caps$nMDM_BYM2,
              KMICS, KDHSu, KDHSr, Qgh, NDRAWS, COVS, regenerate, useInla, truthQuantity)

    # Phase 4: collate
    cat("\n--- Phase 4: collating per-sim scores ---\n")
    .collateScoresFull(model, nsim, outDir, ALL_MODELS)

    cat(sprintf("\nTotal walltime: %.1f hours\n", (proc.time()[3] - t0) / 3600))
    invisible(NULL)
}

# Convenience: score BYM2 then SPDE in one call.
scoreSimStudyFullBoth <- function(...) {
    scoreSimStudyFull(model = "bym2", ...)
    scoreSimStudyFull(model = "spde", ...)
}
