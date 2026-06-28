# ============================================================================
# runValidationFull — callr-parallel FE validation (no SLURM), mirroring the
# simulation study's scoreSimStudyFull.R structure: platform-aware worker caps,
# serial TMB-template precompile, load-balanced phases via callr::r_bg child
# processes, exit/disk auditing, then collation.
#
# Two phases (each balanced across nWorkers):
#   Phase 1 - CLUSTER jobs: the 20-fold CV (collapsed to 11 for single-dataset
#             models) -> runValJobCluster, one task per (model, predictSurvey,
#             fold). 73 tasks.
#   Phase 2 - AREAL jobs: per-area leave-half-out -> runValArealOneArea, one task
#             per (model, area). 5 models x ~37 areas tasks.
#   Phase 3 - collate -> per-model DHS/MICS/combined cluster + areal score CSVs.
#
# Worker caps (from sessionInfo()$running): local (mac/Windows) = 4, else = 16.
# ============================================================================

.detectValWorkerCaps <- function(family = "FE") {
  running <- if(is.null(sessionInfo()$running)) "" else sessionInfo()$running
  isLocal <- grepl("macOS|mac OS|Windows|darwin", running, ignore.case = TRUE)
  # BYM2 fits are heavier (spatial field + INLA-fallback draws); fewer local workers.
  n <- if(isLocal) (if(family == "BYM2") 2L else 4L) else 16L
  list(n = n, label = paste(if(isLocal) "local:" else "cluster:", running))
}

.precompileValTemplates <- function(family = "FE") {
  templates <- if(family == "BYM2")
    c("modM_BYM2_GH_v2", "modD_BYM2_GH_v2", "modMDM_BYM2_GH_v2",   # fitMM / fitMD / fitMDM
      "modM_DSep", "modM_DSepRepar")                               # fitMd (observed coords)
  else
    c("modM_FEnug_GH", "modD_BYM2_GH_v2", "modMDM_BYM2_GH_v2")     # fitFEM / fitFED / fitFEMD
  for(t in templates) {
    if(!file.exists(paste0("code/", t, ".cpp"))) next
    if(any(file.exists(paste0("code/", t, c(".o", ".so", ".dll"))))) next
    cat(sprintf("[precompile] %s ...\n", t)); flush.console()
    TMB::compile(paste0("code/", t, ".cpp"), framework = "TMBad", safebounds = FALSE)
  }
  cat("[precompile] templates ready\n")
}

# Schedule `tasks` across `nWorkers` callr children with LPT (longest-processing-
# time) load balancing on each task's $est. Each worker sources the validation
# stack once, loads fullInp, then runs its task chunk (dispatching by $kind).
.runValPhase <- function(phaseName, tasks, nWorkers, repoRoot, logsDir,
                         nsim, arealNsim, leftOutFolds, foldDir, arealDir) {
  if(length(tasks) == 0) { cat(sprintf("[%s] no tasks.\n", phaseName)); return(invisible()) }
  if(!requireNamespace("callr", quietly = TRUE)) stop("runValidationFull needs the `callr` package")
  nW <- min(nWorkers, length(tasks))
  tasks <- tasks[order(-sapply(tasks, `[[`, "est"))]                 # LPT: heaviest first
  chunks <- vector("list", nW); load <- numeric(nW)
  for(t in tasks) { w <- which.min(load); chunks[[w]] <- c(chunks[[w]], list(t)); load[w] <- load[w] + t$est }
  cat(sprintf("[%s] %d tasks across %d workers (est wall %.1f min)\n",
              phaseName, length(tasks), nW, max(load)/60)); flush.console()

  workers <- lapply(seq_len(nW), function(w) {
    callr::r_bg(
      func = function(taskChunk, repoRoot, nsim, arealNsim, leftOutFolds, foldDir, arealDir) {
        setwd(repoRoot)
        suppressMessages({
          source("code/setup.R")
          source("code/modFED.R"); source("code/modFEM.R"); source("code/modFEMD.R")
          # BYM2 fitters (used when family=="BYM2")
          source("code/modMdSep.R"); source("code/modM_DSep.R")
          source("code/modM_DMSep.R"); source("code/modM_MSep.R")
          source("code/validation.R"); source("code/getDirectEsts.R")
          source("code/validationCommon.R"); source("code/validationModels.R")
        })
        load("savedOutput/validation/fullMDMInputs_FE.RData")          # fullInp
        # Per-task error isolation: a single failing fit must NOT kill the worker
        # and lose the rest of its chunk (that silently dropped ~30% of areas).
        for(t in taskChunk) {
          res <- tryCatch({
            if(t$kind == "cluster")
              runValJobCluster(t$model, t$predict, t$fold, fullInp, outDir = foldDir,
                               nsim = nsim, family = t$family)
            else
              runValArealOneArea(t$model, t$areaIdx, t$area, fullInp, outDir = arealDir,
                                 leftOutFolds = leftOutFolds, nsim = arealNsim, family = t$family)
            "ok"
          }, error = function(e) paste("ERR:", conditionMessage(e)))
          cat(sprintf("    [%s] %-26s %s\n", t$kind, t$tag, res))
        }
        TRUE
      },
      args = list(taskChunk = chunks[[w]], repoRoot = repoRoot, nsim = nsim,
                  arealNsim = arealNsim, leftOutFolds = leftOutFolds,
                  foldDir = foldDir, arealDir = arealDir),
      stdout = file.path(logsDir, sprintf("valFull_%s_w%d.log", phaseName, w)), stderr = "2>&1")
  })
  for(w in workers) w$wait()
  bad <- sum(sapply(workers, function(w) { ec <- tryCatch(w$get_exit_status(), error=function(e) NA); is.na(ec) || ec != 0 }))
  cat(sprintf("[%s] phase done; %d worker(s) exited non-zero%s\n", phaseName, bad,
              if(bad) paste0(" -- check ", logsDir, "/valFull_", phaseName, "_w*.log") else ""))
}

runValidationFull <- function(family = c("FE", "BYM2"),
                              nsim = 10000, arealNsim = 2000, leftOutFolds = 6:10,
                              repoRoot = getwd(), outDir = NULL,
                              KMICS = 100, KDHSu = 16, KDHSr = 21) {
  family <- match.arg(family)
  source("code/setup.R")
  source("code/modFED.R"); source("code/modFEM.R"); source("code/modFEMD.R")
  source("code/modMdSep.R"); source("code/modM_DSep.R"); source("code/modM_DMSep.R"); source("code/modM_MSep.R")
  source("code/validation.R"); source("code/getDirectEsts.R")
  source("code/validationCommon.R"); source("code/validationModels.R")

  # family-suffixed output dirs so FE and BYM2 runs don't clobber each other
  if(is.null(outDir)) outDir <- file.path("savedOutput/validation", tolower(family))
  foldDir  <- file.path(outDir, "folds")
  arealDir <- file.path(outDir, "areal")
  logsDir  <- file.path(outDir, "logs")
  for(d in c(foldDir, arealDir, logsDir)) if(!dir.exists(d)) dir.create(d, recursive = TRUE)

  caps <- .detectValWorkerCaps(family)
  cat(sprintf("\n=============== runValidationFull family=%s (%s) ===============\n", family, caps$label))
  cat(sprintf("workers = %d   nsim = %d   arealNsim = %d   outDir = %s\n", caps$n, nsim, arealNsim, outDir))
  .precompileValTemplates(family)

  fullCache <- "savedOutput/validation/fullMDMInputs_FE.RData"
  if(file.exists(fullCache)) load(fullCache) else {
    cat("[prep] building", fullCache, "...\n"); fullInp <- buildFullMDMInputs(KMICS=KMICS, KDHSurb=KDHSu, KDHSrur=KDHSr)
    save(fullInp, file = fullCache)
  }
  de <- arealDirectEsts(fullInp, leftOutFolds = leftOutFolds)
  areas <- de$combined$area
  t0 <- proc.time()[3]

  # ---- Phase 1: CLUSTER (both families; BYM2 adds the per-cluster field term) ----
  estCluster <- if(family=="BYM2") c(Md=60, M_D=90, M_dm=120, M_DM=240, M_M=110) else c(Md=12, M_D=15, M_dm=25, M_DM=75, M_M=45)
  clTasks <- list()
  for(j in .buildValJobList()) { if(j$type != "cluster") next
    clTasks[[length(clTasks)+1]] <- list(kind="cluster", model=j$model, predict=j$predict, fold=j$fold,
      family=family, est=unname(estCluster[j$model]),
      tag=paste0(j$model, "_pred", j$predict, if(is.na(j$fold)) "_all" else paste0("_fold", j$fold))) }
  cat(sprintf("\n--- Phase 1: CLUSTER (%d tasks, family=%s) ---\n", length(clTasks), family))
  .runValPhase("cluster", clTasks, caps$n, repoRoot, logsDir, nsim, arealNsim, leftOutFolds, foldDir, arealDir)

  # ---- Phase 2: AREAL (both families; predGrid handles the BYM2 field, useInla='auto') ----
  estAreal <- if(family=="BYM2") c(Md=60, M_D=80, M_dm=90, M_DM=180, M_M=110) else c(Md=12, M_D=18, M_dm=22, M_DM=85, M_M=50)
  arTasks <- list()
  for(m in VALMODELS) for(ai in seq_along(areas))
    arTasks[[length(arTasks)+1]] <- list(kind="areal", model=m, areaIdx=ai, area=areas[ai],
      family=family, est=unname(estAreal[m]), tag=sprintf("%s_area%03d", m, ai))
  cat(sprintf("\n--- Phase 2: AREAL (%d tasks, family=%s) ---\n", length(arTasks), family))
  .runValPhase("areal", arTasks, caps$n, repoRoot, logsDir, nsim, arealNsim, leftOutFolds, foldDir, arealDir)

  # ---- Phase 3: collate ----
  cat("\n--- Phase 3: collate ---\n")
  clRows <- list(); arRows <- list()
  for(m in VALMODELS) {
    sc <- .scoreModelCluster(m, foldDir = foldDir); if(!is.null(sc)) clRows[[length(clRows)+1]] <- sc
    ar <- scoreModelAreal(m, fullInp, arealDir = arealDir, leftOutFolds = leftOutFolds, de = de)
    if(!is.null(ar)) arRows[[length(arRows)+1]] <- ar
  }
  if(length(clRows)) {
    cl <- do.call(rbind, clRows); write.csv(cl, file.path(outDir, "scoreSummary_cluster.csv"), row.names=FALSE)
    cat("\n[cluster scores]\n"); print(cl[order(match(cl$side,c("combined","DHS","MICS")), match(cl$model,VALMODELS)),
        c("model","side","CRPS","MSE","Coverage80")], row.names=FALSE, digits=4)
  }
  if(length(arRows)) {
    ar <- do.call(rbind, arRows); write.csv(ar, file.path(outDir, "scoreSummary_areal.csv"), row.names=FALSE)
    cat("\n[areal scores]\n"); print(ar[order(match(ar$directEst,c("combined","DHS","MICS")), match(ar$model,VALMODELS)),
        c("model","directEst","nAreas","CRPS","MSE","Coverage80")], row.names=FALSE, digits=4)
  }
  cat(sprintf("\nTotal wall: %.1f min. CSVs in %s/.\n", (proc.time()[3]-t0)/60, outDir))
  invisible(NULL)
}
# This file only DEFINES functions (like scoreSimStudyFull.R). Launch via the
# driver:  Rscript code/runValFull.R   (mirrors runScoresFull_*.R for the sim study)
