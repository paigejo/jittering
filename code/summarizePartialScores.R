# Aggregate whatever score files exist so far in BYM2/ and SPDE/ subdirs.
# Read-only — safe to run while the scoring processes are still writing.

MODELS <- c("Md_FE",  "Md",
            "M_D_FE", "M_M_FE",  "M_DM_FE",
            "M_D_BYM2", "M_M_BYM2", "M_DM_BYM2")

TRUE_FE_NAMES <- c("alpha", "beta_urban", "beta_access",
                   "beta_elev", "beta_distRiversLakes", "beta_normPop")

loadModelScores <- function(modName, outDir, nsim = 10) {
    out <- list(FE = list(), Hyper = list(), Area = list(), Subarea = list(),
                fitTimes = numeric(0), scoreTimes = numeric(0))
    for(simIdx in seq_len(nsim)) {
        f <- sprintf("%s/scores_%s_sim%d.RData", outDir, modName, simIdx)
        if(!file.exists(f)) next
        e <- new.env(); load(f, envir = e)
        out$FE     [[length(out$FE)      + 1]] <- e$scoresFE
        out$Hyper  [[length(out$Hyper)   + 1]] <- if(exists("scoresHyper", envir = e)) e$scoresHyper else NULL
        out$Area   [[length(out$Area)    + 1]] <- e$scoresArea
        out$Subarea[[length(out$Subarea) + 1]] <- e$scoresSubarea
        if(exists("fitTime",   envir = e)) out$fitTimes   <- c(out$fitTimes,   e$fitTime)
        if(exists("scoreTime", envir = e)) out$scoreTimes <- c(out$scoreTimes, e$scoreTime)
    }
    out
}
avgScoreList <- function(lst) {
    lst <- Filter(function(x) !is.null(x), lst)
    if(length(lst) == 0) return(NULL)
    arr <- simplify2array(lapply(lst, as.matrix))
    if(length(dim(arr)) == 2) arr else apply(arr, c(1, 2), mean, na.rm = TRUE)
}
buildAggTable <- function(modelToScoreMat) {
    rows <- list()
    for(m in names(modelToScoreMat)) {
        a <- modelToScoreMat[[m]]
        if(is.null(a)) next
        rows[[length(rows) + 1]] <- cbind(model = m, as.data.frame(a))
    }
    if(length(rows) == 0) NULL else do.call(rbind, rows)
}

countAvailable <- function(outDir, nsim = 10) {
    sapply(MODELS, function(m) {
        sum(sapply(seq_len(nsim),
                   function(i) file.exists(sprintf("%s/scores_%s_sim%d.RData", outDir, m, i))))
    })
}

summarizeOne <- function(model) {
    tag    <- toupper(model)
    outDir <- sprintf("savedOutput/simStudy1/scores/%s", tag)
    if(!dir.exists(outDir)) { cat("[", tag, "] dir does not exist yet\n", sep=""); return() }

    avail <- countAvailable(outDir)
    cat(sprintf("\n========== [%s] sims completed per model (out of 10) ==========\n", tag))
    print(avail)

    modelData <- setNames(lapply(MODELS, loadModelScores, outDir = outDir, nsim = 10), MODELS)

    # Area-level
    areaTab <- buildAggTable(lapply(modelData, function(d) avgScoreList(d$Area)))
    cat(sprintf("\n----- [%s] Area-level scores (mean across available sims) -----\n", tag))
    if(!is.null(areaTab)) print(areaTab, row.names = FALSE, digits = 4)

    # Subarea-level
    subTab <- buildAggTable(lapply(modelData, function(d) avgScoreList(d$Subarea)))
    cat(sprintf("\n----- [%s] Subarea-level scores (mean across available sims) -----\n", tag))
    if(!is.null(subTab)) print(subTab, row.names = FALSE, digits = 4)

    # Fixed-effect scores (one block per model)
    cat(sprintf("\n----- [%s] Fixed-effect scores (mean across available sims) -----\n", tag))
    for(m in MODELS) {
        a <- avgScoreList(modelData[[m]]$FE)
        if(is.null(a)) next
        if(nrow(a) == length(TRUE_FE_NAMES)) rownames(a) <- TRUE_FE_NAMES
        cat("\n[", m, "]\n", sep = "")
        print(round(a, 4))
    }

    # Hyperparameter scores (mean across available sims)
    cat(sprintf("\n----- [%s] Hyperparameter scores (mean across available sims) -----\n", tag))
    anyHyper <- FALSE
    for(m in MODELS) {
        a <- avgScoreList(modelData[[m]]$Hyper)
        if(is.null(a)) next
        anyHyper <- TRUE
        cat("\n[", m, "]\n", sep = "")
        print(round(a, 4))
    }
    if(!anyHyper) cat("(no scoresHyper fields present yet — they were added after these score files were written)\n")

    # Wall-time per model (mean across available sims)
    cat(sprintf("\n----- [%s] Wall-time per model (minutes; mean across available sims) -----\n", tag))
    timeRows <- list()
    for(m in MODELS) {
        ft <- modelData[[m]]$fitTimes
        st <- modelData[[m]]$scoreTimes
        if(length(ft) == 0 && length(st) == 0) next
        timeRows[[length(timeRows) + 1]] <- data.frame(
            model    = m,
            nSim     = length(ft),
            fitMin   = round(mean(ft),   2),
            scoreMin = round(mean(st),   2),
            totalMin = round(mean(ft) + mean(st), 2))
    }
    if(length(timeRows) > 0) print(do.call(rbind, timeRows), row.names = FALSE)
    else cat("(no timing fields present yet)\n")
}

setwd("c:/Users/jpaige/git/jittering")
summarizeOne("bym2")
summarizeOne("spde")
cat("\nDone.\n")
