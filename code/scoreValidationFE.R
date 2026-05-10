# Run the existing getScores() (from code/scores.R) over all per-fold preds files
# produced by runValidationFE.R, for every FE-family model. Standalone driver:
# safe to run after the SLURM array completes (or partially — it skips missing
# preds files).
#
# Output: savedOutput/validation/folds/scores<tag>_fold<N>.RData per fold, plus a
# combined CSV summary at savedOutput/validation/scoreSummary_FE.csv.

source("setup.R")
options(error=traceback)
source("code/validation.R")           # for getScores (via scores.R sourced from setup)
source("code/validationCommon.R")

models = c("M_D", "M_M", "M_DM", "M_DM_comm")

for(m in models) {
    cat("[scoreValidationFE] scoring", m, "...\n")
    scoreAllFoldsFE(m)
}

# Aggregate into a single summary table.
rows = list()
foldDir = "savedOutput/validation/folds"
for(m in models) {
    for(j in 1:20) {
        side = if(j <= 10) "DHS" else "MICS"
        tag = paste0(m, ifelse(side == "DHS", "_dhsFold", "_micsFold"),
                     ifelse(side == "DHS", j, j - 10))
        scoreFile = file.path(foldDir, paste0("scores", tag, "_fold", j, ".RData"))
        if(!file.exists(scoreFile)) next
        e = new.env(); load(scoreFile, envir=e)
        if(!is.null(e$scores)) {
            df = as.data.frame(e$scores)
            df$model = m
            df$side  = side
            df$fold  = j
            df$tag   = tag
            rows[[length(rows)+1]] = df
        }
    }
}

if(length(rows) > 0) {
    summary = do.call(rbind, rows)
    write.csv(summary, file="savedOutput/validation/scoreSummary_FE.csv", row.names=FALSE)
    cat("\nSummary written to savedOutput/validation/scoreSummary_FE.csv (",
        nrow(summary), " rows)\n")
} else {
    cat("\nNo score files found yet — re-run after the SLURM array completes.\n")
}
