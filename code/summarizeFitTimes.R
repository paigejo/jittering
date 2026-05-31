# Pull fitTime / scoreTime from any score files that have them (post-patch),
# aggregate per model.
setwd("c:/Users/jpaige/git/jittering")
MODELS <- c("Md_FE","Md","M_D_FE","M_M_FE","M_DM_FE","M_D_BYM2","M_M_BYM2","M_DM_BYM2")

oneTag <- function(tag) {
    cat("\n========== [", tag, "] timing (minutes) ==========\n", sep="")
    rows <- list()
    for(m in MODELS) {
        fts <- c(); sts <- c()
        for(i in 1:10) {
            f <- sprintf("savedOutput/simStudy1/scores/%s/scores_%s_sim%d.RData", tag, m, i)
            if(!file.exists(f)) next
            e <- new.env(); load(f, envir = e)
            if(exists("fitTime",   envir = e)) fts <- c(fts, e$fitTime)
            if(exists("scoreTime", envir = e)) sts <- c(sts, e$scoreTime)
        }
        if(length(fts) == 0 && length(sts) == 0) next
        rows[[length(rows) + 1]] <- data.frame(
            model     = m,
            nSim      = length(fts),
            fitMean   = if(length(fts)) round(mean(fts),   2) else NA,
            fitMax    = if(length(fts)) round(max(fts),    2) else NA,
            scoreMean = if(length(sts)) round(mean(sts),   2) else NA,
            scoreMax  = if(length(sts)) round(max(sts),    2) else NA,
            totMean   = if(length(fts) && length(sts)) round(mean(fts) + mean(sts), 2) else NA)
    }
    if(length(rows)) print(do.call(rbind, rows), row.names = FALSE)
    else             cat("(no files with timing fields yet)\n")
}
oneTag("BYM2"); oneTag("SPDE")
