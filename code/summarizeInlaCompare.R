# Aggregate Gaussian-vs-INLA score comparison.
setwd("c:/Users/jpaige/git/jittering")

MODELS    <- c("M_M_BYM2", "M_D_BYM2")
FE_NAMES  <- c("alpha","beta_urban","beta_access","beta_elev","beta_distRiversLakes","beta_normPop")
HYP_NAMES <- c("sigmaBYM2","phi","sigmaEps")

loadSet <- function(tag, modName) {
    out <- list(FE_g=list(), FE_i=list(), Hyp_g=list(), Hyp_i=list(),
                Area_g=list(), Area_i=list(),
                fit=numeric(), score_g=numeric(), score_i=numeric())
    dir <- sprintf("savedOutput/simStudy1/scores_inla_compare/%s", tag)
    for(i in 1:10) {
        f <- sprintf("%s/scores_%s_sim%d.RData", dir, modName, i)
        if(!file.exists(f)) next
        e <- new.env(); load(f, envir = e)
        out$FE_g    <- c(out$FE_g,    list(e$scoresFE_gauss))
        out$FE_i    <- c(out$FE_i,    list(e$scoresFE_inla))
        out$Hyp_g   <- c(out$Hyp_g,   list(e$scoresHyper_gauss))
        out$Hyp_i   <- c(out$Hyp_i,   list(e$scoresHyper_inla))
        out$Area_g  <- c(out$Area_g,  list(e$scoresArea_gauss))
        out$Area_i  <- c(out$Area_i,  list(e$scoresArea_inla))
        if(exists("fitTime",        envir=e)) out$fit     <- c(out$fit,     e$fitTime)
        if(exists("scoreTime_gauss",envir=e)) out$score_g <- c(out$score_g, e$scoreTime_gauss)
        if(exists("scoreTime_inla", envir=e)) out$score_i <- c(out$score_i, e$scoreTime_inla)
    }
    out
}
avgMat <- function(lst) {
    lst <- Filter(function(x) !is.null(x), lst)
    if(length(lst) == 0) return(NULL)
    arr <- simplify2array(lapply(lst, as.matrix))
    if(length(dim(arr)) == 2) arr else apply(arr, c(1,2), mean, na.rm = TRUE)
}

cmpRow <- function(label, gaussV, inlaV, fmt = "%.4f") {
    fmtNA <- function(v) if(is.na(v)||is.null(v)) "    NA " else sprintf(fmt, v)
    cat(sprintf("  %-26s  G=%s   I=%s\n", label, fmtNA(gaussV), fmtNA(inlaV)))
}

reportArea <- function(g, i) {
    cat("  metric                       Gaussian       INLA-style\n")
    keep <- c("Bias","Var","MSE","RMSE","CRPS",
              "IntervalScore50","IntervalScore80","IntervalScore90","IntervalScore95",
              "Coverage50","Coverage80","Coverage90","Coverage95",
              "Width50","Width80","Width90","Width95")
    g <- as.matrix(g); i <- as.matrix(i)
    for(k in keep) {
        if(!(k %in% colnames(g))) next
        cmpRow(k, g[1, k], i[1, k])
    }
}

reportFEorHyp <- function(label, g, i, rowNames = NULL) {
    cat(sprintf("\n  [%s]\n", label))
    cat("                     Bias_G   Bias_I    sd_G   sd_I    Cov95_G Cov95_I    Width95_G W95_I\n")
    g <- as.matrix(g); i <- as.matrix(i)
    if(!is.null(rowNames)) rownames(g) <- rownames(i) <- rowNames
    for(r in rownames(g)) {
        cat(sprintf("  %-19s  %+.3f   %+.3f    %.3f  %.3f    %.2f    %.2f       %.3f    %.3f\n",
                    r,
                    g[r,"Bias"],   i[r,"Bias"],
                    sqrt(if("Var" %in% colnames(g)) g[r,"Var"] else NA),
                    sqrt(if("Var" %in% colnames(i)) i[r,"Var"] else NA),
                    g[r,"Coverage95"], i[r,"Coverage95"],
                    g[r,"Width95"], i[r,"Width95"]))
    }
}

cat("\n========== BYM2-simulated, 10 sims, M_M_BYM2 ==========\n")
mm <- loadSet("BYM2","M_M_BYM2")
cat(sprintf("\n  Sims aggregated: FE_g=%d FE_i=%d  Area_g=%d Area_i=%d\n",
            length(mm$FE_g), length(mm$FE_i), length(mm$Area_g), length(mm$Area_i)))
cat(sprintf("  Mean walltime/sim: fit=%.1f min, scoring G=%.1f min, scoring I=%.1f min\n",
            mean(mm$fit), mean(mm$score_g), mean(mm$score_i)))

cat("\n--- Fixed-effect scores (mean across 10 sims) ---\n")
reportFEorHyp("FE (Gauss vs INLA)",
              avgMat(mm$FE_g), avgMat(mm$FE_i), rowNames = FE_NAMES)
cat("\n--- Hyperparameter scores (mean across 10 sims) ---\n")
reportFEorHyp("Hyper", avgMat(mm$Hyp_g), avgMat(mm$Hyp_i))
cat("\n--- Area-level prevalence scores (mean across 10 sims) ---\n")
reportArea(avgMat(mm$Area_g), avgMat(mm$Area_i))


cat("\n\n========== BYM2-simulated, 10 sims, M_D_BYM2 ==========\n")
md <- loadSet("BYM2","M_D_BYM2")

doSPDE <- function() {
  cat("\n\n========== SPDE-simulated, 10 sims, M_M_BYM2 ==========\n")
  mm <- loadSet("SPDE","M_M_BYM2")
  cat(sprintf("\n  Sims aggregated: FE_g=%d FE_i=%d  Area_g=%d Area_i=%d\n",
              length(mm$FE_g), length(mm$FE_i), length(mm$Area_g), length(mm$Area_i)))
  cat(sprintf("  Mean walltime/sim: fit=%.1f min, scoring G=%.1f min, scoring I=%.1f min\n",
              mean(mm$fit), mean(mm$score_g), mean(mm$score_i)))
  cat("\n--- Fixed-effect scores (mean across 10 sims) ---\n")
  reportFEorHyp("FE", avgMat(mm$FE_g), avgMat(mm$FE_i), rowNames = FE_NAMES)
  cat("\n--- Hyperparameter scores (mean across 10 sims) ---\n")
  reportFEorHyp("Hyper", avgMat(mm$Hyp_g), avgMat(mm$Hyp_i))
  cat("\n--- Area-level prevalence scores (mean across 10 sims) ---\n")
  reportArea(avgMat(mm$Area_g), avgMat(mm$Area_i))

  cat("\n\n========== SPDE-simulated, 10 sims, M_D_BYM2 ==========\n")
  md <- loadSet("SPDE","M_D_BYM2")
  cat(sprintf("\n  Sims aggregated: FE_g=%d FE_i=%d  Area_g=%d Area_i=%d\n",
              length(md$FE_g), length(md$FE_i), length(md$Area_g), length(md$Area_i)))
  cat(sprintf("  Mean walltime/sim: fit=%.1f min, scoring G=%.1f min, scoring I=%.1f min\n",
              mean(md$fit), mean(md$score_g), mean(md$score_i)))
  cat("\n--- Fixed-effect scores (mean across 10 sims) ---\n")
  reportFEorHyp("FE", avgMat(md$FE_g), avgMat(md$FE_i), rowNames = FE_NAMES)
  cat("\n--- Hyperparameter scores (mean across 10 sims) ---\n")
  reportFEorHyp("Hyper", avgMat(md$Hyp_g), avgMat(md$Hyp_i))
  cat("\n--- Area-level prevalence scores (mean across 10 sims) ---\n")
  reportArea(avgMat(md$Area_g), avgMat(md$Area_i))
}
cat(sprintf("\n  Sims aggregated: FE_g=%d FE_i=%d  Area_g=%d Area_i=%d\n",
            length(md$FE_g), length(md$FE_i), length(md$Area_g), length(md$Area_i)))
cat(sprintf("  Mean walltime/sim: fit=%.1f min, scoring G=%.1f min, scoring I=%.1f min\n",
            mean(md$fit), mean(md$score_g), mean(md$score_i)))

cat("\n--- Fixed-effect scores (mean across 10 sims) ---\n")
reportFEorHyp("FE", avgMat(md$FE_g), avgMat(md$FE_i), rowNames = FE_NAMES)
cat("\n--- Hyperparameter scores (mean across 10 sims) ---\n")
reportFEorHyp("Hyper", avgMat(md$Hyp_g), avgMat(md$Hyp_i))
cat("\n--- Area-level prevalence scores (mean across 10 sims) ---\n")
reportArea(avgMat(md$Area_g), avgMat(md$Area_i))

doSPDE()
cat("\nDone.\n")
