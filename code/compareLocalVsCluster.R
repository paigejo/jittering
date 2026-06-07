# Compare cluster vs laptop score files for the 10 sims that were run on both.
# - Cluster files (copied via scp into OneDrive):
#     c:/Users/jpaige/OneDrive - Norsk Regnesentral/Projects/jittering/jitterScores/scores_full_cluster/{BYM2,SPDE}/
# - Laptop files (from runScoresParallel):
#     c:/Users/jpaige/git/jittering/savedOutput/simStudy1/scores/{BYM2,SPDE}/
#
# For each (model, simIdx) present on both sides we compute element-wise
# max-abs differences on scoresFE, scoresHyper, scoresArea (and report the
# fit/score wallclock differences as a sanity check).  An exact match means
# the cluster run reused cached files; a small numerical difference is
# expected from different RNG seeds in posterior draws; a large difference
# means a real divergence between machines.

clustBase <- "c:/Users/jpaige/OneDrive - Norsk Regnesentral/Projects/jittering/jitterScores/scores_full_cluster"
localBase <- "c:/Users/jpaige/git/jittering/savedOutput/simStudy1/scores"

MODELS <- c("Md_FE","Md","M_D_FE","M_M_FE","M_DM_FE",
            "M_D_BYM2","M_M_BYM2","M_DM_BYM2")
SIMS   <- 1:10

loadScores <- function(path) {
    if(!file.exists(path)) return(NULL)
    e <- new.env(); load(path, envir = e)
    list(FE    = if(exists("scoresFE",    envir=e)) e$scoresFE    else NULL,
         Hyper = if(exists("scoresHyper", envir=e)) e$scoresHyper else NULL,
         Area  = if(exists("scoresArea",  envir=e)) e$scoresArea  else NULL,
         fitT  = if(exists("fitTime",     envir=e)) e$fitTime     else NA,
         scrT  = if(exists("scoreTime",   envir=e)) e$scoreTime   else NA)
}

maxAbsDiff <- function(a, b) {
    if(is.null(a) || is.null(b))               return(NA)
    a <- as.matrix(a); b <- as.matrix(b)
    if(!all(dim(a) == dim(b)))                 return(Inf)  # shape mismatch
    max(abs(a - b), na.rm = TRUE)
}

compareOne <- function(tag) {
    rows <- list()
    for(sim in SIMS) {
        for(m in MODELS) {
            cPath <- sprintf("%s/%s/scores_%s_sim%d.RData", clustBase, tag, m, sim)
            lPath <- sprintf("%s/%s/scores_%s_sim%d.RData", localBase, tag, m, sim)
            cHas <- file.exists(cPath);  lHas <- file.exists(lPath)
            if(!cHas && !lHas) next
            cS <- if(cHas) loadScores(cPath) else NULL
            lS <- if(lHas) loadScores(lPath) else NULL
            rows[[length(rows)+1]] <- data.frame(
                tag    = tag,
                sim    = sim,
                model  = m,
                cluster= cHas,
                local  = lHas,
                dFE    = maxAbsDiff(cS$FE,    lS$FE),
                dHyper = maxAbsDiff(cS$Hyper, lS$Hyper),
                dArea  = maxAbsDiff(cS$Area,  lS$Area),
                fit_c  = if(!is.null(cS)) round(cS$fitT, 2) else NA,
                fit_l  = if(!is.null(lS)) round(lS$fitT, 2) else NA,
                scr_c  = if(!is.null(cS)) round(cS$scrT, 2) else NA,
                scr_l  = if(!is.null(lS)) round(lS$scrT, 2) else NA,
                stringsAsFactors = FALSE
            )
        }
    }
    if(length(rows) == 0) return(NULL)
    do.call(rbind, rows)
}

cat("##### BYM2 generative #####\n")
bym2 <- compareOne("BYM2")
print(bym2, row.names = FALSE)

cat("\n##### SPDE generative #####\n")
spde <- compareOne("SPDE")
print(spde, row.names = FALSE)

# Verdict per model: do the FE / Hyper / Area numbers match exactly (cluster
# reused laptop files), are they close (different RNG, same fit), or far?
verdict <- function(d, label) {
    cat(sprintf("\n##### %s: per-model summary across the 10 sims #####\n", label))
    for(m in unique(d$model)) {
        sub <- d[d$model == m, ]
        both <- sub[sub$cluster & sub$local, ]
        if(nrow(both) == 0) {
            cat(sprintf("  [%-10s] no overlap (cluster=%d, local=%d of 10)\n",
                        m, sum(sub$cluster), sum(sub$local)))
            next
        }
        cat(sprintf("  [%-10s] overlap n=%2d   max|dFE|=%.3g   max|dHyper|=%.3g   max|dArea|=%.3g\n",
                    m, nrow(both),
                    max(both$dFE,    na.rm=TRUE),
                    max(both$dHyper, na.rm=TRUE),
                    max(both$dArea,  na.rm=TRUE)))
    }
}
verdict(bym2, "BYM2")
verdict(spde, "SPDE")
