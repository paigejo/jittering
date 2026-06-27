# Aggregate validation scores produced by the runValidationFE.sbatch array, over
# all VALMODELS and BOTH validation types. Standalone: safe to run after the
# array completes (or partially — missing files are skipped). Reuses the exact
# tag / file-name convention of runValFoldFE & runValArealFE (validationModels.R),
# so this only collates; the per-fold scores were already computed at run time.
#
# Output:
#   savedOutput/validation/scoreSummary_FE_cluster.csv   (cluster-level, per fold)
#   savedOutput/validation/scoreSummary_FE_areal.csv     (areal, per direct-est type)

source("setup.R")
options(error=traceback)
source("code/validation.R")           # getScores deps
source("code/getDirectEsts.R")
source("code/validationCommon.R")     # scoreTagFold, ...
source("code/validationModels.R")     # VALMODELS, valModelValid

foldDir  = "savedOutput/validation/folds"
arealDir = "savedOutput/validation/areal"

# ---- cluster-level: pool each model's predsJob_* files -> DHS/MICS/combined ----
# .scoreModelCluster gathers all of a model's per-job held-out predictions (DHS
# folds 1..10 + MICS folds 1..10, or the collapsed predict-all job) and scores
# per-survey + the denominator-weighted combined (weights = cluster size n).
clRows = list()
for(m in VALMODELS) {
    sc = .scoreModelCluster(m, foldDir = foldDir)
    if(!is.null(sc)) clRows[[length(clRows)+1]] = sc
}
if(length(clRows)) {
    cl = do.call(rbind, clRows)                      # already one row per (model, side)
    write.csv(cl, file="savedOutput/validation/scoreSummary_FE_cluster.csv", row.names=FALSE)
    cat("[scoreValidationFE] cluster summary:", nrow(cl), "rows ->",
        "savedOutput/validation/scoreSummary_FE_cluster.csv\n")
    metrics = intersect(c("CRPS","MSE","IntervalScore80","Coverage80","Width80"), names(cl))
    cl = cl[order(match(cl$side, c("combined","DHS","MICS")), match(cl$model, VALMODELS)), ]
    print(cl[, c("model","side", metrics)], row.names=FALSE, digits=4)   # combined = headline
} else {
    cat("[scoreValidationFE] no cluster score files found yet.\n")
}

# ---- areal: one row per (model, direct-est type) --------------------------------
arRows = list()
for(m in VALMODELS) {
    f = file.path(arealDir, paste0("areal_", m, "_areal_fold6-10.RData"))
    if(!file.exists(f)) next
    e = new.env(); load(f, envir=e)
    if(!is.null(e$scores)) arRows[[length(arRows)+1]] = e$scores
}
if(length(arRows)) {
    ar = do.call(rbind, arRows)
    write.csv(ar, file="savedOutput/validation/scoreSummary_FE_areal.csv", row.names=FALSE)
    cat("\n[scoreValidationFE] areal summary:", nrow(ar), "rows ->",
        "savedOutput/validation/scoreSummary_FE_areal.csv\n")
    print(ar[ar$directEst=="combined", c("model","Bias","MSE","CRPS","Coverage80")], row.names=FALSE)
} else {
    cat("\n[scoreValidationFE] no areal score files found yet.\n")
}
