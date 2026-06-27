# SLURM driver for FE+nug+GH cross-validation over the 7 VALMODELS, both
# validation types. Thin wrapper: all logic lives in validationModels.R
# (runValFoldFE, runValArealFE) so the local tests and the cluster array run the
# SAME code.
#
#   srun --partition=CPUQ --time=02:00:00 --mem-per-cpu=10000 --pty bash
#
# Array index layout (decodeValJob over .buildValJobList, models = VALMODELS):
#   cluster jobs : 20-fold CV (10 DHS + 10 MICS held-out folds), collapsed for
#                  single-dataset models (the unused survey -> 1 "predict-all" job).
#                  Md/M_D/M_M = 11 jobs each, M_dm/M_DM = 20 each  -> 73 cluster.
#   areal jobs   : one per model (leave-out folds 6:10)            -> 5 areal.
# Total array size = 78.  (See decodeValJob for the exact ordering.)

source("setup.R")
options(error=traceback)
source("code/modFED.R")
source("code/modFEM.R")
source("code/modFEMD.R")
source("code/validation.R")          # getScores / direct-est deps
source("code/getDirectEsts.R")       # areal direct estimates
source("code/validationCommon.R")    # buildFullMDMInputs, subsetMDMInputs, predFEclusters, ...
source("code/validationModels.R")    # VALMODELS, decodeValJob, runValJobCluster, runValArealFE

index = as.numeric(commandArgs(trailingOnly=TRUE))
job   = decodeValJob(index)
cat("[runValidationFE] index =", index, " type =", job$type, " model =", job$model,
    " predict =", if(is.null(job$predict)) NA else job$predict,
    " fold =", job$fold, " tag =", job$tag, "\n")

if(!isTRUE(job$isValid)) {
    cat("[runValidationFE] invalid index — exiting cleanly.\n")
    quit(save="no", status=0)
}

# ---- canonical full-data inputs (cached on disk by prepValidationFE.R) ----------
fullCache = "savedOutput/validation/fullMDMInputs_FE.RData"
if(file.exists(fullCache)) {
    load(fullCache)                  # fullInp
} else {
    if(!dir.exists("savedOutput/validation")) dir.create("savedOutput/validation", recursive=TRUE)
    fullInp = buildFullMDMInputs(KMICS=100, KDHSurb=16, KDHSrur=21)
    save(fullInp, file=fullCache)
}

# ---- dispatch on validation type ------------------------------------------------
if(job$type == "cluster") {
    cat("[runValidationFE] CLUSTER job:", job$tag, "\n")
    runValJobCluster(job$model, job$predict, job$fold, fullInp,
                     outDir="savedOutput/validation/folds", nsim=10000)
} else {
    cat("[runValidationFE] AREAL:", job$tag, "(leave-out folds 6:10)\n")
    sc = runValArealFE(job$model, fullInp,
                       outDir="savedOutput/validation/areal",
                       leftOutFolds=6:10, nsim=10000)
    print(sc, row.names=FALSE)
}

cat("[runValidationFE] done. tag:", job$tag, "\n")
