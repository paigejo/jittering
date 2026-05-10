# SLURM driver for FE+nug+GH cross-validation.
# srun --partition=CPUQ --time=02:00:00 --mem-per-cpu=10000 --pty bash
#
# Index decoding (i = model, j = fold; uniform maxJ = 20):
#   i: 1..4 corresponding to M_D, M_M, M_DM, M_DM_comm
#   j: 1..20  ->  j <= 10 is a DHS held-out fold, j > 10 is a MICS held-out fold
# Models that are DHS-only (M_D) skip MICS-fold jobs; MICS-only (M_M) skips DHS folds.

source("setup.R")
options(error=traceback)
source("code/modFEMD_comm.R")
source("code/modFED.R")
source("code/modFEMD.R")
source("code/validation.R")          # for stratKmics(), stratKdhs(), getScores deps, scoreValidationPreds
source("code/validationCommon.R")

index = as.numeric(commandArgs(trailingOnly=TRUE))

job = decodeFEjob(index, maxJ=20)
cat("[runValidationFE] index =", index, "\n")
cat("                   model =", job$model, "  fold =", job$fold,
    "  side =", job$side, "  tag =", job$tag, "\n")

if(!job$isValid) {
    cat("[runValidationFE] (model, side) combination invalid — exiting cleanly.\n")
    quit(save="no", status=0)
}

# ---- Build canonical full-data inputs (cached on disk) -------------------------
fullCache = "savedOutput/validation/fullMDMInputs_FE.RData"
if(file.exists(fullCache)) {
    load(fullCache)
} else {
    if(!dir.exists("savedOutput/validation")) dir.create("savedOutput/validation", recursive=TRUE)
    fullInp = buildFullMDMInputs(KMICS=100, KDHSurb=16, KDHSrur=21)
    save(fullInp, file=fullCache)
}

# ---- Build keep-vectors for this fold ------------------------------------------
# Folds were assigned inside buildFullMDMInputs and travel as columns on the data.
foldDHS_vec  = fullInp$datDHS$fold
foldMICS_vec = fullInp$datMICS$fold

if(job$side == "DHS") {
    dhsKeep  = foldDHS_vec != job$foldDHS
    micsKeep = rep(TRUE, length(foldMICS_vec))
    dhsHeldOut  = !dhsKeep
    micsHeldOut = rep(FALSE, length(foldMICS_vec))
} else {  # MICS fold
    dhsKeep  = rep(TRUE, length(foldDHS_vec))
    micsKeep = foldMICS_vec != job$foldMICS
    dhsHeldOut  = rep(FALSE, length(foldDHS_vec))
    micsHeldOut = !micsKeep
}

inSampleInp  = subsetMDMInputs(fullInp, dhsKeep,  micsKeep)
heldOutInp   = subsetMDMInputs(fullInp, dhsHeldOut, micsHeldOut)

cat("[runValidationFE] in-sample sizes: nDHSurb=", inSampleInp$nDHSurb,
    " nDHSrur=", inSampleInp$nDHSrur,
    " nMICSurb=", inSampleInp$nMICSurb, " nMICSrur=", inSampleInp$nMICSrur, "\n", sep="")
cat("[runValidationFE] held-out sizes : nDHSurb=", heldOutInp$nDHSurb,
    " nDHSrur=", heldOutInp$nDHSrur,
    " nMICSurb=", heldOutInp$nMICSurb, " nMICSrur=", heldOutInp$nMICSrur, "\n", sep="")

# ---- Fit the in-sample model ---------------------------------------------------
fitDir = "savedOutput/validation/folds"
if(!dir.exists(fitDir)) dir.create(fitDir, recursive=TRUE)
fitFile = file.path(fitDir, paste0("fit", job$tag, ".RData"))

cat("[runValidationFE] fitting", job$model, "on in-sample fold subset ...\n")
startTime = proc.time()[3]

if(job$model == "M_D") {
    fitRes = fitFED(datDHS=fullInp$datDHS, inputsMDM=inSampleInp,
                    KDHSurb=fullInp$KDHSurb, KDHSrur=fullInp$KDHSrur,
                    Qgh=10, getSDs=TRUE, verbose=FALSE)
} else if(job$model == "M_M") {
    fitRes = fitFEM(datDHS=fullInp$datDHS, datMICS=fullInp$datMICS,
                    inputsMDM=inSampleInp,
                    KMICS=fullInp$KMICS,
                    Qgh=10, getSDs=TRUE, verbose=FALSE)
} else if(job$model == "M_DM") {
    fitRes = fitFEMD(datDHS=fullInp$datDHS, datMICS=fullInp$datMICS,
                     inputsMDM=inSampleInp,
                     KMICS=fullInp$KMICS,
                     KDHSurb=fullInp$KDHSurb, KDHSrur=fullInp$KDHSrur,
                     Qgh=10, getSDs=TRUE, verbose=FALSE)
} else if(job$model == "M_DM_comm") {
    fitRes = fitFEMD_comm(datDHS=fullInp$datDHS, datMICS=fullInp$datMICS,
                          inputsMDM=inSampleInp,
                          KMICS=fullInp$KMICS,
                          KDHSurb=fullInp$KDHSurb, KDHSrur=fullInp$KDHSrur,
                          Qgh=10, getSDs=TRUE, verbose=FALSE)
}

totalTime = proc.time()[3] - startTime
cat("[runValidationFE] fit done in", round(totalTime/60, 2), "min ; saving fit ...\n")
save(fitRes, totalTime, file=fitFile)

# ---- Predict held-out clusters --------------------------------------------------
# Build the per-side out-of-sample piece needed by predFEclusters.
if(job$side == "DHS") {
    nUrb = heldOutInp$nDHSurb;  nRur = heldOutInp$nDHSrur
    outOfSample = list(
        covsUrb = heldOutInp$intPtsDHS$covsUrb,
        covsRur = heldOutInp$intPtsDHS$covsRur,
        wUrban  = heldOutInp$intPtsDHS$wUrban,
        wRural  = heldOutInp$intPtsDHS$wRural
    )
    ysUrb = heldOutInp$ysUrbDHS; nsUrb = heldOutInp$nsUrbDHS
    ysRur = heldOutInp$ysRurDHS; nsRur = heldOutInp$nsRurDHS
} else {  # MICS held-out
    nUrb = heldOutInp$nMICSurb; nRur = heldOutInp$nMICSrur
    outOfSample = list(
        covsUrb = heldOutInp$intPtsMICS$XUrb,
        covsRur = heldOutInp$intPtsMICS$XRur,
        wUrban  = heldOutInp$intPtsMICS$wUrban,
        wRural  = heldOutInp$intPtsMICS$wRural
    )
    ysUrb = heldOutInp$ysUrbMICS; nsUrb = heldOutInp$nsUrbMICS
    ysRur = heldOutInp$ysRurMICS; nsRur = heldOutInp$nsRurMICS
}

cat("[runValidationFE] predicting", nUrb, "urban /", nRur, "rural held-out clusters ...\n")
predObj = predFEclusters(fitRes, outOfSample, nUrb, nRur, nsim=10000, side="DHS")

formatAndSavePreds(predObj, ysUrb, nsUrb, ysRur, nsRur,
                   fold=job$fold, tag=job$tag, totalTime=totalTime,
                   outDir=fitDir)

cat("[runValidationFE] done. Saved fit and preds tagged:", job$tag, "\n")
