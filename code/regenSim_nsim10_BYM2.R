# Regenerate simStudy1 from the BYM2 generative model (nsim = 10) and build
# K=16/21 DHS integration points for sims 1-10. Also writes the BYM2-simulated
# survey data to the canonical `simPopsSurveys.RData` path so existing test
# scripts pick it up without modification (the _BYM2 suffixed file is kept too).
#
# DGP (simData1BYM2 defaults):
#   alpha = -1.25,  gamma_urban = 1.0,  betaRest = (0,0,0,0.5)
#   sigmaEpsilon = sqrt(1.5)
#   sigmaBYM2    = sqrt(0.5) ~ 0.707
#   phi          = 0.8
#
# Outputs (savedOutput/simStudy1/):
#   simPopsSurveys_BYM2.RData    (written by simData1BYM2)
#   simPopsSurveys.RData         (copy of the above, canonical path)
#   intPtsDHS_simStudy1_<i>_K16_21.RData  for i = 1..10
#
# Wall time ~ 45 min total (sim ~15 min, int pts ~30 min).

source("setup.R")
options(error=traceback)
source("code/simData.R")
source("code/makeIntegrationPoints.R")

# simData1BYM2 has an interactive browser() near its save call; override.
browser <- function(...) invisible(NULL)

NSIM <- 10

# ---- Step 1: simulate NSIM BYM2 populations + surveys ----------------------
cat("=== Step 1: simData1BYM2(nsim=", NSIM, ") -> simPopsSurveys_BYM2.RData ===\n", sep="")
t0 <- proc.time()[3]
simData1BYM2(nsim = NSIM, seed = 123)
cat(sprintf("Step 1 done in %.1f min\n", (proc.time()[3] - t0)/60))

# Mirror the file to the canonical path used by downstream test scripts.
src <- "savedOutput/simStudy1/simPopsSurveys_BYM2.RData"
dst <- "savedOutput/simStudy1/simPopsSurveys.RData"
file.copy(src, dst, overwrite = TRUE)
cat("Copied", src, "->", dst, "\n")

# ---- Step 2: build K=16/21 DHS integration points for all NSIM sims --------
cat("\n=== Step 2: K=16/21 DHS integration points for sims 1..", NSIM, " ===\n", sep="")
load(dst)   # populates surveysDHS, surveysMICS, areaPops, subareaPops, ...
stopifnot(length(surveysDHS) >= NSIM)

KDHSu       <- 16
KDHSr       <- 21
JInnerUrban <- 4
JInnerRural <- 4
JOuterRural <- 1

step2Start <- proc.time()[3]
for(i in 1:NSIM) {
    outFile <- sprintf("savedOutput/simStudy1/intPtsDHS_simStudy1_%d_K16_21.RData", i)
    cat("[sim ", i, "/", NSIM, "] building ", outFile, " ...\n", sep="")
    datDHS <- surveysDHS[[i]]
    if(!("n" %in% names(datDHS)) && "N" %in% names(datDHS)) datDHS$n <- datDHS$N
    if(!("y" %in% names(datDHS)) && "Z" %in% names(datDHS)) datDHS$y <- datDHS$Z

    t0 <- proc.time()[3]
    intPtsDHS <- makeAllIntegrationPointsDHS(
        cbind(datDHS$east, datDHS$north), datDHS$urban,
        areaNames=datDHS$subarea, popPrior=TRUE,
        numPointsUrban=KDHSu, numPointsRural=KDHSr,
        JInnerUrban=JInnerUrban, JInnerRural=JInnerRural,
        JOuterRural=JOuterRural, adminMap=adm2Full,
        outFile=outFile, saveOutput=FALSE)
    cat(sprintf("[sim %d] built in %.1f min - wUrban=%d x %d, wRural=%d x %d\n",
                i, (proc.time()[3]-t0)/60,
                nrow(intPtsDHS$wUrban), ncol(intPtsDHS$wUrban),
                nrow(intPtsDHS$wRural), ncol(intPtsDHS$wRural)))
    save(intPtsDHS, file=outFile)
}
cat(sprintf("Step 2 done in %.1f min\n", (proc.time()[3] - step2Start)/60))

# ---- Step 3: invalidate stale inputsMDM and per-sim score caches -----------
cat("\n=== Step 3: invalidate stale caches ===\n")
staleInputs <- Sys.glob("savedOutput/simStudy1/inputs_*KMICS100_KDHS16_21.RData")
for(f in staleInputs) {
    cat("Removing stale inputsMDM cache: ", f, "\n", sep="")
    file.remove(f)
}
staleScores <- Sys.glob("savedOutput/simStudy1/scores/scores_*_sim*.RData")
for(f in staleScores) {
    cat("Removing stale score file: ", f, "\n", sep="")
    file.remove(f)
}

cat("\nDone.\n")
cat("Now rerun runSimStudyScores.R against the new BYM2-simulated data.\n")
