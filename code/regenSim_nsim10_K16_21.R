# Regenerate simStudy1 (nsim = 10 SPDE populations) + K=16/21 DHS integration
# points for sims 1-10. Then invalidate any stale inputs_SPDE_sim* caches.
#
# Outputs:
#   savedOutput/simStudy1/simPopsSurveys.RData              (overwritten, nsim=10)
#   savedOutput/simStudy1/intPtsDHS_simStudy1_<i>_K16_21.RData  for i = 1..10
#
# Wall time: ~45 min total (~15 min for sims, ~30 min for int pts).

source("setup.R")
options(error=traceback)
source("code/simData.R")
source("code/makeIntegrationPoints.R")

# simData1 has an interactive browser() near its save call; override so the
# script runs non-interactively.
browser <- function(...) invisible(NULL)

NSIM <- 10

# ---- Step 1: simulate NSIM SPDE populations + surveys ----
cat("=== Step 1: simData1(nsim=", NSIM, ")  ->  simPopsSurveys.RData ===\n", sep="")
t0 <- proc.time()[3]
simData1(nsim=NSIM, seed=123)
cat(sprintf("Step 1 done in %.1f min\n", (proc.time()[3] - t0)/60))

# ---- Step 2: build K=16/21 DHS integration points for all NSIM sims ----
cat("\n=== Step 2: K=16/21 DHS integration points for sims 1..", NSIM, " ===\n", sep="")
load("savedOutput/simStudy1/simPopsSurveys.RData")
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
    cat(sprintf("[sim %d] built in %.1f min — wUrban=%d x %d, wRural=%d x %d\n",
                i, (proc.time()[3]-t0)/60,
                nrow(intPtsDHS$wUrban), ncol(intPtsDHS$wUrban),
                nrow(intPtsDHS$wRural), ncol(intPtsDHS$wRural)))
    save(intPtsDHS, file=outFile)
}
cat(sprintf("Step 2 done in %.1f min\n", (proc.time()[3] - step2Start)/60))

# ---- Step 3: invalidate stale inputsMDM caches ----
cat("\n=== Step 3: invalidate stale inputs_SPDE_* caches ===\n")
staleCaches <- c(
    Sys.glob("savedOutput/simStudy1/inputs_SPDE_sim*_KMICS100_KDHS16_21.RData"),
    "savedOutput/simStudy1/inputs_sim1_KMICS100_KDHS16_21.RData"
)
for(f in unique(staleCaches)) {
    if(file.exists(f)) {
        cat("Removing stale cache: ", f, "\n", sep="")
        file.remove(f)
    }
}

cat("\nDone.\n")
cat("Now you can rerun test_sim_SPDE_FE_then_BYM2.R locally to confirm,\n")
cat("then copy simPopsSurveys.RData + the 10 intPtsDHS_simStudy1_<i>_K16_21.RData files to the cluster.\n")
