# Regenerate simStudy1 (nsim = 2) + K=16/21 DHS integration points for sims 1 & 2.
# Then invalidate the stale `inputs_SPDE_sim1_KMICS100_KDHS16_21.RData` cache
# so the next test_sim_SPDE_FE_then_BYM2.R run actually rebuilds from the new
# survey data + new integration points.
#
# Outputs:
#   savedOutput/simStudy1/simPopsSurveys.RData               (overwritten, nsim=2)
#   savedOutput/simStudy1/intPtsDHS_simStudy1_1_K16_21.RData (overwritten)
#   savedOutput/simStudy1/intPtsDHS_simStudy1_2_K16_21.RData (new)
#
# Wall time: ~15-25 min total.

source("setup.R")
options(error=traceback)
source("code/simData.R")
source("code/makeIntegrationPoints.R")

# simData1 has an interactive browser() call right before save (line ~168).
# Override it so this can run non-interactively.
browser <- function(...) invisible(NULL)

# ---- Step 1: simulate 2 SPDE populations + surveys ----
cat("=== Step 1: simData1(nsim=2)  ->  simPopsSurveys.RData ===\n")
t0 <- proc.time()[3]
simData1(nsim=2, seed=123)
cat(sprintf("Step 1 done in %.1f min\n", (proc.time()[3] - t0)/60))

# ---- Step 2: build K=16/21 DHS integration points for both sims ----
cat("\n=== Step 2: K=16/21 DHS integration points for sims 1, 2 ===\n")
load("savedOutput/simStudy1/simPopsSurveys.RData")
stopifnot(length(surveysDHS) >= 2)
KDHSu       <- 16
KDHSr       <- 21
JInnerUrban <- 4
JInnerRural <- 4
JOuterRural <- 1

for(i in 1:2) {
    outFile <- sprintf("savedOutput/simStudy1/intPtsDHS_simStudy1_%d_K16_21.RData", i)
    cat("[sim ", i, "] building ", outFile, " ...\n", sep="")
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

# ---- Step 3: invalidate stale inputsMDM caches (so the test rebuilds) ----
staleCaches <- c(
    "savedOutput/simStudy1/inputs_SPDE_sim1_KMICS100_KDHS16_21.RData",
    "savedOutput/simStudy1/inputs_SPDE_sim2_KMICS100_KDHS16_21.RData",
    "savedOutput/simStudy1/inputs_sim1_KMICS100_KDHS16_21.RData"
)
for(f in staleCaches) {
    if(file.exists(f)) {
        cat("Removing stale inputsMDM cache: ", f, "\n", sep="")
        file.remove(f)
    }
}

cat("\nDone.\n")
cat("Now you can rerun test_sim_SPDE_FE_then_BYM2.R locally to confirm K=16/21 on fresh data,\n")
cat("then copy simPopsSurveys.RData + the two intPtsDHS_simStudy1_*_K16_21.RData files to the cluster.\n")
