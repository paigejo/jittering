# Regenerate DHS integration points for sims 2-10 of simStudy1 at K=16/21.
# Sim 1 already has intPtsDHS_simStudy1_1_K16_21.RData on disk; this batch
# fills in sims 2 through 10 so we can ship a 10-replicate set to the cluster.
#
# Output files: savedOutput/simStudy1/intPtsDHS_simStudy1_<i>_K16_21.RData
# Wall time: ~5 min/sim x 9 sims = ~45 min.

source("setup.R")
options(error=traceback)
source("code/makeIntegrationPoints.R")   # makeAllIntegrationPointsDHS

# Reuse settings from fit_GH_SPDE_FEnug.R and the existing _K16_21 build.
KDHSu       <- 16
KDHSr       <- 21
JInnerUrban <- 4
JInnerRural <- 4
JOuterRural <- 1
simIndices  <- 2:10

cat("Loading SPDE-simulated surveys ...\n")
load("savedOutput/simStudy1/simPopsSurveys.RData")
stopifnot(length(surveysDHS) >= max(simIndices))

for(i in simIndices) {
    outFile <- sprintf("savedOutput/simStudy1/intPtsDHS_simStudy1_%d_K16_21.RData", i)
    if(file.exists(outFile)) {
        cat("[sim ", i, "] already exists — skipping\n", sep="")
        next
    }

    cat("\n[sim ", i, "] building K=", KDHSu, "/", KDHSr, " integration points ...\n", sep="")
    datDHS <- surveysDHS[[i]]
    # column-name normalization to match the application data DHS naming
    if(!("n" %in% names(datDHS)) && "N" %in% names(datDHS)) datDHS$n <- datDHS$N
    if(!("y" %in% names(datDHS)) && "Z" %in% names(datDHS)) datDHS$y <- datDHS$Z

    t0 <- proc.time()[3]
    intPtsDHS <- makeAllIntegrationPointsDHS(
        cbind(datDHS$east, datDHS$north), datDHS$urban,
        areaNames=datDHS$subarea, popPrior=TRUE,
        numPointsUrban=KDHSu, numPointsRural=KDHSr,
        JInnerUrban=JInnerUrban, JInnerRural=JInnerRural,
        JOuterRural=JOuterRural,
        adminMap=adm2Full, outFile=outFile, saveOutput=FALSE)
    elapsed <- proc.time()[3] - t0
    cat(sprintf("[sim %d] built in %.1f min — wUrban=%d x %d, wRural=%d x %d\n",
                i, elapsed/60,
                nrow(intPtsDHS$wUrban), ncol(intPtsDHS$wUrban),
                nrow(intPtsDHS$wRural), ncol(intPtsDHS$wRural)))

    save(intPtsDHS, file=outFile)
    cat("[sim ", i, "] saved to ", outFile, "\n", sep="")
}

cat("\nAll done.\n")
