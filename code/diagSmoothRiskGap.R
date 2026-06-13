# Trace the +0.11 areal "smooth risk" gap to its source.
#
# Replicates simData1BYM2's setup for sim 1 (seed=123, same RNG path as the
# production simPopsSurveys_BYM2.RData), then compares four quantities:
#   (1) no-nugget risk:        pop-wtd mean of expit(logitRiskDraws)
#   (2) nugget-integrated:     pop-wtd mean of logitNormMean(logitRiskDraws, sigmaEps)
#   (3) stored pSmoothRisk:    simPop$areaPop$aggregationResults$pSmoothRisk
#   (4) stored prevalence:     simPop$areaPop$aggregationResults$pFineScalePrevalence
#
# If (3) == (1): the pipeline is NOT integrating the nugget (bug in what's
#   passed to logitNormMean — e.g. sigma 0).
# If (3) == (2) but (4) >> (3): the prevalence side is inflated instead
#   (e.g. nugget applied twice on the EA path).
# Also checks the regenerated sim matches the stored truth column 1 exactly
# (confirms we reproduced the production RNG path).

source("code/setup.R")

cat(sprintf("\n=== %s | smooth-risk gap diagnostic (BYM2 sim 1) ===\n", format(Sys.time())))

# ---- replicate simData1BYM2 defaults ----
sigmaBYM2    <- sqrt(0.5); phi <- 0.8; sigmaEpsilon <- sqrt(1.5)
beta0        <- -1.25; gamma <- 1; betaRest <- c(0, 0, 0, .5)
easpaDat     <- easpaNGAed
popMat       <- popMatNGAThresh
targetPopMat <- popMatNGAedThresh
poppsub      <- poppsubNGAThresh
nHHDHS       <- 25

set.seed(123)

popMat  <- popMat[order(popMat$subarea), ]
poppsub <- poppsub[order(poppsub$subarea), ]

cat("Constructing offset (same as simData1BYM2)...\n")
LLcoords   <- cbind(popMat$lon, popMat$lat)
tempDesMat <- getDesignMat(LLcoords, normalized = TRUE, useThreshPopMat = TRUE)
load("savedOutput/global/covariatesNorm.RData")
popVals <- extract(pop, LLcoords, method = "bilinear")
load("savedOutput/global/popMeanSDCal.RData")
normPop <- (log1p(popVals) - popMeanCal) / popSDCal
normPop[is.na(normPop)] <- min(normPop, na.rm = TRUE)
covRestMat <- tempDesMat[, -c(1:3, 7)]
covRestMat <- cbind(covRestMat, normPop = normPop)
offset     <- covRestMat %*% betaRest

out <- load("savedOutput/global/admFinalMat.RData")
graphObj <- admFinalMat
if(!("stratumMICS" %in% names(popMat)))
    popMat$stratumMICS <- adm2ToStratumMICS(popMat$subarea)

cat("Simulating population 1 (same RNG path as production run)...\n")
.silenceBrowser <- function() assign("browser", function(...) invisible(NULL), envir = globalenv())
.silenceBrowser()
simPop <- simPopBYM2(nsim = 1, easpa = easpaDat, popMat = popMat,
                     targetPopMat = targetPopMat, poppsub = poppsub,
                     graphObj = graphObj, areaCol = "stratumMICS",
                     sigmaBYM2 = sigmaBYM2, phi = phi,
                     sigmaEpsilon = sigmaEpsilon, gamma = gamma,
                     beta0 = beta0, seed = NULL,
                     nHHSampled = nHHDHS, stratifyByUrban = TRUE,
                     subareaLevel = TRUE, doFineScaleRisk = FALSE,
                     doSmoothRisk = TRUE, doSmoothRiskLogisticApprox = FALSE,
                     min1PerSubarea = TRUE, offset = offset, verbose = FALSE)

# ---- the four quantities ----
lrd <- as.numeric(simPop$logitRiskDraws)          # nPixels, includes offset?
w   <- targetPopMat$pop
w   <- w / sum(w)

noNug  <- sum(w * expit(lrd))
integ  <- sum(w * SUMMER::logitNormMean(
                  cbind(lrd, rep(sigmaEpsilon, length(lrd))),
                  logisticApprox = FALSE))

aggRes  <- simPop$areaPop$aggregationResults
storedSR <- mean(aggRes$pSmoothRisk,            na.rm = TRUE)
storedPv <- mean(aggRes$pFineScalePrevalence,   na.rm = TRUE)

cat("\n================ RESULTS (national means) ================\n")
cat(sprintf("(1) no-nugget        expit(eta):                 %.4f\n", noNug))
cat(sprintf("(2) nugget-integrated logitNormMean(eta, %.3f):  %.4f\n", sigmaEpsilon, integ))
cat(sprintf("(3) stored pSmoothRisk (mean over areas):        %.4f\n", storedSR))
cat(sprintf("(4) stored pFineScalePrevalence (mean):          %.4f\n", storedPv))

# Did we reproduce the production sim? Compare to saved truth col 1.
e <- new.env(); load("savedOutput/simStudy1/simPopsSurveys_BYM2.RData", envir = e)
repro <- max(abs(aggRes$pSmoothRisk - e$areaPops_smoothRisk[, 1]), na.rm = TRUE)
cat(sprintf("\nmax |fresh pSmoothRisk - stored col1| = %.3e  (%s)\n",
            repro, ifelse(repro < 1e-12, "EXACT reproduction",
                          ifelse(repro < 1e-6, "near-exact", "DIFFERENT RNG PATH"))))

cat("\nVerdict:\n")
if(abs(storedSR - noNug) < abs(storedSR - integ)) {
    cat("  stored pSmoothRisk matches the NO-NUGGET value -> nugget is NOT\n")
    cat("  being integrated in the truth pipeline (bug on the sigma passed\n")
    cat("  to logitNormMean, or logitRiskDraws path).\n")
} else {
    cat("  stored pSmoothRisk matches the INTEGRATED value -> truth side OK;\n")
    cat("  the prevalence/prediction side must be inflated instead.\n")
}
cat(sprintf("\n=== %s | done ===\n", format(Sys.time())))
