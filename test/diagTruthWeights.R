# Why does stored pSmoothRisk (~0.34) differ from pop-weighted integrated
# risk (~0.46) computed from the SAME field? Candidates:
#   W1 truth aggregation effectively flat-weighted
#   W2 popMat sorted by subarea but targetPopMat NOT sorted -> row misalignment
#   W3 urban-share mismatch: EA frame (easpa) vs pixel classification
#
# No model fits — just regenerate sim 1's risk field (same RNG path as
# diagSmoothRiskGap.R) and tabulate.

source("code/setup.R")
cat(sprintf("\n=== %s | truth-weighting diagnostic (BYM2 sim 1) ===\n", format(Sys.time())))

# ---- replicate simData1BYM2 setup ----
sigmaBYM2 <- sqrt(0.5); phi <- 0.8; sigmaEpsilon <- sqrt(1.5)
beta0 <- -1.25; gamma <- 1; betaRest <- c(0, 0, 0, .5)
popMat0       <- popMatNGAThresh          # UNsorted originals
targetPopMat0 <- popMatNGAedThresh

set.seed(123)
popMat  <- popMat0[order(popMat0$subarea), ]
poppsub <- poppsubNGAThresh[order(poppsubNGAThresh$subarea), ]

# W2 check: is popMat already sorted by subarea (=> sorting is a no-op)?
perm <- order(popMat0$subarea)
cat(sprintf("\npopMat pre-sorted by subarea?  %s\n",
            ifelse(identical(perm, seq_len(nrow(popMat0))), "YES (no-op)",
                   sprintf("NO — sorting PERMUTES rows (first moved row: %d)",
                           which(perm != seq_len(nrow(popMat0)))[1]))))
cat(sprintf("rows where sorted popMat$subarea != targetPopMat0$subarea: %d / %d\n",
            sum(popMat$subarea != targetPopMat0$subarea), nrow(popMat)))

cat("Constructing offset...\n")
LLcoords   <- cbind(popMat$lon, popMat$lat)
tempDesMat <- getDesignMat(LLcoords, normalized = TRUE, useThreshPopMat = TRUE)
load("savedOutput/global/covariatesNorm.RData")
popVals <- extract(pop, LLcoords, method = "bilinear")
load("savedOutput/global/popMeanSDCal.RData")
normPop <- (log1p(popVals) - popMeanCal) / popSDCal
normPop[is.na(normPop)] <- min(normPop, na.rm = TRUE)
covRestMat <- cbind(tempDesMat[, -c(1:3, 7)], normPop = normPop)
offset     <- covRestMat %*% betaRest

out <- load("savedOutput/global/admFinalMat.RData")
if(!("stratumMICS" %in% names(popMat)))
    popMat$stratumMICS <- adm2ToStratumMICS(popMat$subarea)

assign("browser", function(...) invisible(NULL), envir = globalenv())
cat("Simulating population 1...\n")
simPop <- simPopBYM2(nsim = 1, easpa = easpaNGAed, popMat = popMat,
                     targetPopMat = targetPopMat0, poppsub = poppsub,
                     graphObj = admFinalMat, areaCol = "stratumMICS",
                     sigmaBYM2 = sigmaBYM2, phi = phi,
                     sigmaEpsilon = sigmaEpsilon, gamma = gamma,
                     beta0 = beta0, seed = NULL,
                     nHHSampled = 25, stratifyByUrban = TRUE,
                     subareaLevel = TRUE, doFineScaleRisk = FALSE,
                     doSmoothRisk = TRUE, doSmoothRiskLogisticApprox = FALSE,
                     min1PerSubarea = TRUE, offset = offset, verbose = FALSE)

lrd   <- as.numeric(simPop$logitRiskDraws)     # aligned with SORTED popMat
integ <- SUMMER::logitNormMean(cbind(lrd, rep(sigmaEpsilon, length(lrd))),
                               logisticApprox = FALSE)

# targetPopMat re-sorted to MATCH sorted popMat (the alignment SUMMER needs)
targetSorted <- targetPopMat0[perm, ]

wEq   <- rep(1, length(lrd))
wPop  <- popMat$pop                      # sorted popMat's own pop
wTarM <- targetPopMat0$pop               # MISaligned (as production passes it)
wTarA <- targetSorted$pop                # aligned

wm <- function(w) sum(w * integ, na.rm = TRUE) / sum(w[!is.na(integ)])
cat("\n========== national integrated-risk means under weightings ==========\n")
cat(sprintf("equal weights                : %.4f\n", wm(wEq)))
cat(sprintf("popMat$pop (sorted, aligned) : %.4f\n", wm(wPop)))
cat(sprintf("targetPop UNALIGNED (as prod): %.4f\n", wm(wTarM)))
cat(sprintf("targetPop ALIGNED            : %.4f\n", wm(wTarA)))

aggRes <- simPop$areaPop$aggregationResults
cat(sprintf("\nstored truth: mean pSmoothRisk = %.4f   mean prevalence = %.4f\n",
            mean(aggRes$pSmoothRisk, na.rm = TRUE),
            mean(aggRes$pFineScalePrevalence, na.rm = TRUE)))

# W3: urban shares
cat("\n========== urban-share comparison ==========\n")
cat(sprintf("popMat:    pop share urban = %.3f   pixel share urban = %.3f\n",
            sum(popMat$pop[popMat$urban])/sum(popMat$pop),
            mean(popMat$urban)))
cat(sprintf("targetPop: pop share urban = %.3f\n",
            sum(targetSorted$pop[targetSorted$urban])/sum(targetSorted$pop)))
cat(sprintf("easpaNGAed: EA share urban = %.3f   target-pop share urban = %.3f\n",
            sum(easpaNGAed$EAUrb)/sum(easpaNGAed$EATotal),
            sum(easpaNGAed$popUrb)/sum(easpaNGAed$popTotal)))

# EA-level reality for this sim
ea <- simPop$eaPop$eaDatList[[1]]
cat(sprintf("\nsimulated EAs: n = %d   urban share of EAs = %.3f   of EA target pop = %.3f\n",
            nrow(ea), mean(ea$urban), sum(ea$N[ea$urban])/sum(ea$N)))
cat(sprintf("EA-pop-weighted mean EA risk (Z/N national) = %.4f\n",
            sum(ea$Z)/sum(ea$N)))

cat(sprintf("\n=== %s | done ===\n", format(Sys.time())))
