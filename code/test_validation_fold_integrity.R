# Per-(stratum x urbanicity) fold integrity check.
#
# For each DHS fold k = 1..10:
#   - subset to in-sample (fold != k) and held-out (fold == k)
#   - compute per-(stratum, urban) cluster counts in the subset using areaidxloc*DHS
#   - compute the expected per-(stratum, urban) counts from fullInp$datDHS$fold
#   - assert they agree exactly
# Same for MICS folds.
#
# This catches bugs where datMICS row order, makeInputsMDM internal sort, area index
# mapping, or subsetMDMInputs slicing get out of sync.

source("setup.R")
options(error=traceback)
source("code/validationCommon.R")

cat("Loading full inputs ...\n")
fullCache = "savedOutput/validation/fullMDMInputs_FE.RData"
if(file.exists(fullCache)) {
    load(fullCache)
} else {
    fullInp = buildFullMDMInputs(KMICS=100, KDHSurb=16, KDHSrur=21)
    save(fullInp, file=fullCache)
}

# Map area index (0-based) -> stratum name.
strataNames = fullInp$intPtsMICS$strataMICS
nStrata     = length(strataNames)
cat("nStrata =", nStrata, "\n")

# Pre-compute stratum names for each cluster in fullInp$datDHS / datMICS.
dhsStratum  = adm2ToStratumMICS(fullInp$datDHS$subarea)
micsStratum = fullInp$datMICS$Stratum

# Sanity: every DHS cluster's stratum should be one of strataNames.
stopifnot(all(dhsStratum %in% strataNames))
stopifnot(all(micsStratum %in% strataNames))

# Reference table: full per-(stratum, urban, fold) counts from datDHS / datMICS.
fullDhsTab  = with(fullInp$datDHS,
                   table(stratum=factor(dhsStratum, levels=strataNames),
                         urban=urban, fold=fold))
fullMicsTab = with(fullInp$datMICS,
                   table(stratum=factor(micsStratum, levels=strataNames),
                         urban=urban, fold=fold))

# ---- helper: per-(stratum,urban) counts derived from a subset ----------------
subsetCounts = function(subset, side=c("DHS","MICS")) {
    side = match.arg(side)
    if(side == "DHS") {
        nU = length(subset$areaidxlocUrbanDHS)
        nR = length(subset$areaidxlocRuralDHS)
        # area index is 0-based; convert to factor with all strata levels for stability
        cuU = if(nU > 0) tabulate(subset$areaidxlocUrbanDHS + 1, nbins=nStrata) else integer(nStrata)
        cuR = if(nR > 0) tabulate(subset$areaidxlocRuralDHS + 1, nbins=nStrata) else integer(nStrata)
    } else {
        nU = length(subset$areaidxlocUrbanMICS)
        nR = length(subset$areaidxlocRuralMICS)
        cuU = if(nU > 0) tabulate(subset$areaidxlocUrbanMICS + 1, nbins=nStrata) else integer(nStrata)
        cuR = if(nR > 0) tabulate(subset$areaidxlocRuralMICS + 1, nbins=nStrata) else integer(nStrata)
    }
    list(urb=cuU, rur=cuR)
}

# ---- DHS folds ---------------------------------------------------------------
cat("\n========== DHS fold integrity ==========\n")
nDhsErr = 0
for(k in 1:10) {
    dhsKeep    = fullInp$datDHS$fold != k
    dhsHeld    = !dhsKeep
    micsKeep   = rep(TRUE,  nrow(fullInp$datMICS))
    micsHeld   = rep(FALSE, nrow(fullInp$datMICS))

    inSamp     = subsetMDMInputs(fullInp, dhsKeep,  micsKeep)
    heldOut    = subsetMDMInputs(fullInp, dhsHeld,  micsHeld)

    # In-sample DHS counts
    expIn      = fullDhsTab
    expIn[,,k] = 0  # zero out the held-out fold
    expInUrb   = apply(expIn[, "TRUE",  ], 1, sum)
    expInRur   = apply(expIn[, "FALSE", ], 1, sum)

    cIn  = subsetCounts(inSamp,  "DHS")
    cOut = subsetCounts(heldOut, "DHS")

    expOutUrb = fullDhsTab[, "TRUE",  k]
    expOutRur = fullDhsTab[, "FALSE", k]

    okInUrb  = all(cIn$urb  == expInUrb)
    okInRur  = all(cIn$rur  == expInRur)
    okOutUrb = all(cOut$urb == expOutUrb)
    okOutRur = all(cOut$rur == expOutRur)

    status = if(okInUrb && okInRur && okOutUrb && okOutRur) "OK" else "MISMATCH"
    cat(sprintf("  fold %2d  inUrb=%4d inRur=%4d outUrb=%3d outRur=%3d  %s\n",
                k, sum(cIn$urb), sum(cIn$rur), sum(cOut$urb), sum(cOut$rur), status))

    if(status == "MISMATCH") {
        nDhsErr = nDhsErr + 1
        if(!okOutUrb) {
            diff = cOut$urb - expOutUrb
            badI = which(diff != 0)
            cat("    held-out urban mismatches at strata:",
                paste(strataNames[badI], "(", diff[badI], ")", collapse=", "), "\n")
        }
        if(!okOutRur) {
            diff = cOut$rur - expOutRur
            badI = which(diff != 0)
            cat("    held-out rural mismatches at strata:",
                paste(strataNames[badI], "(", diff[badI], ")", collapse=", "), "\n")
        }
    }
}

# ---- MICS folds ---------------------------------------------------------------
cat("\n========== MICS fold integrity ==========\n")
nMicsErr = 0
for(k in 1:10) {
    dhsKeep    = rep(TRUE,  nrow(fullInp$datDHS))
    dhsHeld    = rep(FALSE, nrow(fullInp$datDHS))
    micsKeep   = fullInp$datMICS$fold != k
    micsHeld   = !micsKeep

    inSamp     = subsetMDMInputs(fullInp, dhsKeep,  micsKeep)
    heldOut    = subsetMDMInputs(fullInp, dhsHeld,  micsHeld)

    expIn      = fullMicsTab
    expIn[,,k] = 0
    expInUrb   = apply(expIn[, "TRUE",  ], 1, sum)
    expInRur   = apply(expIn[, "FALSE", ], 1, sum)

    cIn  = subsetCounts(inSamp,  "MICS")
    cOut = subsetCounts(heldOut, "MICS")

    expOutUrb = fullMicsTab[, "TRUE",  k]
    expOutRur = fullMicsTab[, "FALSE", k]

    okInUrb  = all(cIn$urb  == expInUrb)
    okInRur  = all(cIn$rur  == expInRur)
    okOutUrb = all(cOut$urb == expOutUrb)
    okOutRur = all(cOut$rur == expOutRur)

    status = if(okInUrb && okInRur && okOutUrb && okOutRur) "OK" else "MISMATCH"
    cat(sprintf("  fold %2d  inUrb=%4d inRur=%4d outUrb=%3d outRur=%3d  %s\n",
                k, sum(cIn$urb), sum(cIn$rur), sum(cOut$urb), sum(cOut$rur), status))

    if(status == "MISMATCH") {
        nMicsErr = nMicsErr + 1
        if(!okOutUrb) {
            diff = cOut$urb - expOutUrb
            badI = which(diff != 0)
            cat("    held-out urban mismatches at strata:",
                paste(strataNames[badI], "(", diff[badI], ")", collapse=", "), "\n")
        }
        if(!okOutRur) {
            diff = cOut$rur - expOutRur
            badI = which(diff != 0)
            cat("    held-out rural mismatches at strata:",
                paste(strataNames[badI], "(", diff[badI], ")", collapse=", "), "\n")
        }
    }
}

cat("\nSummary: DHS errors =", nDhsErr, " / 10   MICS errors =", nMicsErr, " / 10\n")
if(nDhsErr == 0 && nMicsErr == 0) {
    cat("All per-(stratum, urbanicity) fold counts match. Pipeline integrity confirmed.\n")
} else {
    cat("MISMATCH detected — see per-fold diffs above.\n")
}
