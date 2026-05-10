# Diagnostic 2: subset MDM inputs and try fitFED on the subset.
# Trace the subset structure and dimensions.

source("setup.R")
options(error=traceback)
source("code/modFED.R")
source("code/validationCommon.R")

cat("Building full inputs ...\n")
fullCache = "savedOutput/validation/fullMDMInputs_FE.RData"
if(file.exists(fullCache)) {
    load(fullCache)
} else {
    fullInp = buildFullMDMInputs(KMICS=100, KDHSurb=16, KDHSrur=21)
    if(!dir.exists("savedOutput/validation")) dir.create("savedOutput/validation", recursive=TRUE)
    save(fullInp, file=fullCache)
}

# Make a fake DHS keep vector — drop first ~10% of urban DHS clusters and first ~10% of rural.
set.seed(1)
dhsKeep  = rep(TRUE, nrow(fullInp$datDHS))
dhsKeep[sample(which(fullInp$datDHS$urban),  60)]  = FALSE
dhsKeep[sample(which(!fullInp$datDHS$urban), 80)]  = FALSE
micsKeep = rep(TRUE, nrow(fullInp$datMICS))

cat("Subsetting ...\n")
sub = subsetMDMInputs(fullInp, dhsKeep, micsKeep)

cat("\n=== subset structure ===\n")
cat("nDHSurb (subset):", sub$nDHSurb, "\n")
cat("nDHSrur (subset):", sub$nDHSrur, "\n")
cat("ysUrbDHS length:", length(sub$ysUrbDHS), "\n")
cat("ysRurDHS length:", length(sub$ysRurDHS), "\n")
cat("intPtsDHS$covsUrb dim:", dim(sub$intPtsDHS$covsUrb),
    " (expected", sub$nDHSurb, "*", fullInp$KDHSurb, " =", sub$nDHSurb*fullInp$KDHSurb, ")\n")
cat("intPtsDHS$covsRur dim:", dim(sub$intPtsDHS$covsRur),
    " (expected", sub$nDHSrur, "*", fullInp$KDHSrur, " =", sub$nDHSrur*fullInp$KDHSrur, ")\n")
cat("intPtsDHS$wUrban dim:", dim(sub$intPtsDHS$wUrban), "\n")
cat("intPtsDHS$wRural dim:", dim(sub$intPtsDHS$wRural), "\n")
cat("class(covsUrb):",  class(sub$intPtsDHS$covsUrb), "\n")
cat("class(wUrban):",   class(sub$intPtsDHS$wUrban), "\n")
cat("areaidxlocUrbanDHS length:", length(sub$areaidxlocUrbanDHS), "\n")
cat("areaidxlocRuralDHS length:", length(sub$areaidxlocRuralDHS), "\n")

cat("\nMICS pieces (untouched):\n")
cat("intPtsMICS$XUrb dim:", dim(sub$intPtsMICS$XUrb), "\n")
cat("intPtsMICS$XRur dim:", dim(sub$intPtsMICS$XRur), "\n")
cat("intPtsMICS$wUrban dim:", dim(sub$intPtsMICS$wUrban), "\n")

cat("\nNow trying fitFED on subset ...\n")
res = fitFED(datDHS=fullInp$datDHS, inputsMDM=sub,
             KDHSurb=fullInp$KDHSurb, KDHSrur=fullInp$KDHSrur,
             Qgh=10, getSDs=FALSE, verbose=FALSE)
cat("fitFED done. opt$convergence =", res$opt$convergence, "\n")
