# Diagnostic: try fitFED with the full MDM inputs (no fold subsetting) to isolate
# whether the segfault comes from buildFullMDMInputs or from subsetMDMInputs.

source("setup.R")
options(error=traceback)
source("code/modFED.R")
source("code/validationCommon.R")

cat("Building full MDM inputs ...\n")
fullInp = buildFullMDMInputs(KMICS=100, KDHSurb=16, KDHSrur=21)
cat("Done. Structure:\n")
cat("  nDHSurb =", fullInp$nDHSurb, "  nDHSrur =", fullInp$nDHSrur, "\n")
cat("  nMICSurb =", fullInp$nMICSurb, "  nMICSrur =", fullInp$nMICSrur, "\n")

cat("\nintPtsDHS components:\n")
cat("  covsUrb dim:", dim(fullInp$intPtsDHS$covsUrb), "  class:", class(fullInp$intPtsDHS$covsUrb), "\n")
cat("  covsRur dim:", dim(fullInp$intPtsDHS$covsRur), "  class:", class(fullInp$intPtsDHS$covsRur), "\n")
cat("  wUrban dim:", dim(fullInp$intPtsDHS$wUrban), "  class:", class(fullInp$intPtsDHS$wUrban), "\n")
cat("  wRural dim:", dim(fullInp$intPtsDHS$wRural), "  class:", class(fullInp$intPtsDHS$wRural), "\n")

cat("\nintPtsMICS components:\n")
cat("  XUrb dim:", dim(fullInp$intPtsMICS$XUrb), "  class:", class(fullInp$intPtsMICS$XUrb), "\n")
cat("  XRur dim:", dim(fullInp$intPtsMICS$XRur), "  class:", class(fullInp$intPtsMICS$XRur), "\n")
cat("  wUrban dim:", dim(fullInp$intPtsMICS$wUrban), "  class:", class(fullInp$intPtsMICS$wUrban), "\n")
cat("  wRural dim:", dim(fullInp$intPtsMICS$wRural), "  class:", class(fullInp$intPtsMICS$wRural), "\n")

cat("\ncovsUrb columns:", colnames(fullInp$intPtsDHS$covsUrb), "\n")
cat("covsRur columns:", colnames(fullInp$intPtsDHS$covsRur), "\n")
cat("XUrb columns:",   colnames(fullInp$intPtsMICS$XUrb), "\n")
cat("XRur columns:",   colnames(fullInp$intPtsMICS$XRur), "\n")

cat("\nareaidx vectors:\n")
cat("  areaidxlocUrbanDHS  length:", length(fullInp$areaidxlocUrbanDHS),
    "  range:", range(fullInp$areaidxlocUrbanDHS), "\n")
cat("  areaidxlocRuralDHS  length:", length(fullInp$areaidxlocRuralDHS),
    "  range:", range(fullInp$areaidxlocRuralDHS), "\n")
cat("  areaidxlocUrbanMICS length:", length(fullInp$areaidxlocUrbanMICS),
    "  range:", range(fullInp$areaidxlocUrbanMICS), "\n")
cat("  areaidxlocRuralMICS length:", length(fullInp$areaidxlocRuralMICS),
    "  range:", range(fullInp$areaidxlocRuralMICS), "\n")

cat("\nNow trying fitFED with full inputs (no subsetting) ...\n")
res = fitFED(datDHS=fullInp$datDHS, inputsMDM=fullInp,
             KDHSurb=fullInp$KDHSurb, KDHSrur=fullInp$KDHSrur,
             Qgh=10, getSDs=TRUE, verbose=FALSE)
cat("fitFED done. opt$convergence =", res$opt$convergence, "\n")
cat("alpha estimate:", res$opt$par[grep("alpha", names(res$opt$par))], "\n")
