# Quick inspection of the K values in the simStudy1 DHS integration-point files.
source("setup.R")

f1 <- "savedOutput/simStudy1/intPtsDHS_simStudy1_1.RData"
f2 <- "savedOutput/simStudy1/intPtsDHS_simStudy1_1_K16_21.RData"

e <- new.env(); load(f1, envir=e); ip1 <- e$intPtsDHS
e <- new.env(); load(f2, envir=e); ip2 <- e$intPtsDHS

cat("default sim1:\n")
cat("  wUrban:", dim(ip1$wUrban), "  wRural:", dim(ip1$wRural), "\n")
cat("  covsUrb rows:", nrow(ip1$covsUrb), "  covsRur rows:", nrow(ip1$covsRur), "\n")
cat("K16_21 sim1:\n")
cat("  wUrban:", dim(ip2$wUrban), "  wRural:", dim(ip2$wRural), "\n")
cat("  covsUrb rows:", nrow(ip2$covsUrb), "  covsRur rows:", nrow(ip2$covsRur), "\n")

# Also peek at sim 2 default
f3 <- "savedOutput/simStudy1/intPtsDHS_simStudy1_2.RData"
e <- new.env(); load(f3, envir=e); ip3 <- e$intPtsDHS
cat("default sim2:\n")
cat("  wUrban:", dim(ip3$wUrban), "  wRural:", dim(ip3$wRural), "\n")
