#!/usr/bin/env Rscript
# Run quick FE model tests: fitFEM (MICS-only), fitFED (DHS-only), fitFEMD (combined)
# Saves results to savedOutput/global/FE_models_test_results.RData

library(TMB)
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
source("code/modFEM.R")
source("code/modFED.R")
source("code/modFEMD.R")

Qgh <- 10
KMICS <- 100
KDHSu <- 16
KDHSr <- 21

out <- list()

cat("\nRunning fitFEM (MICS-only)...\n")
res_mics <- fitFEM(datDHS=ed, datMICS=edMICS, KMICS=KMICS, Qgh=Qgh,
                   fixedEffectsOnly=TRUE, getSDs=TRUE, verbose=FALSE)
cat("MODEL: fitFEM | CONV:", res_mics$opt$convergence, "| NLL:", sprintf("%.6f", res_mics$opt$objective), "| TIME:", sprintf("%.2f", res_mics$totalTime), "| SDTIME:", sprintf("%.2f", res_mics$sdTime), "\n")
if(!is.null(res_mics$TMBsd)) {
  sumtab <- summary(res_mics$TMBsd, "fixed")
  for(i in 1:nrow(sumtab)) cat("fitFEM |", rownames(sumtab)[i], "|", sprintf("%.6f", sumtab[i,1]), "|", sprintf("%.6f", sumtab[i,2]), "\n")
}
out$fitFEM <- res_mics

cat("\nRunning fitFED (DHS-only)...\n")
res_dhs <- fitFED(datDHS=ed, inputsMDM=NULL, KDHSurb=KDHSu, KDHSrur=KDHSr, Qgh=Qgh,
                  fixedEffectsOnly=TRUE, getSDs=TRUE, verbose=FALSE)
cat("MODEL: fitFED | CONV:", res_dhs$opt$convergence, "| NLL:", sprintf("%.6f", res_dhs$opt$objective), "| TIME:", sprintf("%.2f", res_dhs$totalTime), "| SDTIME:", sprintf("%.2f", res_dhs$sdTime), "\n")
if(!is.null(res_dhs$TMBsd)) {
  sumtab <- summary(res_dhs$TMBsd, "fixed")
  for(i in 1:nrow(sumtab)) cat("fitFED |", rownames(sumtab)[i], "|", sprintf("%.6f", sumtab[i,1]), "|", sprintf("%.6f", sumtab[i,2]), "\n")
}
out$fitFED <- res_dhs

cat("\nRunning fitFEMD (combined DHS+MICS)...\n")
res_comb <- fitFEMD(datDHS=ed, datMICS=edMICS, inputsMDM=NULL, KMICS=KMICS, KDHSurb=KDHSu, KDHSrur=KDHSr, Qgh=Qgh,
                    fixedEffectsOnly=TRUE, getSDs=TRUE, verbose=FALSE)
cat("MODEL: fitFEMD | CONV:", res_comb$opt$convergence, "| NLL:", sprintf("%.6f", res_comb$opt$objective), "| TIME:", sprintf("%.2f", res_comb$totalTime), "| SDTIME:", sprintf("%.2f", res_comb$sdTime), "\n")
if(!is.null(res_comb$TMBsd)) {
  sumtab <- summary(res_comb$TMBsd, "fixed")
  for(i in 1:nrow(sumtab)) cat("fitFEMD |", rownames(sumtab)[i], "|", sprintf("%.6f", sumtab[i,1]), "|", sprintf("%.6f", sumtab[i,2]), "\n")
}
out$fitFEMD <- res_comb

# Save
if(!dir.exists("savedOutput/global")) dir.create("savedOutput/global", recursive=TRUE)
save(out, file="savedOutput/global/FE_models_test_results.RData")
cat("\nSaved results to savedOutput/global/FE_models_test_results.RData\n")
cat("Done.\n")
