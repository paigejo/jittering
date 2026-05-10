#!/usr/bin/env Rscript
# Run BYM2 spatial models: fitMM (MICS-only), fitMD (DHS-only), fitMDM (combined)
# With GH nuggets, K=16/21 for DHS, K=100 for MICS, Q=10

library(TMB)
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")

Qgh <- 10
KMICS <- 100
KDHSu <- 16
KDHSr <- 21

out <- list()

cat("\n========== Running fitMM (MICS-only BYM2+GH) ==========\n")
t0 <- proc.time()
res_mm <- fitMM(datDHS=ed, datMICS=edMICS, KMICS=KMICS, Qgh=Qgh, getSDs=TRUE, verbose=FALSE)
t_mm <- (proc.time() - t0)[3]
cat("MODEL: fitMM | CONV:", res_mm$opt$convergence, "| NLL:", sprintf("%.6f", res_mm$opt$objective), "| WALLTIME:", sprintf("%.2f", t_mm), "\n")
if(!is.null(res_mm$TMBsd)) {
  sumtab <- summary(res_mm$TMBsd, "fixed")
  for(i in 1:nrow(sumtab)) cat("fitMM |", rownames(sumtab)[i], "|", sprintf("%.6f", sumtab[i,1]), "|", sprintf("%.6f", sumtab[i,2]), "\n")
}
out$fitMM <- res_mm

cat("\n========== Running fitMD (DHS-only BYM2+GH) ==========\n")
t0 <- proc.time()
res_md <- fitMD(datDHS=ed, datMICS=edMICS, KDHSurb=KDHSu, KDHSrur=KDHSr, Qgh=Qgh, getSDs=TRUE, verbose=FALSE)
t_md <- (proc.time() - t0)[3]
cat("MODEL: fitMD | CONV:", res_md$opt$convergence, "| NLL:", sprintf("%.6f", res_md$opt$objective), "| WALLTIME:", sprintf("%.2f", t_md), "\n")
if(!is.null(res_md$TMBsd)) {
  sumtab <- summary(res_md$TMBsd, "fixed")
  for(i in 1:nrow(sumtab)) cat("fitMD |", rownames(sumtab)[i], "|", sprintf("%.6f", sumtab[i,1]), "|", sprintf("%.6f", sumtab[i,2]), "\n")
}
out$fitMD <- res_md

cat("\n========== Running fitMDM (DHS+MICS BYM2+GH) ==========\n")
t0 <- proc.time()
res_mdm <- fitMDM(datDHS=ed, datMICS=edMICS, KMICS=KMICS, KDHSurb=KDHSu, KDHSrur=KDHSr, Qgh=Qgh, getSDs=TRUE, verbose=FALSE)
t_mdm <- (proc.time() - t0)[3]
cat("MODEL: fitMDM | CONV:", res_mdm$opt$convergence, "| NLL:", sprintf("%.6f", res_mdm$opt$objective), "| WALLTIME:", sprintf("%.2f", t_mdm), "\n")
if(!is.null(res_mdm$TMBsd)) {
  sumtab <- summary(res_mdm$TMBsd, "fixed")
  for(i in 1:nrow(sumtab)) cat("fitMDM |", rownames(sumtab)[i], "|", sprintf("%.6f", sumtab[i,1]), "|", sprintf("%.6f", sumtab[i,2]), "\n")
}
out$fitMDM <- res_mdm

# Save
if(!dir.exists("savedOutput/global")) dir.create("savedOutput/global", recursive=TRUE)
save(out, file="savedOutput/global/BYM2_models_test_results.RData")
cat("\nSaved results to savedOutput/global/BYM2_models_test_results.RData\n")
cat("Done.\n")
