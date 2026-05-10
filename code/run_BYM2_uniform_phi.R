#!/usr/bin/env Rscript
# Run fitMM and fitMD with uniform prior on phi

library(TMB)
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")

Qgh <- 10; KMICS <- 100; KDHSu <- 16; KDHSr <- 21

cat("\n========== fitMM uniform phi prior ==========\n")
t0 <- proc.time()
res_mm_unif <- fitMM(datDHS=ed, datMICS=edMICS, KMICS=KMICS, Qgh=Qgh,
                     uniform_phi_prior=TRUE, getSDs=TRUE, verbose=FALSE)
t_mm <- (proc.time() - t0)[3]
cat("MODEL: fitMM_unif | CONV:", res_mm_unif$opt$convergence,
    "| NLL:", sprintf("%.6f", res_mm_unif$opt$objective),
    "| WALLTIME:", sprintf("%.2f", t_mm), "\n")
if(!is.null(res_mm_unif$TMBsd)) {
  st <- summary(res_mm_unif$TMBsd, "fixed")
  for(i in 1:nrow(st)) cat("fitMM_unif |", rownames(st)[i], "|",
                             sprintf("%.6f", st[i,1]), "|", sprintf("%.6f", st[i,2]), "\n")
}

cat("\n========== fitMD uniform phi prior ==========\n")
t0 <- proc.time()
res_md_unif <- fitMD(datDHS=ed, datMICS=edMICS, KDHSurb=KDHSu, KDHSrur=KDHSr, Qgh=Qgh,
                     uniform_phi_prior=TRUE, getSDs=TRUE, verbose=FALSE)
t_md <- (proc.time() - t0)[3]
cat("MODEL: fitMD_unif | CONV:", res_md_unif$opt$convergence,
    "| NLL:", sprintf("%.6f", res_md_unif$opt$objective),
    "| WALLTIME:", sprintf("%.2f", t_md), "\n")
if(!is.null(res_md_unif$TMBsd)) {
  st <- summary(res_md_unif$TMBsd, "fixed")
  for(i in 1:nrow(st)) cat("fitMD_unif |", rownames(st)[i], "|",
                             sprintf("%.6f", st[i,1]), "|", sprintf("%.6f", st[i,2]), "\n")
}

out_unif <- list(fitMM=res_mm_unif, fitMD=res_md_unif)
if(!dir.exists("savedOutput/global")) dir.create("savedOutput/global", recursive=TRUE)
save(out_unif, file="savedOutput/global/BYM2_uniform_phi_results.RData")
cat("\nSaved to savedOutput/global/BYM2_uniform_phi_results.RData\n")
cat("Done.\n")
