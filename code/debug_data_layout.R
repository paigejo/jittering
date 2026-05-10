#!/usr/bin/env Rscript
# Quick diagnostic of intPtsDHS data structure
setwd("c:/Users/jpaige/git/jittering")
load("savedOutput/simStudy1/simPopsSurveys_BYM2.RData")
ed = surveysDHS[[1]]
load("savedOutput/simStudy1/intPtsDHS_simStudy1_1.RData")

nU = sum(ed$urban)
nR = sum(!ed$urban)
cat("nUrb:", nU, "nRur:", nR, "\n\n")

cat("=== covsUrb ===\n")
cat("class:", class(intPtsDHS$covsUrb), "\n")
cat("dim:", nrow(intPtsDHS$covsUrb), "x", ncol(intPtsDHS$covsUrb), "\n")
cat("colnames:", colnames(intPtsDHS$covsUrb), "\n")
cat("First 3 rows:\n")
print(intPtsDHS$covsUrb[1:3,])
cat("\nRows 11-13 (if obs-major, these are obs 2 int pts):\n")
print(intPtsDHS$covsUrb[11:13,])
cat("\nRow 573 (if intpt-major, this is intpt 1, obs 1):\n")
print(intPtsDHS$covsUrb[573,])

cat("\n=== wUrban ===\n")
cat("class:", class(intPtsDHS$wUrban), "\n")
cat("dim:", nrow(intPtsDHS$wUrban), "x", ncol(intPtsDHS$wUrban), "\n")
cat("wUrban[1:5,]:\n")
print(as.matrix(intPtsDHS$wUrban)[1:5,])
cat("\nRow sums (first 10):", apply(as.matrix(intPtsDHS$wUrban), 1, sum)[1:10], "\n")

cat("\n=== wRural ===\n")
cat("class:", class(intPtsDHS$wRural), "\n")
cat("dim:", nrow(intPtsDHS$wRural), "x", ncol(intPtsDHS$wRural), "\n")
cat("Row sums (first 10):", apply(as.matrix(intPtsDHS$wRural), 1, sum)[1:10], "\n")

# Check: does the 'urb' covariate vary across integration points for a single obs?
cat("\n=== Covariate variation across int pts ===\n")
# For obs 1, the covariates at different int pts should differ (different jittered locations)
# If covsUrb is stacked by int pt: rows 0, 572, 1144, ... are obs 1 at different int pts  
# If covsUrb is stacked by obs: rows 0..10 are obs 1 at 11 int pts
cat("Assuming intpt-major (C++ expectation): obs 1 across int pts:\n")
idx_intpt_major = seq(1, by=nU, length.out=11)
cat("  Indices:", idx_intpt_major, "\n")
cat("  'urb' col values:", intPtsDHS$covsUrb[idx_intpt_major, 2], "\n")
cat("  'pop' col values:", intPtsDHS$covsUrb[idx_intpt_major, 6], "\n")

cat("\nAssuming obs-major: obs 1 across int pts:\n")
idx_obs_major = 1:11
cat("  Indices:", idx_obs_major, "\n")
cat("  'urb' col values:", intPtsDHS$covsUrb[idx_obs_major, 2], "\n")
cat("  'pop' col values:", intPtsDHS$covsUrb[idx_obs_major, 6], "\n")

# Also check areas — same observation should map to nearby areas
cat("\n=== areasUrban ===\n")
cat("First 20:", intPtsDHS$areasUrban[1:20], "\n")
cat("length:", length(intPtsDHS$areasUrban), "\n")

# Check how Laplace template handles the same data  
cat("\n=== Compare with Laplace template's AprojUrbanDHS ===\n")
library(Matrix)
source("code/utils.R")
load("savedOutput/global/admFinal.RData")
AUrbDHS = makeApointToArea(intPtsDHS$areasUrban, admFinal$NAME_FINAL)
cat("AUrbDHS dim:", nrow(AUrbDHS), "x", ncol(AUrbDHS), "\n")
AUrbDHS_t = t(AUrbDHS)
cat("t(AUrbDHS) dim:", nrow(AUrbDHS_t), "x", ncol(AUrbDHS_t), "\n")
