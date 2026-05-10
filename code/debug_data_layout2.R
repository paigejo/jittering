#!/usr/bin/env Rscript
# Debug: check data structure of ed and intPtsDHS
setwd("c:/Users/jpaige/git/jittering")

load("savedOutput/simStudy1/simPopsSurveys_BYM2.RData")
ed = surveysDHS[[1]]

cat("class(ed):", class(ed), "\n")
cat("names(ed):", names(ed), "\n")
cat("class(ed[[1]]):", class(ed[[1]]), "\n")
cat("class(ed[['urban']]):", class(ed[["urban"]]), "\n")

# ed is likely a data.frame where $urban is itself a data.frame (nested)
# Let's flatten
if(is.data.frame(ed[["urban"]])) {
  cat("urban is a nested data.frame! Flattening...\n")
  # urban column is actually the full ed repeated as a sub-df
  isU = ed[["urban"]][["urban"]]  
} else {
  isU = ed[["urban"]]
}
cat("isU[1:10]:", isU[1:10], "\n")
cat("sum(isU):", sum(isU), "\n")

# Get east/north at true locations
trueEast = ed[["east"]]
if(is.data.frame(trueEast)) trueEast = trueEast[,1]
trueNorth = ed[["north"]]
if(is.data.frame(trueNorth)) trueNorth = trueNorth[,1]

cat("\nTrue east (urban 1:5):", trueEast[isU][1:5], "\n")

load("savedOutput/simStudy1/intPtsDHS_simStudy1_1.RData")
cat("IntPt1 x (urban 1:5):", intPtsDHS$xUrban[1:5,1], "\n")
cat("Diff east:", trueEast[isU][1:5] - intPtsDHS$xUrban[1:5,1], "\n")

cat("\nIntPt colnames xUrban:", colnames(intPtsDHS$xUrban), "\n")
# First column might be "pts1" which is the cluster center
cat("wUrban[1:5,1]:", intPtsDHS$wUrban[1:5,1], "\n")

# Check covariate at true location vs at first intpt
# covsUrb is (nUrb*K) x 6 in intpt-major order
# First intpt for obs i is at row i (rows 1:nUrb are intpt 1)
nU = sum(isU)
cat("\ncovsUrb[1:5,] (intpt 1, obs 1-5):\n")
print(intPtsDHS$covsUrb[1:5,])
cat("\ncovsUrb[nU+1:5,] (intpt 2, obs 1-5): nU=", nU, "\n")
print(intPtsDHS$covsUrb[(nU+1):(nU+5),])
