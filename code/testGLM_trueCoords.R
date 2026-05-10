#!/usr/bin/env Rscript
# GLM sanity check: extract covariates at true (un-jittered) cluster locations
# and fit a simple binomial GLM with intercept + urban + normPop.
# This bypasses the integration-point machinery entirely.

load("savedOutput/simStudy1/simPopsSurveys.RData")
load("savedOutput/global/popMatNGAThresh.RData")
load("savedOutput/global/popMeanSDCal.RData")

popEN = cbind(popMatNGAThresh[["east"]], popMatNGAThresh[["north"]])

matchPixels = function(ed) {
  ENcoords = cbind(ed[["east"]], ed[["north"]])
  apply(ENcoords, 1, function(pt) {
    closeE = (pt[1] > popEN[,1] - 2.5) & (pt[1] <= popEN[,1] + 2.5)
    closeN = (pt[2] > popEN[,2] - 2.5) & (pt[2] <= popEN[,2] + 2.5)
    w = which(closeE & closeN)
    if(length(w) == 1) return(w)
    if(length(w) > 1) return(w[1])
    dists = (pt[1] - popEN[,1])^2 + (pt[2] - popEN[,2])^2
    return(which.min(dists))
  })
}

fitGLM = function(ed, label) {
  cat("Matching", nrow(ed), "clusters to pixels...\n")
  pixelIs = matchPixels(ed)
  
  urbanTrue = as.numeric(popMatNGAThresh[["urban"]][pixelIs])
  popTrue = popMatNGAThresh[["pop"]][pixelIs]
  normPop = (log1p(popTrue) - popMeanCalThresh) / popSDCalThresh
  normPop[is.na(normPop)] = min(normPop, na.rm=TRUE)
  
  dat = data.frame(y = ed[["Z"]], n = ed[["N"]], urban = urbanTrue, normPop = normPop)
  fit = glm(cbind(y, n - y) ~ urban + normPop, data = dat, family = binomial)
  
  cat("\n=== GLM at true cluster locations -", label, "===\n")
  cat("Truth: alpha=-1.25, urban=1.00, normPop=0.50\n\n")
  print(summary(fit)[["coefficients"]])
  cat("\n")
  invisible(fit)
}

fitGLM(surveysMICS[[1]], "Sim 1")
fitGLM(surveysMICS[[2]], "Sim 2")
