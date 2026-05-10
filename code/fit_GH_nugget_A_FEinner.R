#!/usr/bin/env Rscript
# Config A with all FEs as inner (random) parameters
# GH Q=10 for nuggets, outer = log_tauEps only

library(TMB)
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
source("code/modM_MIIDonly.R")

load("savedOutput/simStudy1/simPopsSurveys_BYM2.RData")
edMICS = surveysMICS[[1]]

nameTab = rbind(c("N", "ns"), c("N", "n"), c("Z", "ys"), c("Z", "y"))
for(i in 1:nrow(nameTab)) {
  fromN = nameTab[i,1]; toN = nameTab[i,2]
  if(!(toN %in% names(edMICS))) edMICS[[toN]] = edMICS[[fromN]]
  if(!(toN %in% names(ed))) ed[[toN]] = ed[[fromN]]
}
if(!("Stratum" %in% names(edMICS)))
  edMICS$Stratum = adm2ToStratumMICS(edMICS$subarea)

inputsMDM = makeInputsMDM(datDHS=ed, datMICS=edMICS,
                          KMICS=100, admMICS=admFinal, adm2DHS=adm2Full)

allCovNames = colnames(inputsMDM$intPtsMICS$XUrb)
keepIdx = which(allCovNames %in% c("urban", "normPop"))
inputsMDM$intPtsMICS$XUrb = inputsMDM$intPtsMICS$XUrb[, keepIdx, drop=FALSE]
inputsMDM$intPtsMICS$XRur = inputsMDM$intPtsMICS$XRur[, keepIdx, drop=FALSE]

thisEnv = environment()
list2env(inputsMDM, envir=thisEnv)

lambdaTau = getLambdaPCprec(u=1, alpha=0.1)
lambdaTauEps = getLambdaPCprec(u=1, alpha=0.1)
nAreas = max(c(areaidxlocUrbanMICS, areaidxlocRuralMICS)) + 1

gh = fastGHQuad::gaussHermiteData(10)

initUrbP = sum(ysUrbMICS)/sum(nsUrbMICS)
initRurP = sum(ysRurMICS)/sum(nsRurMICS)
initAlpha = log(initRurP/(1-initRurP))
initBeta1 = log(initUrbP/(1-initUrbP)) - initAlpha

data_gh = list(
  y_iUrbanMICS=ysUrbMICS, y_iRuralMICS=ysRurMICS,
  n_iUrbanMICS=nsUrbMICS, n_iRuralMICS=nsRurMICS,
  areaidxlocUrbanMICS=areaidxlocUrbanMICS, areaidxlocRuralMICS=areaidxlocRuralMICS,
  X_betaUrbanMICS=intPtsMICS$XUrb, X_betaRuralMICS=intPtsMICS$XRur,
  wUrbanMICS=intPtsMICS$wUrban, wRuralMICS=intPtsMICS$wRural,
  nAreas=nAreas,
  gh_nodes=gh$x, gh_weights=gh$w,
  alpha_pri=c(0, 100^2), beta_pri=c(0, sqrt(1000)),
  lambdaTau=lambdaTau, lambdaTauEps=lambdaTauEps
)

unloadDynlibs()
dyn.load(dynlib("code/modM_MIIDonly_GH"))

params_A = list(
  log_tau=0, log_tauEps=0,
  alpha=initAlpha, beta_urban=initBeta1, beta_normPop=0,
  u_spatial=rep(0, nAreas)
)
# Map out spatial effects and log_tau (no spatial in this config)
mapA = list(u_spatial=factor(rep(NA, nAreas)), log_tau=factor(NA))

cat("\n============================================================\n")
cat("FE+NUG(GH Q=10): all FEs inner, outer = log_tauEps only\n")
cat("============================================================\n")

obj = MakeADFun(data=data_gh, parameters=params_A,
                random=c('alpha', 'beta_urban', 'beta_normPop'),
                map=mapA, DLL='modM_MIIDonly_GH', silent=TRUE)

cat("Outer params:", names(obj$par), "\n")
cat("Random effects:", length(obj$env$random), "params\n\n")

t0 = proc.time()[3]
opt = nlminb(obj$par, obj$fn, obj$gr,
             control=list(eval.max=1000, iter.max=500, trace=1))
time_opt = proc.time()[3] - t0

cat(sprintf("\nConvergence: %d (%s)\n", opt$convergence, opt$message))
cat(sprintf("NLL: %.4f\n", opt$objective))
cat(sprintf("Time: %.1f s\n", time_opt))

sd_rep = sdreport(obj)
est_fix = summary(sd_rep, "fixed")
est_ran = summary(sd_rep, "random")
cat("\nFixed (outer) estimates:\n"); print(est_fix)
cat("\nRandom (inner) estimates:\n"); print(est_ran)

pe_f = est_fix[,"Estimate"]; nms_f = rownames(est_fix)
pe_r = est_ran[,"Estimate"]; se_r = est_ran[,"Std. Error"]
nms_r = rownames(est_ran)

sigEps = exp(-0.5 * pe_f[nms_f=="log_tauEps"])

cat(sprintf("\n--- All-FE-inner Results ---\n"))
cat(sprintf("  alpha:        %7.4f (SE %.4f)  truth: -1.2500\n",
            pe_r[nms_r=="alpha"], se_r[nms_r=="alpha"]))
cat(sprintf("  beta_urban:   %7.4f (SE %.4f)  truth:  1.0000\n",
            pe_r[nms_r=="beta_urban"], se_r[nms_r=="beta_urban"]))
cat(sprintf("  beta_normPop: %7.4f (SE %.4f)  truth:  0.5000\n",
            pe_r[nms_r=="beta_normPop"], se_r[nms_r=="beta_normPop"]))
cat(sprintf("  sigmaEps:     %7.4f           truth:  1.2247\n", sigEps))

cat("\n--- Comparison ---\n")
cat(sprintf("  %-30s %8s %8s %8s %10s %8s\n",
            "Method", "normPop", "alpha", "urban", "sigmaEps", "Time"))
cat(sprintf("  %-30s %8.4f %8.4f %8.4f %10.4f %8s\n",
            "Truth", 0.50, -1.25, 1.00, 1.2247, "--"))
cat(sprintf("  %-30s %8.4f %8.4f %8.4f %10.4f %8.1f\n",
            "GH FE+nug all-outer nlminb", 0.6769, -1.5182, 0.8390, 1.2908, 15.8))
cat(sprintf("  %-30s %8.4f %8.4f %8.4f %10.4f %8.1f\n",
            "GH FE+nug all-outer BFGS", 0.6769, -1.5182, 0.8389, 1.2908, 10.4))
cat(sprintf("  %-30s %8.4f %8.4f %8.4f %10.4f %8.1f\n",
            "GH FE+nug FEs-inner",
            pe_r[nms_r=="beta_normPop"],
            pe_r[nms_r=="alpha"],
            pe_r[nms_r=="beta_urban"],
            sigEps, time_opt))
cat("Done.\n")
