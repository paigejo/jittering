#!/usr/bin/env Rscript
# Config A FE+nug GH Q=10 with BFGS optimizer (instead of nlminb)

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
mapA = list(u_spatial=factor(rep(NA, nAreas)), log_tau=factor(NA))

obj_A = MakeADFun(data=data_gh, parameters=params_A,
                  map=mapA, DLL='modM_MIIDonly_GH', silent=TRUE)

cat("Starting BFGS with gradient...\n")
t0 = proc.time()[3]
opt_bfgs = optim(obj_A$par, obj_A$fn, obj_A$gr, method="BFGS",
                 control=list(trace=1, REPORT=1, maxit=500))
time_bfgs = proc.time()[3] - t0

cat(sprintf("\nConvergence: %d\nNLL: %.4f\nTime: %.1f s\n",
            opt_bfgs$convergence, opt_bfgs$value, time_bfgs))
cat("Estimates:\n"); print(opt_bfgs$par)
sigEps = exp(-0.5 * opt_bfgs$par[["log_tauEps"]])

cat(sprintf("\n--- BFGS Results ---\n"))
cat(sprintf("  normPop = %.4f   alpha = %.4f   urban = %.4f   sigmaEps = %.4f\n",
            opt_bfgs$par[["beta_normPop"]], opt_bfgs$par[["alpha"]],
            opt_bfgs$par[["beta_urban"]], sigEps))
cat(sprintf("  Time: %.1f s\n", time_bfgs))

cat("\n--- Comparison ---\n")
cat(sprintf("  %-20s %8s %8s %8s %10s %8s\n", "Method", "normPop", "alpha", "urban", "sigmaEps", "Time"))
cat(sprintf("  %-20s %8.4f %8.4f %8.4f %10.4f %8s\n", "Truth", 0.50, -1.25, 1.00, 1.2247, "--"))
cat(sprintf("  %-20s %8.4f %8.4f %8.4f %10.4f %8.1f\n", "nlminb Q=10", 0.6769, -1.5182, 0.8390, 1.2908, 15.8))
cat(sprintf("  %-20s %8.4f %8.4f %8.4f %10.4f %8.1f\n", "BFGS Q=10",
            opt_bfgs$par[["beta_normPop"]], opt_bfgs$par[["alpha"]],
            opt_bfgs$par[["beta_urban"]], sigEps, time_bfgs))
cat("Done.\n")
