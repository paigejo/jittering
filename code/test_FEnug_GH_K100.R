#!/usr/bin/env Rscript
# Test FE+nugget GH model with K=100 (verifying the weight-matrix fix)
# Compiles with safebounds=TRUE to confirm no Eigen bounds errors

library(TMB)
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
source("code/modM_MIIDonly.R")

# ── Load BYM2 simulated data ─────────────────────────────────────
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

# Subset to urban + normPop covariates
allCovNames = colnames(inputsMDM$intPtsMICS$XUrb)
keepIdx = which(allCovNames %in% c("urban", "normPop"))
inputsMDM$intPtsMICS$XUrb = inputsMDM$intPtsMICS$XUrb[, keepIdx, drop=FALSE]
inputsMDM$intPtsMICS$XRur = inputsMDM$intPtsMICS$XRur[, keepIdx, drop=FALSE]

cat("X columns:", colnames(inputsMDM$intPtsMICS$XUrb), "\n")
thisEnv = environment()
list2env(inputsMDM, envir=thisEnv)

lambdaTau = getLambdaPCprec(u=1, alpha=0.1)
lambdaTauEps = getLambdaPCprec(u=1, alpha=0.1)
nAreas = max(c(areaidxlocUrbanMICS, areaidxlocRuralMICS)) + 1

Q = 25
gh = fastGHQuad::gaussHermiteData(Q)

# ── Verify weight matrix dimensions ─────────────────────────────
cat(sprintf("\nWeight matrix shapes:\n"))
cat(sprintf("  wUrban: %d x %d  (expect nObs x K = %d x 100)\n",
            nrow(intPtsMICS$wUrban), ncol(intPtsMICS$wUrban), length(ysUrbMICS)))
cat(sprintf("  wRural: %d x %d  (expect nObs x K = %d x 100)\n",
            nrow(intPtsMICS$wRural), ncol(intPtsMICS$wRural), length(ysRurMICS)))
stopifnot(ncol(intPtsMICS$wUrban) == 100)
stopifnot(ncol(intPtsMICS$wRural) == 100)

# ── TMB Data ────────────────────────────────────────────────────
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

# ── Compile with safebounds=TRUE to verify no Eigen errors ──────
unloadDynlibs()
cat("\nCompiling modM_MIIDonly_GH.cpp with safebounds=TRUE...\n")
compile("code/modM_MIIDonly_GH.cpp")  # default: safebounds=TRUE
dyn.load(dynlib("code/modM_MIIDonly_GH"))

# ── True values ──────────────────────────────────────────────────
trueAlpha = -1.25; trueUrban = 1.00; trueNormPop = 0.50
trueSigmaEps = sqrt(1.5)

# ── Initial values ───────────────────────────────────────────────
initUrbP = sum(ysUrbMICS)/sum(nsUrbMICS)
initRurP = sum(ysRurMICS)/sum(nsRurMICS)
initAlpha = log(initRurP/(1-initRurP))
initBeta1 = log(initUrbP/(1-initUrbP)) - initAlpha

params_A = list(
  log_tau=0, log_tauEps=0,
  alpha=initAlpha, beta_urban=initBeta1, beta_normPop=0,
  u_spatial=rep(0, nAreas)
)
mapA = list(u_spatial=factor(rep(NA, nAreas)), log_tau=factor(NA))

cat("\nBuilding AD tape (FE+nugget, K=100, Q=25)...\n")
obj_A = MakeADFun(data=data_gh, parameters=params_A,
                  map=mapA, DLL='modM_MIIDonly_GH', silent=TRUE)
cat("MakeADFun succeeded!\n")
cat("Parameters to optimize:", names(obj_A$par), "\n")

cat("\nOptimizing...\n"); flush.console()
t0 = proc.time()[3]
opt_A = nlminb(obj_A$par, obj_A$fn, obj_A$gr,
               control=list(eval.max=1000, iter.max=500, trace=1))
time_A = proc.time()[3] - t0

cat(sprintf("\nConvergence: %d (%s)\n", opt_A$convergence, opt_A$message))
cat(sprintf("NLL: %.4f\n", opt_A$objective))
cat(sprintf("Time: %.1f s\n", time_A))

sd_A = sdreport(obj_A)
est_A = summary(sd_A, "fixed")
cat("\nFixed parameter estimates:\n")
print(est_A)

pe = est_A[,"Estimate"]; se = est_A[,"Std. Error"]
nms = rownames(est_A)

cat(sprintf("\n=== FE+nugget GH (K=100, Q=25) Results ===\n"))
cat(sprintf("  alpha:        %7.4f (SE %.4f)  truth: %7.4f  |bias|/SE: %.2f\n",
            pe[nms=="alpha"], se[nms=="alpha"], trueAlpha,
            abs(pe[nms=="alpha"]-trueAlpha)/se[nms=="alpha"]))
cat(sprintf("  beta_urban:   %7.4f (SE %.4f)  truth: %7.4f  |bias|/SE: %.2f\n",
            pe[nms=="beta_urban"], se[nms=="beta_urban"], trueUrban,
            abs(pe[nms=="beta_urban"]-trueUrban)/se[nms=="beta_urban"]))
cat(sprintf("  beta_normPop: %7.4f (SE %.4f)  truth: %7.4f  |bias|/SE: %.2f\n",
            pe[nms=="beta_normPop"], se[nms=="beta_normPop"], trueNormPop,
            abs(pe[nms=="beta_normPop"]-trueNormPop)/se[nms=="beta_normPop"]))
sigmaEps_est = exp(-0.5 * pe[nms=="log_tauEps"])
cat(sprintf("  sigmaEps:     %7.4f              truth: %7.4f\n", sigmaEps_est, trueSigmaEps))
cat(sprintf("\nDone.\n"))
