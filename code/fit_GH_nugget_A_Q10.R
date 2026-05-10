#!/usr/bin/env Rscript
# Config (A) ONLY: FE + nugget(GH) with Q=10 nodes
# Compare Q=10 accuracy against Q=25 result

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
if(!("Stratum" %in% names(edMICS))) {
  edMICS$Stratum = adm2ToStratumMICS(edMICS$subarea)
}

inputsMDM = makeInputsMDM(datDHS=ed, datMICS=edMICS,
                          KMICS=100, admMICS=admFinal, adm2DHS=adm2Full)

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

# ── Gauss-Hermite nodes — Q=10 ──────────────────────────────────
Q = 10
gh = fastGHQuad::gaussHermiteData(Q)
gh_nodes   = gh$x
gh_weights = gh$w
cat(sprintf("Using %d Gauss-Hermite quadrature nodes\n", Q))

# ── True values ──────────────────────────────────────────────────
trueAlpha = -1.25; trueUrban = 1.00; trueNormPop = 0.50
trueSigmaEps = sqrt(1.5)
trueLogTauEps = -2*log(trueSigmaEps)

# ── Initial values ───────────────────────────────────────────────
initUrbP = sum(ysUrbMICS)/sum(nsUrbMICS)
initRurP = sum(ysRurMICS)/sum(nsRurMICS)
initAlpha = log(initRurP/(1-initRurP))
initBeta1 = log(initUrbP/(1-initUrbP)) - initAlpha

# ── TMB Data ─────────────────────────────────────────────────────
data_gh = list(
  y_iUrbanMICS=ysUrbMICS, y_iRuralMICS=ysRurMICS,
  n_iUrbanMICS=nsUrbMICS, n_iRuralMICS=nsRurMICS,
  areaidxlocUrbanMICS=areaidxlocUrbanMICS, areaidxlocRuralMICS=areaidxlocRuralMICS,
  X_betaUrbanMICS=intPtsMICS$XUrb, X_betaRuralMICS=intPtsMICS$XRur,
  wUrbanMICS=intPtsMICS$wUrban, wRuralMICS=intPtsMICS$wRural,
  nAreas=nAreas,
  gh_nodes=gh_nodes, gh_weights=gh_weights,
  alpha_pri=c(0, 100^2), beta_pri=c(0, sqrt(1000)),
  lambdaTau=lambdaTau, lambdaTauEps=lambdaTauEps
)

# ── Compile ──────────────────────────────────────────────────────
unloadDynlibs()
cat("\nCompiling modM_MIIDonly_GH.cpp...\n")
compile("code/modM_MIIDonly_GH.cpp", framework="TMBad", safebounds=FALSE)
dyn.load(dynlib("code/modM_MIIDonly_GH"))

# ── Config A: FE+nug, Q=10 ──────────────────────────────────────
cat("\n============================================================\n")
cat("(A) FE+NUGGET(GH Q=10): pure optimization (no random effects)\n")
cat("============================================================\n")

params_A = list(
  log_tau=0, log_tauEps=0,
  alpha=initAlpha, beta_urban=initBeta1, beta_normPop=0,
  u_spatial=rep(0, nAreas)
)
mapA = list(u_spatial=factor(rep(NA, nAreas)), log_tau=factor(NA))

obj_A = MakeADFun(data=data_gh, parameters=params_A,
                  map=mapA, DLL='modM_MIIDonly_GH', silent=TRUE)

cat("Parameters to optimize:", names(obj_A$par), "\n")
cat("Number of params:", length(obj_A$par), "\n\n")

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
sigEps_est = exp(-0.5*pe[nms=="log_tauEps"])

cat(sprintf("\n--- FE+nugget(GH Q=10) Results ---\n"))
cat(sprintf("  alpha:       %7.4f (SE %.4f)  truth: %7.4f\n",
            pe[nms=="alpha"], se[nms=="alpha"], trueAlpha))
cat(sprintf("  beta_urban:  %7.4f (SE %.4f)  truth: %7.4f\n",
            pe[nms=="beta_urban"], se[nms=="beta_urban"], trueUrban))
cat(sprintf("  beta_normPop:%7.4f (SE %.4f)  truth: %7.4f\n",
            pe[nms=="beta_normPop"], se[nms=="beta_normPop"], trueNormPop))
cat(sprintf("  log_tauEps:  %7.4f (SE %.4f)  truth: %7.4f\n",
            pe[nms=="log_tauEps"], se[nms=="log_tauEps"], trueLogTauEps))
cat(sprintf("  sigmaEps:    %7.4f  (truth: %.4f)\n", sigEps_est, trueSigmaEps))

# ── Comparison with Q=25 ────────────────────────────────────────
cat("\n\n════════════════════════════════════════════════════════════\n")
cat("COMPARISON: Q=10 vs Q=25 and other methods\n")
cat("════════════════════════════════════════════════════════════\n\n")

cat(sprintf("%-35s %8s %8s %8s %10s %8s\n",
            "Method", "normPop", "alpha", "urban", "sigmaEps", "Time(s)"))
cat(paste0(rep("-", 80), collapse=""), "\n")
cat(sprintf("%-35s %8.3f %8.3f %8.3f %10.4f %8s\n",
            "Truth", 0.50, -1.25, 1.00, trueSigmaEps, "--"))
cat(sprintf("%-35s %8.3f %8.3f %8.3f %10s %8.0f\n",
            "Laplace all-random BFGS", 0.130, -0.940, 1.300, "--", 103))
cat(sprintf("%-35s %8.3f %8.3f %8.3f %10s %8.0f\n",
            "Laplace normPop-outer NM", 0.890, -1.680, 0.630, "--", 239))
cat(sprintf("%-35s %8.3f %8.3f %8.3f %10s %8.0f\n",
            "MCMC laplace=FALSE", 0.680, -1.520, 0.860, "--", 24127))
cat(sprintf("%-35s %8.4f %8.4f %8.4f %10.4f %8.1f\n",
            "GH FE+nug Q=25", 0.6854, -1.5267, 0.8531, 1.2927, 60.7))
cat(sprintf("%-35s %8.4f %8.4f %8.4f %10.4f %8.1f\n",
            "GH FE+nug Q=10",
            pe[nms=="beta_normPop"], pe[nms=="alpha"],
            pe[nms=="beta_urban"], sigEps_est, time_A))

# ── Save ─────────────────────────────────────────────────────────
resultA_Q10 = list(opt=opt_A, sd=sd_A, time=time_A, Q=Q)
save(resultA_Q10, file="savedOutput/fitGH_FEnug_Q10.RData")
cat("\nSaved to savedOutput/fitGH_FEnug_Q10.RData\n")
cat("Done.\n")
