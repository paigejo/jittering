#!/usr/bin/env Rscript
# Config (B) only: IID + nugget(GH) with fewer GH nodes (Q=10)
# to avoid the memory crash from Q=25

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
nUrb = length(ysUrbMICS)
nRur = length(ysRurMICS)

# ── Gauss-Hermite nodes — Q=10 to fit in memory ─────────────────
Q = 10
gh = fastGHQuad::gaussHermiteData(Q)
gh_nodes   = gh$x
gh_weights = gh$w
cat(sprintf("Using %d Gauss-Hermite quadrature nodes\n", Q))

# ── True values ──────────────────────────────────────────────────
trueAlpha = -1.25; trueUrban = 1.00; trueNormPop = 0.50
trueSigmaU = sqrt(0.5); trueSigmaEps = sqrt(1.5)
trueLogTau = -2*log(trueSigmaU); trueLogTauEps = -2*log(trueSigmaEps)

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
TMB::config(tmbad.sparse_hessian_compress = 1)

# ═══════════════════════════════════════════════════════════════
# Also run Config (A) FE+nug with Q=10 for comparison with Q=25
# ═══════════════════════════════════════════════════════════════
cat("\n============================================================\n")
cat("(A-10) FE+NUGGET(GH Q=10): pure optimization\n")
cat("============================================================\n")

params_A = list(
  log_tau=0, log_tauEps=0,
  alpha=initAlpha, beta_urban=initBeta1, beta_normPop=0,
  u_spatial=rep(0, nAreas)
)
mapA = list(u_spatial=factor(rep(NA, nAreas)), log_tau=factor(NA))

obj_A10 = MakeADFun(data=data_gh, parameters=params_A,
                    map=mapA, DLL='modM_MIIDonly_GH', silent=TRUE)

t0 = proc.time()[3]
opt_A10 = nlminb(obj_A10$par, obj_A10$fn, obj_A10$gr,
                 control=list(eval.max=1000, iter.max=500, trace=1))
time_A10 = proc.time()[3] - t0

cat(sprintf("\nConvergence: %d (%s)\n", opt_A10$convergence, opt_A10$message))
cat(sprintf("NLL: %.4f  Time: %.1f s\n", opt_A10$objective, time_A10))
sd_A10 = sdreport(obj_A10)
est_A10 = summary(sd_A10, "fixed")
cat("\nFixed parameter estimates (FE+nug Q=10):\n")
print(est_A10)

pe10 = est_A10[,"Estimate"]; nms10 = rownames(est_A10)
sigEps_A10 = exp(-0.5*pe10[nms10=="log_tauEps"])
cat(sprintf("  sigmaEps: %.4f (truth: %.4f)\n", sigEps_A10, trueSigmaEps))

# Clean up A10 tape to free memory before B
rm(obj_A10); gc()

# ═══════════════════════════════════════════════════════════════
# Config (B): IID + nugget(GH Q=10)
# ═══════════════════════════════════════════════════════════════
cat("\n============================================================\n")
cat("(B) IID+NUGGET(GH Q=10): u_spatial + alpha + beta_urban inner\n")
cat("============================================================\n")

params_B = list(
  log_tau=0, log_tauEps=0,
  alpha=initAlpha, beta_urban=initBeta1, beta_normPop=0,
  u_spatial=rep(0, nAreas)
)

cat("Building AD tape...\n"); flush.console()
obj_B = MakeADFun(data=data_gh, parameters=params_B,
                  random=c('alpha', 'beta_urban', 'u_spatial'),
                  DLL='modM_MIIDonly_GH', silent=TRUE)

cat("Outer params:", names(obj_B$par), "\n")
cat("Random effects:", sum(obj_B$env$random), "params\n\n")

t0 = proc.time()[3]
opt_B = nlminb(obj_B$par, obj_B$fn, obj_B$gr,
               control=list(eval.max=1000, iter.max=500, trace=1))
time_B = proc.time()[3] - t0

cat(sprintf("\nConvergence: %d (%s)\n", opt_B$convergence, opt_B$message))
cat(sprintf("NLL: %.4f\n", opt_B$objective))
cat(sprintf("Time: %.1f s\n", time_B))

sd_B = sdreport(obj_B)
est_B_fix = summary(sd_B, "fixed")
est_B_ran = summary(sd_B, "random")
cat("\nFixed parameter estimates:\n"); print(est_B_fix)

pe_f = est_B_fix[,"Estimate"]; se_f = est_B_fix[,"Std. Error"]
nms_f = rownames(est_B_fix)
pe_r = est_B_ran[,"Estimate"]; se_r = est_B_ran[,"Std. Error"]
nms_r = rownames(est_B_ran)

cat(sprintf("\n--- IID+nugget(GH Q=10) Results ---\n"))
cat(sprintf("  alpha:       %7.4f (SE %.4f)  truth: %7.4f  bias: %+.4f\n",
            pe_r[nms_r=="alpha"][1], se_r[nms_r=="alpha"][1], trueAlpha,
            pe_r[nms_r=="alpha"][1]-trueAlpha))
cat(sprintf("  beta_urban:  %7.4f (SE %.4f)  truth: %7.4f  bias: %+.4f\n",
            pe_r[nms_r=="beta_urban"][1], se_r[nms_r=="beta_urban"][1], trueUrban,
            pe_r[nms_r=="beta_urban"][1]-trueUrban))
cat(sprintf("  beta_normPop:%7.4f (SE %.4f)  truth: %7.4f  bias: %+.4f\n",
            pe_f[nms_f=="beta_normPop"], se_f[nms_f=="beta_normPop"], trueNormPop,
            pe_f[nms_f=="beta_normPop"]-trueNormPop))
cat(sprintf("  log_tau:     %7.4f (SE %.4f)  truth: %7.4f\n",
            pe_f[nms_f=="log_tau"], se_f[nms_f=="log_tau"], trueLogTau))
cat(sprintf("  log_tauEps:  %7.4f (SE %.4f)  truth: %7.4f\n",
            pe_f[nms_f=="log_tauEps"], se_f[nms_f=="log_tauEps"], trueLogTauEps))
sigU_est = exp(-0.5*pe_f[nms_f=="log_tau"])
sigEps_B = exp(-0.5*pe_f[nms_f=="log_tauEps"])
cat(sprintf("  sigmaU:      %7.4f  (truth: %.4f)\n", sigU_est, trueSigmaU))
cat(sprintf("  sigmaEps:    %7.4f  (truth: %.4f)\n", sigEps_B, trueSigmaEps))

resultB = list(opt=opt_B, sd=sd_B, time=time_B)
resultA10 = list(opt=opt_A10, sd=sd_A10, time=time_A10)

# ── Comparison ───────────────────────────────────────────────────
cat("\n\n════════════════════════════════════════════════════════════\n")
cat("COMPARISON TABLE\n")
cat("════════════════════════════════════════════════════════════\n\n")

cat(sprintf("%-35s %8s %8s %8s %10s %10s\n",
            "Method", "normPop", "alpha", "urban", "sigmaEps", "Time(s)"))
cat(paste0(rep("-", 85), collapse=""), "\n")
cat(sprintf("%-35s %8.2f %8.2f %8.2f %10.4f %10s\n",
            "Truth", 0.50, -1.25, 1.00, trueSigmaEps, "--"))
cat(sprintf("%-35s %8.2f %8.2f %8.2f %10s %10.0f\n",
            "Laplace: all random BFGS", 0.13, -0.94, 1.30, "--", 103))
cat(sprintf("%-35s %8.2f %8.2f %8.2f %10s %10.0f\n",
            "Laplace: normPop outer NM", 0.89, -1.68, 0.63, "--", 239))
cat(sprintf("%-35s %8.2f %8.2f %8.2f %10s %10.0f\n",
            "MCMC laplace=FALSE", 0.68, -1.52, 0.86, "--", 24127))
cat(sprintf("%-35s %8.4f %8.4f %8.4f %10.4f %10.1f\n",
            "GH FE+nug Q=25", 0.6854, -1.5267, 0.8531, 1.2927, 60.7))
cat(sprintf("%-35s %8.4f %8.4f %8.4f %10.4f %10.1f\n",
            "GH FE+nug Q=10",
            pe10[nms10=="beta_normPop"], pe10[nms10=="alpha"],
            pe10[nms10=="beta_urban"], sigEps_A10, time_A10))
cat(sprintf("%-35s %8.4f %8.4f %8.4f %10.4f %10.1f\n",
            "GH IID+nug Q=10",
            pe_f[nms_f=="beta_normPop"],
            pe_r[nms_r=="alpha"][1],
            pe_r[nms_r=="beta_urban"][1],
            sigEps_B, time_B))

# ── Save ─────────────────────────────────────────────────────────
save(resultA10, resultB, file="savedOutput/fitGH_IIDnug_Q10.RData")
cat("\nResults saved to savedOutput/fitGH_IIDnug_Q10.RData\n")
cat("\nDone.\n")
