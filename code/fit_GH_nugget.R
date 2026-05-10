#!/usr/bin/env Rscript
# Fit FE+nugget model using Gauss-Hermite quadrature for nuggets
# Nuggets are integrated out in C++ — no random effects for FE-only model
# Compare: (A) FE-only (pure 4-param optimization)
#          (B) IID + nugget (u_spatial + alpha + beta_urban in inner Newton)

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

# Subset to urban + normPop covariates
allCovNames = colnames(inputsMDM$intPtsMICS$XUrb)
keepIdx = which(allCovNames %in% c("urban", "normPop"))
inputsMDM$intPtsMICS$XUrb = inputsMDM$intPtsMICS$XUrb[, keepIdx, drop=FALSE]
inputsMDM$intPtsMICS$XRur = inputsMDM$intPtsMICS$XRur[, keepIdx, drop=FALSE]

cat("X columns:", colnames(inputsMDM$intPtsMICS$XUrb), "\n")
stopifnot(colnames(inputsMDM$intPtsMICS$XUrb)[1] == "urban")
stopifnot(colnames(inputsMDM$intPtsMICS$XUrb)[2] == "normPop")

thisEnv = environment()
list2env(inputsMDM, envir=thisEnv)

lambdaTau = getLambdaPCprec(u=1, alpha=0.1)
lambdaTauEps = getLambdaPCprec(u=1, alpha=0.1)
nAreas = max(c(areaidxlocUrbanMICS, areaidxlocRuralMICS)) + 1
nUrb = length(ysUrbMICS)
nRur = length(ysRurMICS)

# ── Gauss-Hermite nodes ─────────────────────────────────────────
# Using fastGHQuad if available, otherwise manual computation
Q = 25  # number of GH nodes
if(requireNamespace("fastGHQuad", quietly=TRUE)) {
  gh = fastGHQuad::gaussHermiteData(Q)
  gh_nodes   = gh$x
  gh_weights = gh$w
} else {
  # Fallback: use statmod
  if(requireNamespace("statmod", quietly=TRUE)) {
    gh = statmod::gauss.quad(Q, kind="hermite")
    gh_nodes   = gh$nodes
    gh_weights = gh$weights
  } else {
    stop("Need fastGHQuad or statmod package for GH nodes. Install with: install.packages('fastGHQuad')")
  }
}
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

# ── TMB Data (shared) ───────────────────────────────────────────
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

# ── Compile GH template ─────────────────────────────────────────
unloadDynlibs()
cat("\nCompiling modM_MIIDonly_GH.cpp...\n")
compile("code/modM_MIIDonly_GH.cpp", framework="TMBad", safebounds=FALSE)
dyn.load(dynlib("code/modM_MIIDonly_GH"))
TMB::config(tmbad.sparse_hessian_compress = 1)

# ════════════════════════════════════════════════════════════════
# (A) FE + nugget(GH): zero random effects, pure optimization
# ════════════════════════════════════════════════════════════════
cat("\n============================================================\n")
cat("(A) FE+NUGGET(GH): pure optimization (no random effects)\n")
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
cat("Number of params:", length(obj_A$par), "\n")
cat("(No random effects — nuggets integrated out by GH)\n\n")

t0 = proc.time()[3]
opt_A = nlminb(obj_A$par, obj_A$fn, obj_A$gr,
               control=list(eval.max=1000, iter.max=500, trace=1))
time_A = proc.time()[3] - t0

cat(sprintf("\nConvergence: %d (%s)\n", opt_A$convergence, opt_A$message))
cat(sprintf("NLL: %.4f\n", opt_A$objective))
cat(sprintf("Time: %.1f s\n", time_A))

# Get SEs via sdreport
sd_A = sdreport(obj_A)
est_A = summary(sd_A, "fixed")
cat("\nFixed parameter estimates:\n")
print(est_A)

# Extract named estimates
pe = est_A[,"Estimate"]; se = est_A[,"Std. Error"]
nms = rownames(est_A)

cat(sprintf("\n--- FE+nugget(GH) Results ---\n"))
cat(sprintf("  alpha:       %7.4f (SE %.4f)  truth: %7.4f  bias: %+.4f\n",
            pe[nms=="alpha"], se[nms=="alpha"], trueAlpha, pe[nms=="alpha"]-trueAlpha))
cat(sprintf("  beta_urban:  %7.4f (SE %.4f)  truth: %7.4f  bias: %+.4f\n",
            pe[nms=="beta_urban"], se[nms=="beta_urban"], trueUrban, pe[nms=="beta_urban"]-trueUrban))
cat(sprintf("  beta_normPop:%7.4f (SE %.4f)  truth: %7.4f  bias: %+.4f\n",
            pe[nms=="beta_normPop"], se[nms=="beta_normPop"], trueNormPop, pe[nms=="beta_normPop"]-trueNormPop))
cat(sprintf("  log_tauEps:  %7.4f (SE %.4f)  truth: %7.4f\n",
            pe[nms=="log_tauEps"], se[nms=="log_tauEps"], trueLogTauEps))
sigEps_est = exp(-0.5*pe[nms=="log_tauEps"])
cat(sprintf("  sigmaEps:    %7.4f  (truth: %.4f)\n", sigEps_est, trueSigmaEps))

resultA = list(opt=opt_A, sd=sd_A, time=time_A)

# ════════════════════════════════════════════════════════════════
# (B) IID + nugget(GH): u_spatial in inner Newton
# ════════════════════════════════════════════════════════════════
cat("\n============================================================\n")
cat("(B) IID+NUGGET(GH): u_spatial + alpha + beta_urban inner\n")
cat("============================================================\n")

params_B = list(
  log_tau=0, log_tauEps=0,
  alpha=initAlpha, beta_urban=initBeta1, beta_normPop=0,
  u_spatial=rep(0, nAreas)
)

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

# sdreport
sd_B = sdreport(obj_B)
est_B_fix = summary(sd_B, "fixed")
est_B_ran = summary(sd_B, "random")
cat("\nFixed parameter estimates:\n"); print(est_B_fix)
cat("\nRandom effect summaries (first few):\n")
rn = rownames(est_B_ran)
for(nm in unique(rn)) {
  idx = which(rn == nm)
  cat(sprintf("  %s: %d params, mean=%.4f\n", nm, length(idx), mean(est_B_ran[idx,"Estimate"])))
}

# Extract
pe_f = est_B_fix[,"Estimate"]; se_f = est_B_fix[,"Std. Error"]
nms_f = rownames(est_B_fix)
pe_r = est_B_ran[,"Estimate"]; se_r = est_B_ran[,"Std. Error"]
nms_r = rownames(est_B_ran)

cat(sprintf("\n--- IID+nugget(GH) Results ---\n"))
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

# ════════════════════════════════════════════════════════════════
# Comparison table
# ════════════════════════════════════════════════════════════════
cat("\n\n════════════════════════════════════════════════════════════\n")
cat("COMPARISON TABLE (BYM2 simulated data)\n")
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

# FE+nug(GH)
cat(sprintf("%-35s %8.4f %8.4f %8.4f %10.4f %10.1f\n",
            "GH: FE+nug (pure optim)",
            pe[nms=="beta_normPop"], pe[nms=="alpha"], pe[nms=="beta_urban"],
            sigEps_est, time_A))

# IID+nug(GH)
cat(sprintf("%-35s %8.4f %8.4f %8.4f %10.4f %10.1f\n",
            "GH: IID+nug (u_spatial inner)",
            pe_f[nms_f=="beta_normPop"],
            pe_r[nms_r=="alpha"][1],
            pe_r[nms_r=="beta_urban"][1],
            sigEps_B, time_B))

cat("\nUncertainty (SEs):\n")
cat(sprintf("  GH FE+nug:    normPop SE=%.4f  alpha SE=%.4f  urban SE=%.4f\n",
            se[nms=="beta_normPop"], se[nms=="alpha"], se[nms=="beta_urban"]))
cat(sprintf("  GH IID+nug:   normPop SE=%.4f  alpha SE=%.4f  urban SE=%.4f\n",
            se_f[nms_f=="beta_normPop"],
            se_r[nms_r=="alpha"][1],
            se_r[nms_r=="beta_urban"][1]))
cat(sprintf("  MCMC:         normPop SD=0.09   alpha SD=0.09   urban SD=0.13\n"))

# ── Save ─────────────────────────────────────────────────────────
save(resultA, resultB, file="savedOutput/fitGH_FEnug_IIDnug.RData")
cat("\nResults saved to savedOutput/fitGH_FEnug_IIDnug.RData\n")

cat("\nDone.\n")
