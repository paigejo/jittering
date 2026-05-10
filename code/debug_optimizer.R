#!/usr/bin/env Rscript
# Test: alpha/beta as outer (fixed) params + Nelder-Mead + init at truth
# Tests options 1+2+4 simultaneously to see if the normPop estimation improves.

library(TMB)
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
source("code/modM_MIIDonly.R")

# Load BYM2 sim data
load("savedOutput/simStudy1/simPopsSurveys_BYM2.RData")
edMICS = surveysMICS[[1]]

# Name conversions
nameTab = rbind(c("N", "ns"), c("N", "n"), c("Z", "ys"), c("Z", "y"))
for(i in 1:nrow(nameTab)) {
  fromN = nameTab[i,1]; toN = nameTab[i,2]
  if(!(toN %in% names(edMICS))) edMICS[[toN]] = edMICS[[fromN]]
  if(!(toN %in% names(ed))) ed[[toN]] = ed[[fromN]]
}
if(!("Stratum" %in% names(edMICS))) {
  edMICS$Stratum = adm2ToStratumMICS(edMICS$subarea)
}

# Make inputs
inputsMDM = makeInputsMDM(datDHS=ed, datMICS=edMICS,
                          KMICS=100, admMICS=admFinal, adm2DHS=adm2Full)

# Subset covariates
allCovNames = colnames(inputsMDM$intPtsMICS$XUrb)
keepIdx = which(allCovNames %in% c("urban", "normPop"))
inputsMDM$intPtsMICS$XUrb = inputsMDM$intPtsMICS$XUrb[, keepIdx, drop=FALSE]
inputsMDM$intPtsMICS$XRur = inputsMDM$intPtsMICS$XRur[, keepIdx, drop=FALSE]

thisEnv = environment()
list2env(inputsMDM, envir=thisEnv)

lambdaTau = getLambdaPCprec(u=1, alpha=0.1)
lambdaTauEps = getLambdaPCprec(u=1, alpha=0.1)
nAreas = max(c(areaidxlocUrbanMICS, areaidxlocRuralMICS)) + 1

data_full = list(
  y_iUrbanMICS=ysUrbMICS,
  y_iRuralMICS=ysRurMICS,
  n_iUrbanMICS=nsUrbMICS,
  n_iRuralMICS=nsRurMICS,
  areaidxlocUrbanMICS=areaidxlocUrbanMICS,
  areaidxlocRuralMICS=areaidxlocRuralMICS,
  X_betaUrbanMICS=intPtsMICS$XUrb,
  X_betaRuralMICS=intPtsMICS$XRur,
  wUrbanMICS=intPtsMICS$wUrban,
  wRuralMICS=intPtsMICS$wRural,
  nAreas=nAreas,
  alpha_pri=c(0, 100^2),
  beta_pri=c(0, sqrt(1000)),
  lambdaTau=lambdaTau,
  lambdaTauEps=lambdaTauEps
)

nUrb = length(data_full$y_iUrbanMICS)
nRur = length(data_full$y_iRuralMICS)
nBeta = ncol(data_full$X_betaUrbanMICS)

# True parameter values
trueAlpha = -1.25
trueUrban = 1.00
trueNormPop = 0.50
trueSigmaU = sqrt(0.5)  # log_tau = -2*log(sigmaU) = -2*log(sqrt(0.5)) = log(2) ~ 0.693
trueSigmaEps = sqrt(1.5) # log_tauEps = -2*log(sigmaEps) = -2*log(sqrt(1.5)) = -log(1.5) ~ -0.405
trueLogTau = -2*log(trueSigmaU)
trueLogTauEps = -2*log(trueSigmaEps)

cat("True parameters:\n")
cat("  alpha:", trueAlpha, "\n")
cat("  beta_urban:", trueUrban, "\n")
cat("  beta_normPop:", trueNormPop, "\n")
cat("  log_tau:", trueLogTau, "(sigmaU =", trueSigmaU, ")\n")
cat("  log_tauEps:", trueLogTauEps, "(sigmaEps =", trueSigmaEps, ")\n\n")

unloadDynlibs()
if(!any(file.exists(paste0("code/modM_MIIDonly", c(".o", ".so", ".dll"))))) {
  compile("code/modM_MIIDonly.cpp", framework="TMBad", safebounds=FALSE)
}
dyn.load(dynlib("code/modM_MIIDonly"))
TMB::config(tmbad.sparse_hessian_compress = 1)

# ============================================================
# Helper to extract estimates
# ============================================================
extractResults = function(obj, opt, label, fixedBeta=FALSE) {
  cat(sprintf("\n--- %s ---\n", label))
  cat("Convergence:", opt$convergence, "\n")
  cat("NLL:", opt$value, "\n")
  
  pn = names(opt$par)
  
  if(fixedBeta) {
    # alpha, beta are outer (fixed) params
    alpha_est = as.numeric(opt$par[pn == "alpha"])
    beta_est = as.numeric(opt$par[pn == "beta"])
    log_tau_est = as.numeric(opt$par[pn == "log_tau"])
    log_tauEps_est = as.numeric(opt$par[pn == "log_tauEps"])
  } else {
    # alpha, beta are random — get from last.par.best
    log_tau_est = as.numeric(opt$par[pn == "log_tau"])
    log_tauEps_est = as.numeric(opt$par[pn == "log_tauEps"])
    lp = obj$env$last.par.best
    rn = names(lp)
    alpha_est = as.numeric(lp[rn == "alpha"])
    beta_est = as.numeric(lp[rn == "beta"])
  }
  
  sigmaU = exp(-0.5 * log_tau_est)
  sigmaEps = exp(-0.5 * log_tauEps_est)
  
  cat(sprintf("  alpha:       %7.4f  (truth: %7.4f, bias: %+.4f)\n", 
              alpha_est, trueAlpha, alpha_est - trueAlpha))
  cat(sprintf("  beta_urban:  %7.4f  (truth: %7.4f, bias: %+.4f)\n",
              beta_est[1], trueUrban, beta_est[1] - trueUrban))
  cat(sprintf("  beta_normPop:%7.4f  (truth: %7.4f, bias: %+.4f)\n",
              beta_est[2], trueNormPop, beta_est[2] - trueNormPop))
  cat(sprintf("  sigmaU:      %7.4f  (truth: %7.4f)\n", sigmaU, trueSigmaU))
  cat(sprintf("  sigmaEps:    %7.4f  (truth: %7.4f)\n", sigmaEps, trueSigmaEps))
  
  list(alpha=alpha_est, beta=beta_est, log_tau=log_tau_est, 
       log_tauEps=log_tauEps_est, nll=opt$value, conv=opt$convergence)
}

# ============================================================
# (A) BASELINE: alpha/beta random, BFGS, default init
# ============================================================
cat("============================================================\n")
cat("(A) BASELINE: alpha/beta RANDOM, BFGS, default init\n")
cat("============================================================\n")

initUrbP = sum(data_full$y_iUrbanMICS)/sum(data_full$n_iUrbanMICS)
initRurP = sum(data_full$y_iRuralMICS)/sum(data_full$n_iRuralMICS)
initAlpha = log(initRurP/(1-initRurP))
initBeta1 = log(initUrbP/(1-initUrbP)) - initAlpha

tmb_params_A = list(
  log_tau = 0, log_tauEps = 0,
  alpha = initAlpha,
  beta = c(initBeta1, 0),
  u_spatial = rep(0, nAreas),
  nuggetUrbMICS = rep(0, nUrb),
  nuggetRurMICS = rep(0, nRur)
)

obj_A = MakeADFun(data=data_full, parameters=tmb_params_A,
                  random=c('alpha', 'beta', 'u_spatial', 'nuggetUrbMICS', 'nuggetRurMICS'),
                  hessian=TRUE, DLL='modM_MIIDonly', silent=TRUE)

startA = proc.time()[3]
opt_A = optim(par=obj_A$par, fn=obj_A$fn, gr=obj_A$gr, method="BFGS",
              control=list(reltol=1e-6))
timeA = proc.time()[3] - startA
cat(sprintf("Time: %.1f sec\n", timeA))
res_A = extractResults(obj_A, opt_A, "Baseline (random beta, BFGS, default init)", fixedBeta=FALSE)

# ============================================================
# (B) OPTIONS 1+2+4: alpha/beta FIXED, Nelder-Mead, init at truth
# ============================================================
cat("\n============================================================\n")
cat("(B) OPTIONS 1+2+4: alpha/beta FIXED, Nelder-Mead, init at truth\n")
cat("============================================================\n")

tmb_params_B = list(
  log_tau = trueLogTau, log_tauEps = trueLogTauEps,
  alpha = trueAlpha,
  beta = c(trueUrban, trueNormPop),
  u_spatial = rep(0, nAreas),
  nuggetUrbMICS = rep(0, nUrb),
  nuggetRurMICS = rep(0, nRur)
)

# alpha/beta are outer (fixed), only u_spatial and nuggets are random
obj_B = MakeADFun(data=data_full, parameters=tmb_params_B,
                  random=c('u_spatial', 'nuggetUrbMICS', 'nuggetRurMICS'),
                  hessian=TRUE, DLL='modM_MIIDonly', silent=TRUE)

cat("Outer params:", names(obj_B$par), "\n")
cat("Init values:", obj_B$par, "\n")

startB = proc.time()[3]
opt_B = optim(par=obj_B$par, fn=obj_B$fn, method="Nelder-Mead",
              control=list(reltol=1e-8, maxit=5000))
timeB = proc.time()[3] - startB
cat(sprintf("Time: %.1f sec\n", timeB))
res_B = extractResults(obj_B, opt_B, "Options 1+2+4 (fixed beta, NM, truth init)", fixedBeta=TRUE)

# ============================================================
# (C) OPTION 1 only: alpha/beta FIXED, BFGS, default init
# ============================================================
cat("\n============================================================\n")
cat("(C) OPTION 1 ONLY: alpha/beta FIXED, BFGS, default init\n")
cat("============================================================\n")

tmb_params_C = list(
  log_tau = 0, log_tauEps = 0,
  alpha = initAlpha,
  beta = c(initBeta1, 0),
  u_spatial = rep(0, nAreas),
  nuggetUrbMICS = rep(0, nUrb),
  nuggetRurMICS = rep(0, nRur)
)

obj_C = MakeADFun(data=data_full, parameters=tmb_params_C,
                  random=c('u_spatial', 'nuggetUrbMICS', 'nuggetRurMICS'),
                  hessian=TRUE, DLL='modM_MIIDonly', silent=TRUE)

startC = proc.time()[3]
opt_C = optim(par=obj_C$par, fn=obj_C$fn, gr=obj_C$gr, method="BFGS",
              control=list(reltol=1e-6))
timeC = proc.time()[3] - startC
cat(sprintf("Time: %.1f sec\n", timeC))
res_C = extractResults(obj_C, opt_C, "Option 1 only (fixed beta, BFGS, default init)", fixedBeta=TRUE)

# ============================================================
# (D) OPTION 1+2: alpha/beta FIXED, Nelder-Mead, default init
# ============================================================
cat("\n============================================================\n")
cat("(D) OPTIONS 1+2: alpha/beta FIXED, Nelder-Mead, default init\n")
cat("============================================================\n")

obj_D = MakeADFun(data=data_full, parameters=tmb_params_C,
                  random=c('u_spatial', 'nuggetUrbMICS', 'nuggetRurMICS'),
                  hessian=TRUE, DLL='modM_MIIDonly', silent=TRUE)

startD = proc.time()[3]
opt_D = optim(par=obj_D$par, fn=obj_D$fn, method="Nelder-Mead",
              control=list(reltol=1e-8, maxit=5000))
timeD = proc.time()[3] - startD
cat(sprintf("Time: %.1f sec\n", timeD))
res_D = extractResults(obj_D, opt_D, "Options 1+2 (fixed beta, NM, default init)", fixedBeta=TRUE)

# ============================================================
# (E) OPTION 4 only: alpha/beta RANDOM, BFGS, init at truth
# ============================================================
cat("\n============================================================\n")
cat("(E) OPTION 4 ONLY: alpha/beta RANDOM, BFGS, init at truth\n")
cat("============================================================\n")

tmb_params_E = list(
  log_tau = trueLogTau, log_tauEps = trueLogTauEps,
  alpha = trueAlpha,
  beta = c(trueUrban, trueNormPop),
  u_spatial = rep(0, nAreas),
  nuggetUrbMICS = rep(0, nUrb),
  nuggetRurMICS = rep(0, nRur)
)

obj_E = MakeADFun(data=data_full, parameters=tmb_params_E,
                  random=c('alpha', 'beta', 'u_spatial', 'nuggetUrbMICS', 'nuggetRurMICS'),
                  hessian=TRUE, DLL='modM_MIIDonly', silent=TRUE)

startE = proc.time()[3]
opt_E = optim(par=obj_E$par, fn=obj_E$fn, gr=obj_E$gr, method="BFGS",
              control=list(reltol=1e-6))
timeE = proc.time()[3] - startE
cat(sprintf("Time: %.1f sec\n", timeE))
res_E = extractResults(obj_E, opt_E, "Option 4 only (random beta, BFGS, truth init)", fixedBeta=FALSE)

# ============================================================
# Summary table
# ============================================================
cat("\n\n============================================================\n")
cat("SUMMARY\n")
cat("============================================================\n\n")

cat(sprintf("%-45s %8s %8s %8s %10s %5s\n",
            "Configuration", "alpha", "urban", "normPop", "NLL", "conv"))
cat(sprintf("%-45s %8s %8s %8s %10s %5s\n",
            "-------------", "-----", "-----", "-------", "---", "----"))
cat(sprintf("%-45s %8.4f %8.4f %8.4f %10s %5s\n",
            "TRUTH", trueAlpha, trueUrban, trueNormPop, "", ""))

fmtRow = function(label, res) {
  cat(sprintf("%-45s %8.4f %8.4f %8.4f %10.2f %5d\n",
              label, res$alpha, res$beta[1], res$beta[2], res$nll, res$conv))
}

fmtRow("(A) Random beta, BFGS, default init", res_A)
fmtRow("(B) Fixed beta, NM, truth init [1+2+4]", res_B)
fmtRow("(C) Fixed beta, BFGS, default init [1]", res_C)
fmtRow("(D) Fixed beta, NM, default init [1+2]", res_D)
fmtRow("(E) Random beta, BFGS, truth init [4]", res_E)

# ============================================================
# Get SEs for the best configuration
# ============================================================
cat("\n\nGetting SEs for configuration (B)...\n")
tryCatch({
  SD_B = TMB::sdreport(obj_B, getJointPrecision=TRUE,
                       bias.correct=TRUE,
                       bias.correct.control=list(sd=TRUE))
  cat("pdHess:", SD_B$pdHess, "\n")
  
  pf = SD_B$par.fixed
  se_f = sqrt(diag(SD_B$cov.fixed))
  pn = names(pf)
  
  cat("\nParameter estimates with SEs (config B):\n")
  for(i in seq_along(pf)) {
    cat(sprintf("  %s: %.4f (SE %.4f)\n", pn[i], pf[i], se_f[i]))
  }
}, error = function(e) {
  cat("SE computation failed:", conditionMessage(e), "\n")
})

cat("\nGetting SEs for configuration (D)...\n")
tryCatch({
  SD_D = TMB::sdreport(obj_D, getJointPrecision=TRUE,
                       bias.correct=TRUE,
                       bias.correct.control=list(sd=TRUE))
  cat("pdHess:", SD_D$pdHess, "\n")
  
  pf = SD_D$par.fixed
  se_f = sqrt(diag(SD_D$cov.fixed))
  pn = names(pf)
  
  cat("\nParameter estimates with SEs (config D):\n")
  for(i in seq_along(pf)) {
    cat(sprintf("  %s: %.4f (SE %.4f)\n", pn[i], pf[i], se_f[i]))
  }
}, error = function(e) {
  cat("SE computation failed:", conditionMessage(e), "\n")
})

cat("\nDone.\n")
