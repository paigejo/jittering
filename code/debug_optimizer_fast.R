#!/usr/bin/env Rscript
# Lean test: Compare 3 fast configurations using BFGS
# (C) Option 1: fixed beta, BFGS, default init
# (E) Option 4: random beta, BFGS, truth init
# (F) Options 1+4: fixed beta, BFGS, truth init

library(TMB)
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
source("code/modM_MIIDonly.R")

# Load BYM2 sim data
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

thisEnv = environment()
list2env(inputsMDM, envir=thisEnv)

lambdaTau = getLambdaPCprec(u=1, alpha=0.1)
lambdaTauEps = getLambdaPCprec(u=1, alpha=0.1)
nAreas = max(c(areaidxlocUrbanMICS, areaidxlocRuralMICS)) + 1

data_full = list(
  y_iUrbanMICS=ysUrbMICS, y_iRuralMICS=ysRurMICS,
  n_iUrbanMICS=nsUrbMICS, n_iRuralMICS=nsRurMICS,
  areaidxlocUrbanMICS=areaidxlocUrbanMICS, areaidxlocRuralMICS=areaidxlocRuralMICS,
  X_betaUrbanMICS=intPtsMICS$XUrb, X_betaRuralMICS=intPtsMICS$XRur,
  wUrbanMICS=intPtsMICS$wUrban, wRuralMICS=intPtsMICS$wRural,
  nAreas=nAreas,
  alpha_pri=c(0, 100^2), beta_pri=c(0, sqrt(1000)),
  lambdaTau=lambdaTau, lambdaTauEps=lambdaTauEps
)

nUrb = length(data_full$y_iUrbanMICS)
nRur = length(data_full$y_iRuralMICS)

trueAlpha = -1.25; trueUrban = 1.00; trueNormPop = 0.50
trueSigmaU = sqrt(0.5); trueSigmaEps = sqrt(1.5)
trueLogTau = -2*log(trueSigmaU); trueLogTauEps = -2*log(trueSigmaEps)

initUrbP = sum(data_full$y_iUrbanMICS)/sum(data_full$n_iUrbanMICS)
initRurP = sum(data_full$y_iRuralMICS)/sum(data_full$n_iRuralMICS)
initAlpha = log(initRurP/(1-initRurP))
initBeta1 = log(initUrbP/(1-initUrbP)) - initAlpha

unloadDynlibs()
if(!any(file.exists(paste0("code/modM_MIIDonly", c(".o", ".so", ".dll"))))) {
  compile("code/modM_MIIDonly.cpp", framework="TMBad", safebounds=FALSE)
}
dyn.load(dynlib("code/modM_MIIDonly"))
TMB::config(tmbad.sparse_hessian_compress = 1)

extractResults = function(obj, opt, label, fixedBeta=FALSE) {
  cat(sprintf("\n--- %s ---\n", label))
  cat("Convergence:", opt$convergence, "\n")
  cat("NLL:", opt$value, "\n")
  
  pn = names(opt$par)
  if(fixedBeta) {
    alpha_est = as.numeric(opt$par[pn == "alpha"])
    beta_est = as.numeric(opt$par[pn == "beta"])
    log_tau_est = as.numeric(opt$par[pn == "log_tau"])
    log_tauEps_est = as.numeric(opt$par[pn == "log_tauEps"])
  } else {
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
# (C) Option 1: fixed beta, BFGS, default init
# ============================================================
cat("\n============================================================\n")
cat("(C) OPTION 1: alpha/beta FIXED, BFGS, default init\n")
cat("============================================================\n")

tmb_params_C = list(
  log_tau = 0, log_tauEps = 0,
  alpha = initAlpha, beta = c(initBeta1, 0),
  u_spatial = rep(0, nAreas),
  nuggetUrbMICS = rep(0, nUrb), nuggetRurMICS = rep(0, nRur)
)

obj_C = MakeADFun(data=data_full, parameters=tmb_params_C,
                  random=c('u_spatial', 'nuggetUrbMICS', 'nuggetRurMICS'),
                  hessian=TRUE, DLL='modM_MIIDonly', silent=TRUE)

cat("Outer params:", names(obj_C$par), "\n")
startC = proc.time()[3]
opt_C = optim(par=obj_C$par, fn=obj_C$fn, gr=obj_C$gr, method="BFGS",
              control=list(reltol=1e-6, trace=1, REPORT=1))
timeC = proc.time()[3] - startC
cat(sprintf("\nTime: %.1f sec\n", timeC))
res_C = extractResults(obj_C, opt_C, "Option 1 (fixed beta, BFGS, default init)", fixedBeta=TRUE)

# ============================================================
# (E) Option 4: random beta, BFGS, truth init (hyperparams)
# ============================================================
cat("\n============================================================\n")
cat("(E) OPTION 4: alpha/beta RANDOM, BFGS, truth init\n")
cat("============================================================\n")

tmb_params_E = list(
  log_tau = trueLogTau, log_tauEps = trueLogTauEps,
  alpha = trueAlpha, beta = c(trueUrban, trueNormPop),
  u_spatial = rep(0, nAreas),
  nuggetUrbMICS = rep(0, nUrb), nuggetRurMICS = rep(0, nRur)
)

obj_E = MakeADFun(data=data_full, parameters=tmb_params_E,
                  random=c('alpha', 'beta', 'u_spatial', 'nuggetUrbMICS', 'nuggetRurMICS'),
                  hessian=TRUE, DLL='modM_MIIDonly', silent=TRUE)

startE = proc.time()[3]
opt_E = optim(par=obj_E$par, fn=obj_E$fn, gr=obj_E$gr, method="BFGS",
              control=list(reltol=1e-6, trace=1, REPORT=1))
timeE = proc.time()[3] - startE
cat(sprintf("\nTime: %.1f sec\n", timeE))
res_E = extractResults(obj_E, opt_E, "Option 4 (random beta, BFGS, truth init)", fixedBeta=FALSE)

# ============================================================
# (F) Options 1+4: fixed beta, BFGS, truth init
# ============================================================
cat("\n============================================================\n")
cat("(F) OPTIONS 1+4: alpha/beta FIXED, BFGS, truth init\n")
cat("============================================================\n")

tmb_params_F = list(
  log_tau = trueLogTau, log_tauEps = trueLogTauEps,
  alpha = trueAlpha, beta = c(trueUrban, trueNormPop),
  u_spatial = rep(0, nAreas),
  nuggetUrbMICS = rep(0, nUrb), nuggetRurMICS = rep(0, nRur)
)

obj_F = MakeADFun(data=data_full, parameters=tmb_params_F,
                  random=c('u_spatial', 'nuggetUrbMICS', 'nuggetRurMICS'),
                  hessian=TRUE, DLL='modM_MIIDonly', silent=TRUE)

startF = proc.time()[3]
opt_F = optim(par=obj_F$par, fn=obj_F$fn, gr=obj_F$gr, method="BFGS",
              control=list(reltol=1e-6, trace=1, REPORT=1))
timeF = proc.time()[3] - startF
cat(sprintf("\nTime: %.1f sec\n", timeF))
res_F = extractResults(obj_F, opt_F, "Options 1+4 (fixed beta, BFGS, truth init)", fixedBeta=TRUE)

# ============================================================
# Summary
# ============================================================
cat("\n\n============================================================\n")
cat("SUMMARY TABLE\n")
cat("============================================================\n\n")

cat(sprintf("%-45s %8s %8s %8s %10s %5s %8s\n",
            "Configuration", "alpha", "urban", "normPop", "NLL", "conv", "time(s)"))
cat(sprintf("%-45s %8s %8s %8s %10s %5s %8s\n",
            "-------------", "-----", "-----", "-------", "---", "----", "-------"))
cat(sprintf("%-45s %8.4f %8.4f %8.4f\n",
            "TRUTH", trueAlpha, trueUrban, trueNormPop))

# Baseline from previous run
cat(sprintf("%-45s %8.4f %8.4f %8.4f %10.2f %5d %8.0f\n",
            "(A) Random beta, BFGS, default init", -0.838, 1.3945, 0.0277, 4355.53, 0, 692))

fmtRow = function(label, res, time) {
  cat(sprintf("%-45s %8.4f %8.4f %8.4f %10.2f %5d %8.0f\n",
              label, res$alpha, res$beta[1], res$beta[2], res$nll, res$conv, time))
}
fmtRow("(C) Fixed beta, BFGS, default init [1]", res_C, timeC)
fmtRow("(E) Random beta, BFGS, truth init [4]", res_E, timeE)
fmtRow("(F) Fixed beta, BFGS, truth init [1+4]", res_F, timeF)

# ============================================================
# Get SEs for all configurations
# ============================================================
getSEs = function(obj, label, fixedBeta) {
  cat(sprintf("\nSEs for %s:\n", label))
  tryCatch({
    SD = TMB::sdreport(obj, getJointPrecision=TRUE,
                       bias.correct=TRUE, bias.correct.control=list(sd=TRUE))
    cat("  pdHess:", SD$pdHess, "\n")
    
    if(fixedBeta) {
      pf = SD$par.fixed
      se = sqrt(diag(SD$cov.fixed))
      pn = names(pf)
      for(i in seq_along(pf)) {
        cat(sprintf("  %s: %.4f (SE %.4f)\n", pn[i], pf[i], se[i]))
      }
    } else {
      pf = SD$par.fixed
      se_f = sqrt(diag(SD$cov.fixed))
      pfn = names(pf)
      for(i in seq_along(pf)) {
        cat(sprintf("  [fixed] %s: %.4f (SE %.4f)\n", pfn[i], pf[i], se_f[i]))
      }
      sr = summary(SD, select="random")
      rn = names(SD$par.random)
      alpha_i = which(rn == "alpha")
      beta_i = which(rn == "beta")
      cat(sprintf("  [random] alpha: %.4f (SE %.4f)\n", sr[alpha_i, 1], sr[alpha_i, 2]))
      for(j in beta_i) {
        cat(sprintf("  [random] beta: %.4f (SE %.4f)\n", sr[j, 1], sr[j, 2]))
      }
    }
    SD
  }, error = function(e) {
    cat("  FAILED:", conditionMessage(e), "\n")
    NULL
  })
}

SD_C = getSEs(obj_C, "(C)", fixedBeta=TRUE)
SD_E = getSEs(obj_E, "(E)", fixedBeta=FALSE)
SD_F = getSEs(obj_F, "(F)", fixedBeta=TRUE)

cat("\nDone.\n")
