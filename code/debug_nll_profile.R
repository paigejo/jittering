#!/usr/bin/env Rscript
# Profile the marginal NLL over beta_normPop for the integration-point model
# on BYM2 sim data. This tells us whether the LIKELIHOOD itself
# prefers normPop ≈ 0 (data/model structure issue) or ≈ 0.5 (optimization failure).

library(TMB)
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
source("code/modM_MIIDonly.R")

# Load BYM2 sim data
load("savedOutput/simStudy1/simPopsSurveys_BYM2.RData")
edMICS = surveysMICS[[1]]

# Do name conversions
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

# Subset to urban + normPop
allCovNames = colnames(inputsMDM$intPtsMICS$XUrb)
keepIdx = which(allCovNames %in% c("urban", "normPop"))
XUrbSub = inputsMDM$intPtsMICS$XUrb[, keepIdx, drop=FALSE]
XRurSub = inputsMDM$intPtsMICS$XRur[, keepIdx, drop=FALSE]

cat("XUrb cols:", colnames(XUrbSub), "\n")
cat("XRur cols:", colnames(XRurSub), "\n")

# Build TMB data
lambdaTau = getLambdaPCprec(u=1, alpha=0.1)
lambdaTauEps = getLambdaPCprec(u=1, alpha=0.1)
nAreas = max(c(inputsMDM$areaidxlocUrbanMICS, inputsMDM$areaidxlocRuralMICS)) + 1

data_full = list(
  y_iUrbanMICS=inputsMDM$ysUrbMICS,
  y_iRuralMICS=inputsMDM$ysRurMICS,
  n_iUrbanMICS=inputsMDM$nsUrbMICS,
  n_iRuralMICS=inputsMDM$nsRurMICS,
  areaidxlocUrbanMICS=inputsMDM$areaidxlocUrbanMICS,
  areaidxlocRuralMICS=inputsMDM$areaidxlocRuralMICS,
  X_betaUrbanMICS=XUrbSub,
  X_betaRuralMICS=XRurSub,
  wUrbanMICS=inputsMDM$intPtsMICS$wUrban,
  wRuralMICS=inputsMDM$intPtsMICS$wRural,
  nAreas=nAreas,
  alpha_pri=c(0, 100^2),
  beta_pri=c(0, sqrt(1000)),
  lambdaTau=lambdaTau,
  lambdaTauEps=lambdaTauEps
)

cat("nAreas:", nAreas, "\n")
cat("nUrb:", length(data_full$y_iUrbanMICS), "\n")
cat("nRur:", length(data_full$y_iRuralMICS), "\n")

unloadDynlibs()
if(!any(file.exists(paste0("code/modM_MIIDonly", c(".o", ".so", ".dll"))))) {
  compile("code/modM_MIIDonly.cpp", framework="TMBad", safebounds=FALSE)
}
dyn.load(dynlib("code/modM_MIIDonly"))

# ============================================================
# First: Run full optimization to get baseline estimates
# ============================================================
cat("\n=== Full optimization (unrestricted) ===\n")

initUrbP = sum(data_full$y_iUrbanMICS)/sum(data_full$n_iUrbanMICS)
initRurP = sum(data_full$y_iRuralMICS)/sum(data_full$n_iRuralMICS)
initAlpha = log(initRurP/(1-initRurP))
initBeta1 = log(initUrbP/(1-initUrbP)) - initAlpha

tmb_params = list(
  log_tau = 0,
  log_tauEps = 0,
  alpha = initAlpha,
  beta = c(initBeta1, 0),
  u_spatial = rep(0, nAreas),
  nuggetUrbMICS = rep(0, length(data_full$y_iUrbanMICS)),
  nuggetRurMICS = rep(0, length(data_full$y_iRuralMICS))
)

rand_effs = c('alpha', 'beta', 'u_spatial', 'nuggetUrbMICS', 'nuggetRurMICS')

obj = MakeADFun(data=data_full,
                parameters=tmb_params,
                random=rand_effs,
                hessian=TRUE,
                DLL='modM_MIIDonly',
                silent=TRUE)

opt = optim(par=obj$par, fn=obj$fn, gr=obj$gr, method="BFGS",
            control=list(reltol=1e-6))

cat("Converged:", opt$convergence == 0, "\n")
cat("Outer params (log_tau, log_tauEps):", opt$par, "\n")

# Get inner parameter estimates
nll_baseline = obj$fn(opt$par)
lp = obj$env$last.par.best
rn = names(lp)
est_alpha = lp[rn == "alpha"]
est_beta = lp[rn == "beta"]
est_log_tau = opt$par[1]
est_log_tauEps = opt$par[2]

cat("Estimated alpha:", est_alpha, "\n")
cat("Estimated beta (urban, normPop):", est_beta, "\n")
cat("Baseline NLL:", nll_baseline, "\n")

# ============================================================
# Profile NLL: fix beta_normPop, optimize everything else
# ============================================================
cat("\n=== Profile NLL over beta_normPop ===\n")
cat("(fixing beta[2]=normPop, optimizing log_tau, log_tauEps, and inner params)\n\n")

beta_normPop_grid = seq(-0.5, 1.5, by=0.05)
nll_profile = rep(NA, length(beta_normPop_grid))
alpha_profile = rep(NA, length(beta_normPop_grid))
beta1_profile = rep(NA, length(beta_normPop_grid))

for(gi in seq_along(beta_normPop_grid)) {
  this_normPop = beta_normPop_grid[gi]
  
  tmb_params2 = list(
    log_tau = est_log_tau,
    log_tauEps = est_log_tauEps,
    alpha = est_alpha,
    beta = c(est_beta[1], this_normPop),
    u_spatial = rep(0, nAreas),
    nuggetUrbMICS = rep(0, length(data_full$y_iUrbanMICS)),
    nuggetRurMICS = rep(0, length(data_full$y_iRuralMICS))
  )
  
  # Map: beta[2] fixed (NA), beta[1] free (factor 1)
  mapList2 = list(beta = factor(c(1, NA)))
  
  tryCatch({
    obj2 = MakeADFun(data=data_full,
                     parameters=tmb_params2,
                     random=rand_effs,
                     map=mapList2,
                     hessian=FALSE,
                     DLL='modM_MIIDonly',
                     silent=TRUE)
    
    opt2 = optim(par=obj2$par, fn=obj2$fn, gr=obj2$gr, method="BFGS",
                 control=list(reltol=1e-6, maxit=200))
    
    nll_profile[gi] = opt2$value
    
    # Get inner params
    lp2 = obj2$env$last.par.best
    rn2 = names(lp2)
    alpha_profile[gi] = lp2[rn2 == "alpha"]
    beta1_profile[gi] = lp2[rn2 == "beta"]
    
    cat(sprintf("beta_normPop=%6.3f: NLL=%10.4f  alpha=%7.4f  beta_urban=%7.4f  conv=%d\n",
                this_normPop, nll_profile[gi], alpha_profile[gi], beta1_profile[gi],
                opt2$convergence))
    
  }, error = function(e) {
    cat(sprintf("beta_normPop=%6.3f: FAILED - %s\n", this_normPop, conditionMessage(e)))
  })
}

cat("\n=== Profile Summary ===\n")
best_idx = which.min(nll_profile)
cat("Minimum NLL at beta_normPop =", beta_normPop_grid[best_idx], 
    " (NLL =", nll_profile[best_idx], ")\n")

# Find NLL at truth and at ~0
idx_truth = which.min(abs(beta_normPop_grid - 0.5))
idx_zero = which.min(abs(beta_normPop_grid - 0.0))
cat("NLL at truth (0.5):", nll_profile[idx_truth], "\n")
cat("NLL at ~0:", nll_profile[idx_zero], "\n")
cat("Delta NLL (truth - best):", nll_profile[idx_truth] - nll_profile[best_idx], "\n")
cat("Delta NLL (0 - best):", nll_profile[idx_zero] - nll_profile[best_idx], "\n")

# Alpha at profile minimum vs truth
cat("\nAt profile min:  alpha=", alpha_profile[best_idx], 
    " beta_urban=", beta1_profile[best_idx], "\n")
cat("At truth (0.5):  alpha=", alpha_profile[idx_truth],
    " beta_urban=", beta1_profile[idx_truth], "\n")
cat("True values:     alpha= -1.25  beta_urban= 1.0  beta_normPop= 0.5\n")

# ============================================================
# Also do the same profile for the FE-only model (no spatial RE)
# ============================================================
cat("\n\n=== Profile NLL: FE-only model (no spatial RE) ===\n\n")

# FE-only: map out u_spatial, treat alpha/beta as fixed
tmb_params_fe = list(
  log_tau = 0,
  log_tauEps = 0,
  alpha = initAlpha,
  beta = c(initBeta1, 0),
  u_spatial = rep(0, nAreas),
  nuggetUrbMICS = rep(0, length(data_full$y_iUrbanMICS)),
  nuggetRurMICS = rep(0, length(data_full$y_iRuralMICS))
)

mapList_fe = list(u_spatial = factor(rep(NA, nAreas)),
                  log_tau = factor(NA))
rand_effs_fe = c('nuggetUrbMICS', 'nuggetRurMICS')

obj_fe = MakeADFun(data=data_full,
                   parameters=tmb_params_fe,
                   random=rand_effs_fe,
                   map=mapList_fe,
                   hessian=TRUE,
                   DLL='modM_MIIDonly',
                   silent=TRUE)

opt_fe = optim(par=obj_fe$par, fn=obj_fe$fn, gr=obj_fe$gr, method="BFGS",
               control=list(reltol=1e-6))

cat("FE baseline:\n")
cat("  alpha:", opt_fe$par["alpha"], "\n")
cat("  beta:", opt_fe$par[grep("beta", names(opt_fe$par))], "\n")
cat("  NLL:", opt_fe$value, "\n")

# Profile for FE-only
nll_profile_fe = rep(NA, length(beta_normPop_grid))

for(gi in seq_along(beta_normPop_grid)) {
  this_normPop = beta_normPop_grid[gi]
  
  tmb_params_fe2 = list(
    log_tau = 0,
    log_tauEps = opt_fe$par["log_tauEps"],
    alpha = opt_fe$par["alpha"],
    beta = c(opt_fe$par[grep("beta", names(opt_fe$par))][1], this_normPop),
    u_spatial = rep(0, nAreas),
    nuggetUrbMICS = rep(0, length(data_full$y_iUrbanMICS)),
    nuggetRurMICS = rep(0, length(data_full$y_iRuralMICS))
  )
  
  mapList_fe2 = list(u_spatial = factor(rep(NA, nAreas)),
                     log_tau = factor(NA),
                     beta = factor(c(1, NA)))
  
  tryCatch({
    obj_fe2 = MakeADFun(data=data_full,
                        parameters=tmb_params_fe2,
                        random=rand_effs_fe,
                        map=mapList_fe2,
                        hessian=FALSE,
                        DLL='modM_MIIDonly',
                        silent=TRUE)
    
    opt_fe2 = optim(par=obj_fe2$par, fn=obj_fe2$fn, gr=obj_fe2$gr, method="BFGS",
                    control=list(reltol=1e-6, maxit=200))
    
    nll_profile_fe[gi] = opt_fe2$value
    
  }, error = function(e) {
    cat(sprintf("FE beta_normPop=%6.3f: FAILED\n", this_normPop))
  })
}

best_fe = which.min(nll_profile_fe)
cat("\nFE Profile:\n")
cat("  Min NLL at beta_normPop =", beta_normPop_grid[best_fe], 
    " (NLL =", nll_profile_fe[best_fe], ")\n")
cat("  NLL at truth (0.5):", nll_profile_fe[idx_truth], "\n")
cat("  NLL at ~0:", nll_profile_fe[idx_zero], "\n")
cat("  Delta NLL (truth - best):", nll_profile_fe[idx_truth] - nll_profile_fe[best_fe], "\n")

cat("\nDone.\n")
