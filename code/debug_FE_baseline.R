#!/usr/bin/env Rscript
# Baseline comparison: FE+nugget, beta RANDOM (inner), BFGS with gradient
# Same model/data as debug_optimizer_NM.R but with standard TMB optimization

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
nUrb = length(ysUrbMICS)
nRur = length(ysRurMICS)

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

trueAlpha = -1.25; trueUrban = 1.00; trueNormPop = 0.50
trueSigmaEps = sqrt(1.5)

initUrbP = sum(ysUrbMICS)/sum(nsUrbMICS)
initRurP = sum(ysRurMICS)/sum(nsRurMICS)
initAlpha = log(initRurP/(1-initRurP))
initBeta1 = log(initUrbP/(1-initUrbP)) - initAlpha

unloadDynlibs()
if(!any(file.exists(paste0("code/modM_MIIDonly", c(".o", ".so", ".dll"))))) {
  compile("code/modM_MIIDonly.cpp", framework="TMBad", safebounds=FALSE)
}
dyn.load(dynlib("code/modM_MIIDonly"))
TMB::config(tmbad.sparse_hessian_compress = 1)

# ============================================================
# FE + nugget: beta RANDOM (standard TMB), BFGS
# ============================================================
cat("\n============================================================\n")
cat("BASELINE: FE+NUGGET, beta random, alpha random, BFGS\n")
cat("============================================================\n")

tmb_params = list(
  log_tau = 0, log_tauEps = 0,
  alpha = initAlpha, beta = c(initBeta1, 0),
  u_spatial = rep(0, nAreas),
  nuggetUrbMICS = rep(0, nUrb), nuggetRurMICS = rep(0, nRur)
)

mapList = list(
  u_spatial = factor(rep(NA, nAreas)),
  log_tau   = factor(NA)
)

# Standard: alpha and beta both random (inner Newton)
obj = MakeADFun(data=data_full, parameters=tmb_params,
                random=c('alpha', 'beta', 'nuggetUrbMICS', 'nuggetRurMICS'),
                map=mapList,
                hessian=TRUE, DLL='modM_MIIDonly', silent=TRUE)

cat("Outer params:", names(obj$par), "\n")
cat("Number of outer params:", length(obj$par), "\n")
cat("Number of random effects:", length(obj$env$random), "\n\n")

start = proc.time()[3]
opt = optim(par=obj$par, fn=obj$fn, gr=obj$gr, method="BFGS",
            control=list(reltol=1e-6, trace=1, REPORT=1))
elapsed = proc.time()[3] - start
cat(sprintf("\nTime: %.1f sec\n", elapsed))

cat(sprintf("\n--- BASELINE: FE+nugget, beta random, BFGS ---\n"))
cat("Convergence:", opt$convergence, "\n")
cat("NLL:", opt$value, "\n")

log_tauEps_est = as.numeric(opt$par["log_tauEps"])
lp = obj$env$last.par.best
rn = names(lp)
alpha_est = as.numeric(lp[rn == "alpha"])
beta_est = as.numeric(lp[rn == "beta"])
sigmaEps = exp(-0.5 * log_tauEps_est)

cat(sprintf("  alpha:       %7.4f  (truth: %7.4f, bias: %+.4f)\n",
            alpha_est, trueAlpha, alpha_est - trueAlpha))
cat(sprintf("  beta_urban:  %7.4f  (truth: %7.4f, bias: %+.4f)\n",
            beta_est[1], trueUrban, beta_est[1] - trueUrban))
cat(sprintf("  beta_normPop:%7.4f  (truth: %7.4f, bias: %+.4f)\n",
            beta_est[2], trueNormPop, beta_est[2] - trueNormPop))
cat(sprintf("  sigmaEps:    %7.4f  (truth: %7.4f)\n", sigmaEps, trueSigmaEps))

# SEs
cat("\n--- Computing SEs ---\n")
tryCatch({
  sdr = sdreport(obj)
  summ = summary(sdr, "fixed")
  cat("Fixed parameter SEs:\n")
  print(summ)
  summ_rand = summary(sdr, "random")
  rnames = rownames(summ_rand)
  alpha_row = which(rnames == "alpha")
  beta_rows = which(rnames == "beta")
  cat(sprintf("\n  alpha SE:       %.4f\n", summ_rand[alpha_row, 2]))
  cat(sprintf("  beta_urban SE:  %.4f\n", summ_rand[beta_rows[1], 2]))
  cat(sprintf("  beta_normPop SE:%.4f\n", summ_rand[beta_rows[2], 2]))
}, error = function(e) {
  cat("sdreport failed:", conditionMessage(e), "\n")
})

cat("\nDone.\n")
