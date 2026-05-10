#!/usr/bin/env Rscript
# Compare FE+nugget baseline (beta random) and IID+nugget NM (beta outer)
# Then print summary table across all 4 configurations

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
trueSigmaU = sqrt(0.5); trueSigmaEps = sqrt(1.5)
trueLogTau = -2*log(trueSigmaU); trueLogTauEps = -2*log(trueSigmaEps)

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

results = list()

# ============================================================
# (1) FE+nugget, beta RANDOM, alpha RANDOM, BFGS
# ============================================================
cat("\n============================================================\n")
cat("(1) FE+NUGGET: beta random, alpha random, BFGS\n")
cat("============================================================\n")

tmb_params_1 = list(
  log_tau = 0, log_tauEps = 0,
  alpha = initAlpha, beta = c(initBeta1, 0),
  u_spatial = rep(0, nAreas),
  nuggetUrbMICS = rep(0, nUrb), nuggetRurMICS = rep(0, nRur)
)

mapFE = list(u_spatial = factor(rep(NA, nAreas)), log_tau = factor(NA))

obj_1 = MakeADFun(data=data_full, parameters=tmb_params_1,
                  random=c('alpha', 'beta', 'nuggetUrbMICS', 'nuggetRurMICS'),
                  map=mapFE, hessian=TRUE, DLL='modM_MIIDonly', silent=TRUE)

cat("Outer params:", names(obj_1$par), "\n")
start1 = proc.time()[3]
opt_1 = optim(par=obj_1$par, fn=obj_1$fn, gr=obj_1$gr, method="BFGS",
              control=list(reltol=1e-6, trace=1, REPORT=1))
time1 = proc.time()[3] - start1
cat(sprintf("\nTime: %.1f sec\n", time1))

lp1 = obj_1$env$last.par.best
rn1 = names(lp1)
r1 = list(
  alpha = as.numeric(lp1[rn1 == "alpha"]),
  beta_urban = as.numeric(lp1[rn1 == "beta"])[1],
  beta_normPop = as.numeric(lp1[rn1 == "beta"])[2],
  sigmaU = NA,
  sigmaEps = exp(-0.5 * as.numeric(opt_1$par["log_tauEps"])),
  nll = opt_1$value, conv = opt_1$convergence, time = time1
)
cat(sprintf("  alpha: %.4f, urban: %.4f, normPop: %.4f, sigmaEps: %.4f, NLL: %.2f\n",
            r1$alpha, r1$beta_urban, r1$beta_normPop, r1$sigmaEps, r1$nll))

# SEs
tryCatch({
  sdr1 = sdreport(obj_1)
  sr1 = summary(sdr1, "random")
  rn = rownames(sr1)
  r1$se_alpha = sr1[rn == "alpha", 2]
  r1$se_urban = sr1[rn == "beta", 2][1]
  r1$se_normPop = sr1[rn == "beta", 2][2]
  cat(sprintf("  SE: alpha=%.4f, urban=%.4f, normPop=%.4f\n",
              r1$se_alpha, r1$se_urban, r1$se_normPop))
}, error = function(e) cat("  sdreport failed:", conditionMessage(e), "\n"))

results[["FE+nug, beta rand, BFGS"]] = r1

# ============================================================
# (2) IID+nugget, beta OUTER, alpha RANDOM, Nelder-Mead
# ============================================================
cat("\n============================================================\n")
cat("(2) IID+NUGGET: beta outer, alpha random, Nelder-Mead\n")
cat("============================================================\n")

tmb_params_2 = list(
  log_tau = 0, log_tauEps = 0,
  alpha = initAlpha, beta = c(initBeta1, 0),
  u_spatial = rep(0, nAreas),
  nuggetUrbMICS = rep(0, nUrb), nuggetRurMICS = rep(0, nRur)
)

# No map — spatial effects included
obj_2 = MakeADFun(data=data_full, parameters=tmb_params_2,
                  random=c('alpha', 'u_spatial', 'nuggetUrbMICS', 'nuggetRurMICS'),
                  hessian=TRUE, DLL='modM_MIIDonly', silent=TRUE)

cat("Outer params:", names(obj_2$par), "\n")
cat("Number of outer params:", length(obj_2$par), "\n")
cat("Number of random effects:", length(obj_2$env$random), "\n\n")

start2 = proc.time()[3]
opt_2 = optim(par=obj_2$par, fn=obj_2$fn, method="Nelder-Mead",
              control=list(trace=1, maxit=5000, reltol=1e-8))
time2 = proc.time()[3] - start2
cat(sprintf("\nTime: %.1f sec\n", time2))

cat("Convergence:", opt_2$convergence, "\n")
cat("NLL:", opt_2$value, "\n")

pn2 = names(opt_2$par)
lp2 = obj_2$env$last.par.best
rn2 = names(lp2)
r2 = list(
  alpha = as.numeric(lp2[rn2 == "alpha"]),
  beta_urban = as.numeric(opt_2$par[pn2 == "beta"])[1],
  beta_normPop = as.numeric(opt_2$par[pn2 == "beta"])[2],
  sigmaU = exp(-0.5 * as.numeric(opt_2$par[pn2 == "log_tau"])),
  sigmaEps = exp(-0.5 * as.numeric(opt_2$par[pn2 == "log_tauEps"])),
  nll = opt_2$value, conv = opt_2$convergence, time = time2
)
cat(sprintf("  alpha: %.4f, urban: %.4f, normPop: %.4f, sigmaU: %.4f, sigmaEps: %.4f, NLL: %.2f\n",
            r2$alpha, r2$beta_urban, r2$beta_normPop, r2$sigmaU, r2$sigmaEps, r2$nll))

# SEs
tryCatch({
  sdr2 = sdreport(obj_2)
  sf2 = summary(sdr2, "fixed")
  cat("Fixed parameter SEs:\n")
  print(sf2)
  sr2 = summary(sdr2, "random")
  rn = rownames(sr2)
  r2$se_alpha = sr2[rn == "alpha", 2]
  cat(sprintf("  alpha SE: %.4f\n", r2$se_alpha))
}, error = function(e) cat("  sdreport failed:", conditionMessage(e), "\n"))

results[["IID+nug, beta outer, NM"]] = r2

# ============================================================
# SUMMARY TABLE
# ============================================================
cat("\n\n============================================================\n")
cat("SUMMARY TABLE\n")
cat("============================================================\n\n")

cat(sprintf("%-40s %8s %8s %8s %8s %8s %10s %8s\n",
            "Configuration", "alpha", "urban", "normPop", "sigmaU", "sigmaEps", "NLL", "time(s)"))
cat(sprintf("%-40s %8s %8s %8s %8s %8s %10s %8s\n",
            "----------------------------------------", "------", "------", "-------", "------", "--------", "--------", "-------"))
cat(sprintf("%-40s %8.4f %8.4f %8.4f %8.4f %8.4f\n",
            "TRUTH", trueAlpha, trueUrban, trueNormPop, trueSigmaU, trueSigmaEps))

# (1) FE+nug, beta random
cat(sprintf("%-40s %8.4f %8.4f %8.4f %8s %8.4f %10.2f %8.0f\n",
            "(1) FE+nug, beta rand, BFGS",
            r1$alpha, r1$beta_urban, r1$beta_normPop, "N/A", r1$sigmaEps, r1$nll, r1$time))

# (2) FE+nug, beta outer, NM (from earlier run)
cat(sprintf("%-40s %8.4f %8.4f %8.4f %8s %8.4f %10.2f %8.0f\n",
            "(2) FE+nug, beta outer, NM",
            -1.6828, 0.6316, 0.8909, "N/A", 1.2416, 4430.133, 239))

# (3) IID+nug, beta random, BFGS (from earlier config A/E)
cat(sprintf("%-40s %8.4f %8.4f %8.4f %8.4f %8.4f %10.2f %8.0f\n",
            "(3) IID+nug, beta rand, BFGS",
            -0.838, 1.3945, 0.0277, 0.5270, 1.1987, 4355.53, 692))

# (4) IID+nug, beta outer, NM (just ran)
cat(sprintf("%-40s %8.4f %8.4f %8.4f %8.4f %8.4f %10.2f %8.0f\n",
            "(4) IID+nug, beta outer, NM",
            r2$alpha, r2$beta_urban, r2$beta_normPop, r2$sigmaU, r2$sigmaEps, r2$nll, r2$time))

cat("\nDone.\n")
