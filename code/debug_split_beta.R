#!/usr/bin/env Rscript
# Split-beta test: only beta_normPop in outer, everything else random
# Uses modM_MIIDonly_split.cpp with separate beta_urban / beta_normPop
#
# Configs:
# (1) FE+nug, all beta random, BFGS (baseline)
# (2) IID+nug, all beta random, BFGS (baseline)
# (3) IID+nug, only normPop outer, NM (the fix)

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

# Verify column order
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

# Compile split template
cat("\nCompiling modM_MIIDonly_split.cpp...\n")
compile("code/modM_MIIDonly_split.cpp", framework="TMBad", safebounds=FALSE)
dyn.load(dynlib("code/modM_MIIDonly_split"))
TMB::config(tmbad.sparse_hessian_compress = 1)

# Helper: extract and print results, returning beta as c(urban, normPop) vector
printResults = function(label, opt, alpha_est, beta_urban_est, beta_normPop_est,
                        sigmaU_est=NA, sigmaEps_est, nll, time_sec) {
  beta = c(urban=beta_urban_est, normPop=beta_normPop_est)
  cat(sprintf("\n--- %s ---\n", label))
  cat("Convergence:", opt$convergence, "\n")
  cat(sprintf("NLL: %.3f\n", nll))
  cat(sprintf("  alpha:       %7.4f  (truth: %7.4f, bias: %+.4f)\n",
              alpha_est, trueAlpha, alpha_est - trueAlpha))
  cat(sprintf("  beta_urban:  %7.4f  (truth: %7.4f, bias: %+.4f)\n",
              beta[1], trueUrban, beta[1] - trueUrban))
  cat(sprintf("  beta_normPop:%7.4f  (truth: %7.4f, bias: %+.4f)\n",
              beta[2], trueNormPop, beta[2] - trueNormPop))
  if(!is.na(sigmaU_est))
    cat(sprintf("  sigmaU:      %7.4f  (truth: %7.4f)\n", sigmaU_est, trueSigmaU))
  cat(sprintf("  sigmaEps:    %7.4f  (truth: %7.4f)\n", sigmaEps_est, trueSigmaEps))
  cat(sprintf("  Time: %.1f sec\n", time_sec))
  invisible(list(alpha=alpha_est, beta=beta, sigmaU=sigmaU_est,
                 sigmaEps=sigmaEps_est, nll=nll, conv=opt$convergence, time=time_sec))
}

results = list()

# ============================================================
# (1) FE+nugget, all random (baseline with split template)
# ============================================================
cat("\n============================================================\n")
cat("(1) FE+NUGGET: all random, BFGS (baseline)\n")
cat("============================================================\n")

tmb_params_1 = list(
  log_tau=0, log_tauEps=0,
  alpha=initAlpha, beta_urban=initBeta1, beta_normPop=0,
  u_spatial=rep(0, nAreas),
  nuggetUrbMICS=rep(0, nUrb), nuggetRurMICS=rep(0, nRur)
)
mapFE = list(u_spatial=factor(rep(NA, nAreas)), log_tau=factor(NA))

obj_1 = MakeADFun(data=data_full, parameters=tmb_params_1,
                  random=c('alpha','beta_urban','beta_normPop','nuggetUrbMICS','nuggetRurMICS'),
                  map=mapFE, hessian=TRUE, DLL='modM_MIIDonly_split', silent=TRUE)

cat("Outer params:", names(obj_1$par), "\n")
start1 = proc.time()[3]
opt_1 = optim(par=obj_1$par, fn=obj_1$fn, gr=obj_1$gr, method="BFGS",
              control=list(reltol=1e-6, trace=1, REPORT=1))
time1 = proc.time()[3] - start1

lp1 = obj_1$env$last.par.best; rn1 = names(lp1)
r1 = printResults("(1) FE+nug, all random, BFGS", opt_1,
  alpha_est    = as.numeric(lp1[rn1=="alpha"]),
  beta_urban_est = as.numeric(lp1[rn1=="beta_urban"]),
  beta_normPop_est = as.numeric(lp1[rn1=="beta_normPop"]),
  sigmaEps_est = exp(-0.5*as.numeric(opt_1$par["log_tauEps"])),
  nll=opt_1$value, time_sec=time1)
results[["1"]] = r1

# ============================================================
# (2) IID+nugget, all random, BFGS (baseline)
# ============================================================
cat("\n============================================================\n")
cat("(2) IID+NUGGET: all random, BFGS (baseline)\n")
cat("============================================================\n")

tmb_params_2 = list(
  log_tau=0, log_tauEps=0,
  alpha=initAlpha, beta_urban=initBeta1, beta_normPop=0,
  u_spatial=rep(0, nAreas),
  nuggetUrbMICS=rep(0, nUrb), nuggetRurMICS=rep(0, nRur)
)

obj_2 = MakeADFun(data=data_full, parameters=tmb_params_2,
                  random=c('alpha','beta_urban','beta_normPop','u_spatial','nuggetUrbMICS','nuggetRurMICS'),
                  hessian=TRUE, DLL='modM_MIIDonly_split', silent=TRUE)

cat("Outer params:", names(obj_2$par), "\n")
start2 = proc.time()[3]
opt_2 = optim(par=obj_2$par, fn=obj_2$fn, gr=obj_2$gr, method="BFGS",
              control=list(reltol=1e-6, trace=1, REPORT=1))
time2 = proc.time()[3] - start2

lp2 = obj_2$env$last.par.best; rn2 = names(lp2)
r2 = printResults("(2) IID+nug, all random, BFGS", opt_2,
  alpha_est    = as.numeric(lp2[rn2=="alpha"]),
  beta_urban_est = as.numeric(lp2[rn2=="beta_urban"]),
  beta_normPop_est = as.numeric(lp2[rn2=="beta_normPop"]),
  sigmaU_est   = exp(-0.5*as.numeric(opt_2$par["log_tau"])),
  sigmaEps_est = exp(-0.5*as.numeric(opt_2$par["log_tauEps"])),
  nll=opt_2$value, time_sec=time2)
results[["2"]] = r2

# ============================================================
# (3) IID+nugget, only beta_normPop outer, Nelder-Mead
# ============================================================
cat("\n============================================================\n")
cat("(3) IID+NUGGET: only normPop outer, NM\n")
cat("============================================================\n")

tmb_params_3 = list(
  log_tau=0, log_tauEps=0,
  alpha=initAlpha, beta_urban=initBeta1, beta_normPop=0,
  u_spatial=rep(0, nAreas),
  nuggetUrbMICS=rep(0, nUrb), nuggetRurMICS=rep(0, nRur)
)

# beta_normPop is outer; everything else that varies is random
obj_3 = MakeADFun(data=data_full, parameters=tmb_params_3,
                  random=c('alpha','beta_urban','u_spatial','nuggetUrbMICS','nuggetRurMICS'),
                  hessian=TRUE, DLL='modM_MIIDonly_split', silent=TRUE)

cat("Outer params:", names(obj_3$par), "\n")
cat("Number of outer params:", length(obj_3$par), "\n")
cat("Number of random effects:", length(obj_3$env$random), "\n\n")

start3 = proc.time()[3]
opt_3 = optim(par=obj_3$par, fn=obj_3$fn, method="Nelder-Mead",
              control=list(trace=1, maxit=5000, reltol=1e-8))
time3 = proc.time()[3] - start3

pn3 = names(opt_3$par)
lp3 = obj_3$env$last.par.best; rn3 = names(lp3)
r3 = printResults("(3) IID+nug, normPop outer, NM", opt_3,
  alpha_est    = as.numeric(lp3[rn3=="alpha"]),
  beta_urban_est = as.numeric(lp3[rn3=="beta_urban"]),
  beta_normPop_est = as.numeric(opt_3$par[pn3=="beta_normPop"]),
  sigmaU_est   = exp(-0.5*as.numeric(opt_3$par[pn3=="log_tau"])),
  sigmaEps_est = exp(-0.5*as.numeric(opt_3$par[pn3=="log_tauEps"])),
  nll=opt_3$value, time_sec=time3)
results[["3"]] = r3

# SEs for config 3
cat("\n--- Computing SEs for config (3) ---\n")
tryCatch({
  sdr3 = sdreport(obj_3)
  sf3 = summary(sdr3, "fixed")
  cat("Outer/fixed parameter SEs:\n")
  print(sf3)
  sr3 = summary(sdr3, "random")
  rn = rownames(sr3)
  cat(sprintf("\n  alpha SE:      %.4f\n", sr3[rn=="alpha", 2]))
  cat(sprintf("  beta_urban SE: %.4f\n", sr3[rn=="beta_urban", 2]))
}, error = function(e) cat("  sdreport failed:", conditionMessage(e), "\n"))

# ============================================================
# SUMMARY TABLE
# ============================================================
cat("\n\n============================================================\n")
cat("SUMMARY TABLE\n")
cat("============================================================\n\n")

cat(sprintf("%-42s %8s %8s %8s %8s %8s %10s %8s\n",
            "Config", "alpha", "urban", "normPop", "sigmaU", "sigmaEps", "NLL", "time(s)"))
cat(paste0(rep("-", 106), collapse=""), "\n")
cat(sprintf("%-42s %8.4f %8.4f %8.4f %8.4f %8.4f\n",
            "TRUTH", trueAlpha, trueUrban, trueNormPop, trueSigmaU, trueSigmaEps))

for(k in names(results)) {
  r = results[[k]]
  sigU = if(is.na(r$sigmaU)) "N/A" else sprintf("%8.4f", r$sigmaU)
  cat(sprintf("%-42s %8.4f %8.4f %8.4f %8s %8.4f %10.2f %8.0f\n",
              paste0("(", k, ") ", attr(r, "label")),
              r$alpha, r$beta[1], r$beta[2], sigU, r$sigmaEps, r$nll, r$time))
}

cat("\nDone.\n")
