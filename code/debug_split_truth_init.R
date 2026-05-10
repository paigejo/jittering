#!/usr/bin/env Rscript
# Config (3) with truth initialization: IID+nug, only normPop outer, NM
# Compare to default-init result: normPop=0.9468, NLL=4338.10

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

# Truth values
trueAlpha = -1.25; trueUrban = 1.00; trueNormPop = 0.50
trueSigmaU = sqrt(0.5); trueSigmaEps = sqrt(1.5)
trueLogTau = -2*log(trueSigmaU)       # 0.6931
trueLogTauEps = -2*log(trueSigmaEps)  # -0.4055

unloadDynlibs()
dyn.load(dynlib("code/modM_MIIDonly_split"))
TMB::config(tmbad.sparse_hessian_compress = 1)

cat("============================================================\n")
cat("IID+NUGGET: only normPop outer, NM — TRUTH INIT\n")
cat("============================================================\n")
cat(sprintf("Starting outer: log_tau=%.4f, log_tauEps=%.4f, beta_normPop=%.4f\n",
            trueLogTau, trueLogTauEps, trueNormPop))
cat(sprintf("Starting random: alpha=%.4f, beta_urban=%.4f\n", trueAlpha, trueUrban))

tmb_params = list(
  log_tau=trueLogTau, log_tauEps=trueLogTauEps,
  alpha=trueAlpha, beta_urban=trueUrban, beta_normPop=trueNormPop,
  u_spatial=rep(0, nAreas),
  nuggetUrbMICS=rep(0, nUrb), nuggetRurMICS=rep(0, nRur)
)

obj = MakeADFun(data=data_full, parameters=tmb_params,
                random=c('alpha','beta_urban','u_spatial','nuggetUrbMICS','nuggetRurMICS'),
                hessian=TRUE, DLL='modM_MIIDonly_split', silent=TRUE)

cat("Outer params:", names(obj$par), "\n")
cat("Starting values:", obj$par, "\n")
cat("Number of outer params:", length(obj$par), "\n")
cat("Number of random effects:", length(obj$env$random), "\n")
cat("Initial NLL:", obj$fn(obj$par), "\n\n")

start = proc.time()[3]
opt = optim(par=obj$par, fn=obj$fn, method="Nelder-Mead",
            control=list(trace=1, maxit=5000, reltol=1e-8))
elapsed = proc.time()[3] - start

pn = names(opt$par)
lp = obj$env$last.par.best; rn = names(lp)

cat(sprintf("\n--- TRUTH INIT result ---\n"))
cat("Convergence:", opt$convergence, "\n")
cat(sprintf("NLL: %.3f\n", opt$value))
cat(sprintf("  alpha:       %7.4f  (truth: %7.4f, bias: %+.4f)\n",
            as.numeric(lp[rn=="alpha"]), trueAlpha,
            as.numeric(lp[rn=="alpha"]) - trueAlpha))
cat(sprintf("  beta_urban:  %7.4f  (truth: %7.4f, bias: %+.4f)\n",
            as.numeric(lp[rn=="beta_urban"]), trueUrban,
            as.numeric(lp[rn=="beta_urban"]) - trueUrban))
cat(sprintf("  beta_normPop:%7.4f  (truth: %7.4f, bias: %+.4f)\n",
            as.numeric(opt$par[pn=="beta_normPop"]), trueNormPop,
            as.numeric(opt$par[pn=="beta_normPop"]) - trueNormPop))
cat(sprintf("  sigmaU:      %7.4f  (truth: %7.4f)\n",
            exp(-0.5*as.numeric(opt$par[pn=="log_tau"])), trueSigmaU))
cat(sprintf("  sigmaEps:    %7.4f  (truth: %7.4f)\n",
            exp(-0.5*as.numeric(opt$par[pn=="log_tauEps"])), trueSigmaEps))
cat(sprintf("  Time: %.1f sec\n", elapsed))

cat("\n--- Comparison ---\n")
cat("                    Default init    Truth init\n")
cat(sprintf("  normPop:          0.9468          %.4f\n", as.numeric(opt$par[pn=="beta_normPop"])))
cat(sprintf("  NLL:              4338.100        %.3f\n", opt$value))
cat(sprintf("  Time:             3670 sec        %.0f sec\n", elapsed))

cat("\nDone.\n")
