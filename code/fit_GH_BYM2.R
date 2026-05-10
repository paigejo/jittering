#!/usr/bin/env Rscript
# BYM2 (2n-2 constrained) + GH nuggets, Q=10, BFGS outer
# All covariates (not subsetted), all FEs + spatial as inner params
# Outer: log_tau, logit_phi, log_tauEps

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
if(!("Stratum" %in% names(edMICS)))
  edMICS$Stratum = adm2ToStratumMICS(edMICS$subarea)

inputsMDM = makeInputsMDM(datDHS=ed, datMICS=edMICS,
                          KMICS=100, admMICS=admFinal, adm2DHS=adm2Full)

# Use ALL covariates (no subsetting)
cat("X columns:", colnames(inputsMDM$intPtsMICS$XUrb), "\n")
nBeta = ncol(inputsMDM$intPtsMICS$XUrb)
cat("Number of covariates (beta):", nBeta, "\n")

thisEnv = environment()
list2env(inputsMDM, envir=thisEnv)

# ── BYM2 adjacency / precision setup ────────────────────────────
load("savedOutput/global/admFinalMat.RData")
bym2ArgsTMB = prepareBYM2argumentsForTMB(admFinalMat, u=0.5, alpha=2/3,
                                         constr=TRUE, scale.model=TRUE,
                                         matrixType="TsparseMatrix")
nAreas = ncol(bym2ArgsTMB$Q)
cat("nAreas:", nAreas, "\n")

# ── Priors ───────────────────────────────────────────────────────
alpha_pri = c(0, 100^2)
beta_pri = c(0, sqrt(1000))
lambdaTau = getLambdaPCprec(u=1, alpha=0.1)
lambdaTauEps = getLambdaPCprec(u=1, alpha=0.1)

# ── GH nodes, Q=10 ──────────────────────────────────────────────
Q = 10
gh = fastGHQuad::gaussHermiteData(Q)
cat(sprintf("Using %d Gauss-Hermite quadrature nodes\n", Q))

# ── True values ──────────────────────────────────────────────────
trueAlpha = -1.25
trueGamma = 1.00  # urban
trueBetaRest = c(0, 0, 0, 0.5)  # access, elev, distRiversLakes, normPop
trueSigmaBYM2 = sqrt(0.5)
truePhi = 0.8
trueSigmaEps = sqrt(1.5)
trueLogTau = -2*log(trueSigmaBYM2)
trueLogitPhi = log(truePhi/(1-truePhi))
trueLogTauEps = -2*log(trueSigmaEps)

# ── Initial values ───────────────────────────────────────────────
initUrbP = sum(ysUrbMICS)/sum(nsUrbMICS)
initRurP = sum(ysRurMICS)/sum(nsRurMICS)
initAlpha = log(initRurP/(1-initRurP))
initBeta1 = log(initUrbP/(1-initUrbP)) - initAlpha

# ── TMB Data ─────────────────────────────────────────────────────
data_bym2gh = list(
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
  Q_bym2=bym2ArgsTMB$Q,
  gh_nodes=gh$x,
  gh_weights=gh$w,
  alpha_pri=alpha_pri,
  beta_pri=beta_pri,
  tr=bym2ArgsTMB$tr,
  gammaTildesm1=bym2ArgsTMB$gammaTildesm1,
  lambdaPhi=bym2ArgsTMB$lambda,
  lchoose_urban=lchoose(nsUrbMICS, ysUrbMICS),
  lchoose_rural=lchoose(nsRurMICS, ysRurMICS),
  lambdaTau=lambdaTau,
  lambdaTauEps=lambdaTauEps,
  options=0
)

# Clean up to free memory
rm(surveysMICS, inputsMDM, admFinal, adm2Full); gc()

# ── Compile ──────────────────────────────────────────────────────
unloadDynlibs()
cat("\nCompiling modM_BYM2_GH_v2.cpp...\n")
compile("code/modM_BYM2_GH_v2.cpp", framework="TMBad", safebounds=FALSE)
dyn.load(dynlib("code/modM_BYM2_GH_v2"))
TMB::config(tmbad.sparse_hessian_compress = 1)

# ── Parameters ───────────────────────────────────────────────────
nFree = nAreas - 1
tmb_params = list(
  log_tau = 0,
  logit_phi = 0,
  log_tauEps = 0,
  alpha = initAlpha,
  beta = c(initBeta1, rep(0, nBeta - 1)),
  w_bym2Free = rep(0, nFree),
  u_bym2Free = rep(0, nFree)
)

cat("\n============================================================\n")
cat("BYM2(2n-2) + GH nuggets Q=10: all FEs + spatial inner, BFGS\n")
cat("============================================================\n")

# Inner: alpha, beta, w_bym2Free, u_bym2Free
# Outer: log_tau, logit_phi, log_tauEps
rand_effs = c('alpha', 'beta', 'w_bym2Free', 'u_bym2Free')

cat("Building AD tape...\n"); flush.console()
obj = MakeADFun(data=data_bym2gh, parameters=tmb_params,
                random=rand_effs,
                DLL='modM_BYM2_GH_v2', silent=TRUE)

cat("Outer params:", names(obj$par), "\n")
cat("Random effects:", length(obj$env$random), "params\n")
cat("Tape built. Starting BFGS optimization...\n\n"); flush.console()

t0 = proc.time()[3]
opt = optim(obj$par, obj$fn, obj$gr, method="BFGS",
            control=list(trace=1, REPORT=1, maxit=500))
time_opt = proc.time()[3] - t0

cat(sprintf("\nConvergence: %d\n", opt$convergence))
cat(sprintf("NLL: %.4f\n", opt$value))
cat(sprintf("Time: %.1f s\n", time_opt))

# ── sdreport ─────────────────────────────────────────────────────
cat("Computing sdreport...\n"); flush.console()
sd_rep = sdreport(obj)
est_fix = summary(sd_rep, "fixed")
est_ran = summary(sd_rep, "random")
cat("\nFixed (outer) estimates:\n"); print(est_fix)

pe_f = est_fix[,"Estimate"]; se_f = est_fix[,"Std. Error"]
nms_f = rownames(est_fix)
pe_r = est_ran[,"Estimate"]; se_r = est_ran[,"Std. Error"]
nms_r = rownames(est_ran)

# Extract key estimates
est_alpha = pe_r[nms_r=="alpha"][1]
se_alpha = se_r[nms_r=="alpha"][1]
est_beta = pe_r[nms_r=="beta"]
se_beta = se_r[nms_r=="beta"]
est_logTau = pe_f[nms_f=="log_tau"]
est_logitPhi = pe_f[nms_f=="logit_phi"]
est_logTauEps = pe_f[nms_f=="log_tauEps"]

sigmaBYM2_est = exp(-0.5 * est_logTau)
phi_est = 1/(1 + exp(-est_logitPhi))
sigmaEps_est = exp(-0.5 * est_logTauEps)

covNames = colnames(data_bym2gh$X_betaUrbanMICS)
trueBeta = c(trueGamma, trueBetaRest)

cat(sprintf("\n--- BYM2+GH Results ---\n"))
cat(sprintf("  alpha:      %7.4f (SE %.4f)  truth: %7.4f\n",
            est_alpha, se_alpha, trueAlpha))
for(j in 1:length(est_beta)) {
  cat(sprintf("  beta[%s]: %7.4f (SE %.4f)  truth: %7.4f\n",
              covNames[j], est_beta[j], se_beta[j], trueBeta[j]))
}
cat(sprintf("  sigmaBYM2:  %7.4f           truth: %7.4f\n",
            sigmaBYM2_est, trueSigmaBYM2))
cat(sprintf("  phi:        %7.4f           truth: %7.4f\n",
            phi_est, truePhi))
cat(sprintf("  sigmaEps:   %7.4f           truth: %7.4f\n",
            sigmaEps_est, trueSigmaEps))
cat(sprintf("  Time: %.1f s\n", time_opt))

# ── Save ─────────────────────────────────────────────────────────
result_bym2gh = list(opt=opt, sd=sd_rep, time=time_opt, Q=Q,
                     covNames=covNames, trueBeta=trueBeta)
save(result_bym2gh, file="savedOutput/fitGH_BYM2_Q10.RData")
cat("\nSaved to savedOutput/fitGH_BYM2_Q10.RData\n")
cat("Done.\n")
