#!/usr/bin/env Rscript
# DHS-only BYM2 comparison on BYM2-sim data with K=31/41 integration points
# Compares:
#   (1) Laplace (modBYM2JitterDHS.cpp) with nuggets as random effects
#   (2) GH Q=10 (modD_BYM2_GH_v2.cpp) with nuggets marginalized via GH

library(TMB)
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")

Qgh = 10
intFile = "savedOutput/simStudy1/intPtsDHS_simStudy1_1_K31_41.RData"

# -------------------------
# Data
# -------------------------
load("savedOutput/simStudy1/simPopsSurveys_BYM2.RData")
ed = surveysDHS[[1]]
nameTab = rbind(c("N", "n"), c("Z", "y"))
for(i in 1:nrow(nameTab)) {
  fromN = nameTab[i,1]; toN = nameTab[i,2]
  if(!(toN %in% names(ed))) ed[[toN]] = ed[[fromN]]
}

isU = ed$urban
cat("DHS clusters:", nrow(ed), "(urban:", sum(isU), "rural:", sum(!isU), ")\n")

# -------------------------
# Integration points K=31/41
# -------------------------
if(file.exists(intFile)) {
  load(intFile) # loads intPtsDHS
  cat("Loaded", intFile, "\n")
} else {
  cat("Generating K=31/41 integration points...\n")
  cat("Urban: 1 + 3*10 = 31 (JInner=4, JOuter=0)\n")
  cat("Rural: 1 + 3*10 + 1*10 = 41 (JInner=4, JOuter=1)\n")
  intPtsDHS = makeAllIntegrationPointsDHS(
    cbind(ed$east, ed$north), ed$urban,
    areaNames = ed$subarea,
    numPointsUrban = 31,
    numPointsRural = 41,
    JInnerUrban = 4, JOuterUrban = 0,
    JInnerRural = 4, JOuterRural = 1,
    popPrior = TRUE,
    adminMap = adm2Full,
    saveOutput = FALSE
  )
  save(intPtsDHS, file=intFile)
  cat("Saved", intFile, "\n")
}

ysUrbDHS = ed$y[isU]; ysRurDHS = ed$y[!isU]
nsUrbDHS = ed$n[isU]; nsRurDHS = ed$n[!isU]

# Covariates: remove intercept if present
XUrb = as.matrix(intPtsDHS$covsUrb)
XRur = as.matrix(intPtsDHS$covsRur)
if("int" %in% colnames(XUrb)) XUrb = XUrb[, colnames(XUrb) != "int", drop=FALSE]
if("int" %in% colnames(XRur)) XRur = XRur[, colnames(XRur) != "int", drop=FALSE]
nBeta = ncol(XUrb)
covNames = colnames(XUrb)

# BYM2 setup
load("savedOutput/global/admFinalMat.RData")
bym2ArgsTMB = prepareBYM2argumentsForTMB(admFinalMat, u=0.5, alpha=2/3,
                                         constr=TRUE, scale.model=TRUE,
                                         matrixType="TsparseMatrix")
nAreas = ncol(bym2ArgsTMB$Q)
nFree = nAreas - 1

# Area mappings
AUrbDHS = makeApointToArea(intPtsDHS$areasUrban, admFinal$NAME_FINAL)
ARurDHS = makeApointToArea(intPtsDHS$areasRural, admFinal$NAME_FINAL)
AUrbDHS = t(AUrbDHS); mode(AUrbDHS) = "numeric"
ARurDHS = t(ARurDHS); mode(ARurDHS) = "numeric"

# Priors
alpha_pri = c(0, 100^2)
beta_pri = c(0, 10^2)
lambdaTau = getLambdaPCprec(u=1, alpha=0.1)
lambdaTauEps = getLambdaPCprec(u=1, alpha=0.1)

# Initial values
initUrbP = sum(ysUrbDHS)/sum(nsUrbDHS)
initRurP = sum(ysRurDHS)/sum(nsRurDHS)
initAlpha = log(initRurP/(1-initRurP))
initBeta1 = log(initUrbP/(1-initUrbP)) - initAlpha

# Truth
trueAlpha = -1.25
trueBeta = c(1.00, 0, 0, 0, 0.5)
trueSigmaBYM2 = sqrt(0.5)
truePhi = 0.8
trueSigmaEps = sqrt(1.5)

# -------------------------
# (1) Laplace BYM2
# -------------------------
cat("\n============================================================\n")
cat("(1) Laplace BYM2 (K=31/41)\n")
cat("============================================================\n")

unloadDynlibs()
if(!any(file.exists(paste0("code/modBYM2JitterDHS", c(".o", ".so", ".dll"))))) {
  compile("code/modBYM2JitterDHS.cpp", framework="TMBad", safebounds=FALSE)
}
dyn.load(dynlib("code/modBYM2JitterDHS"))

data_lap = list(
  y_iUrbanDHS=ysUrbDHS, y_iRuralDHS=ysRurDHS,
  n_iUrbanDHS=nsUrbDHS, n_iRuralDHS=nsRurDHS,
  AprojUrbanDHS=AUrbDHS, AprojRuralDHS=ARurDHS,
  X_betaUrbanDHS=XUrb, X_betaRuralDHS=XRur,
  wUrbanDHS=intPtsDHS$wUrban, wRuralDHS=intPtsDHS$wRural,
  V_bym2=bym2ArgsTMB$V, Q_bym2=bym2ArgsTMB$Q,
  alpha_pri=alpha_pri, beta_pri=beta_pri,
  tr=bym2ArgsTMB$tr, gammaTildesm1=bym2ArgsTMB$gammaTildesm1,
  lambdaPhi=bym2ArgsTMB$lambda, lambdaTau=lambdaTau, lambdaTauEps=lambdaTauEps,
  options=0
)
params_lap = list(
  alpha=initAlpha,
  beta=c(initBeta1, rep(0, nBeta-1)),
  log_tau=0, logit_phi=0, log_tauEps=0,
  Epsilon_bym2=rep(0, nAreas),
  nuggetUrbDHS=rep(0, length(ysUrbDHS)),
  nuggetRurDHS=rep(0, length(ysRurDHS))
)
obj_lap = MakeADFun(
  data=data_lap, parameters=params_lap,
  random=c("alpha", "beta", "Epsilon_bym2", "nuggetUrbDHS", "nuggetRurDHS"),
  hessian=TRUE,
  DLL="modBYM2JitterDHS", silent=TRUE
)

t0 = proc.time()[3]
opt_lap = optim(obj_lap$par, obj_lap$fn, obj_lap$gr,
                method="BFGS", control=list(maxit=500, reltol=1e-12))
time_lap = proc.time()[3] - t0

sd_lap = sdreport(obj_lap)
fix_lap = summary(sd_lap, "fixed")
ran_lap = summary(sd_lap, "random")

alpha_lap = ran_lap[rownames(ran_lap)=="alpha", "Estimate"][1]
alpha_lap_se = ran_lap[rownames(ran_lap)=="alpha", "Std. Error"][1]
beta_lap = ran_lap[rownames(ran_lap)=="beta", "Estimate"]
beta_lap_se = ran_lap[rownames(ran_lap)=="beta", "Std. Error"]
log_tau_lap = fix_lap[rownames(fix_lap)=="log_tau", "Estimate"][1]
log_tau_lap_se = fix_lap[rownames(fix_lap)=="log_tau", "Std. Error"][1]
logit_phi_lap = fix_lap[rownames(fix_lap)=="logit_phi", "Estimate"][1]
logit_phi_lap_se = fix_lap[rownames(fix_lap)=="logit_phi", "Std. Error"][1]
log_tauEps_lap = fix_lap[rownames(fix_lap)=="log_tauEps", "Estimate"][1]
log_tauEps_lap_se = fix_lap[rownames(fix_lap)=="log_tauEps", "Std. Error"][1]

sigmaBYM2_lap = exp(-0.5 * log_tau_lap)
sigmaBYM2_lap_se = 0.5 * sigmaBYM2_lap * log_tau_lap_se
phi_lap = 1/(1 + exp(-logit_phi_lap))
phi_lap_se = phi_lap * (1 - phi_lap) * logit_phi_lap_se
sigmaEps_lap = exp(-0.5 * log_tauEps_lap)
sigmaEps_lap_se = 0.5 * sigmaEps_lap * log_tauEps_lap_se

cat(sprintf("Laplace: conv=%d, NLL=%.4f, Time=%.1fs\n", opt_lap$convergence, opt_lap$value, time_lap))

# -------------------------
# (2) GH BYM2 (Q=10)
# -------------------------
cat("\n============================================================\n")
cat("(2) GH BYM2 Q=10 (K=31/41)\n")
cat("============================================================\n")

unloadDynlibs()
if(!any(file.exists(paste0("code/modD_BYM2_GH_v2", c(".o", ".so", ".dll"))))) {
  compile("code/modD_BYM2_GH_v2.cpp", framework="TMBad", safebounds=FALSE)
}
dyn.load(dynlib("code/modD_BYM2_GH_v2"))

gh = fastGHQuad::gaussHermiteData(Qgh)

data_gh = list(
  y_iUrbanDHS=ysUrbDHS, y_iRuralDHS=ysRurDHS,
  n_iUrbanDHS=nsUrbDHS, n_iRuralDHS=nsRurDHS,
  AprojUrbanDHS=AUrbDHS, AprojRuralDHS=ARurDHS,
  X_betaUrbanDHS=XUrb, X_betaRuralDHS=XRur,
  wUrbanDHS=intPtsDHS$wUrban, wRuralDHS=intPtsDHS$wRural,
  Q_bym2=bym2ArgsTMB$Q,
  lchoose_urban=lchoose(nsUrbDHS, ysUrbDHS),
  lchoose_rural=lchoose(nsRurDHS, ysRurDHS),
  gh_nodes=gh$x, gh_weights=gh$w,
  alpha_pri=alpha_pri, beta_pri=beta_pri,
  tr=bym2ArgsTMB$tr, gammaTildesm1=bym2ArgsTMB$gammaTildesm1,
  lambdaPhi=bym2ArgsTMB$lambda, lambdaTau=lambdaTau, lambdaTauEps=lambdaTauEps,
  options=0
)

# Step 1: FE+nug GH init (map out spatial + BYM2 hyper)
params_init = list(
  log_tau=0, logit_phi=0, log_tauEps=0,
  alpha=initAlpha,
  beta=c(initBeta1, rep(0, nBeta-1)),
  w_bym2Free=rep(0, nFree),
  u_bym2Free=rep(0, nFree)
)
map_init = list(
  w_bym2Free=factor(rep(NA, nFree)),
  u_bym2Free=factor(rep(NA, nFree)),
  log_tau=factor(NA),
  logit_phi=factor(NA)
)
obj_init = MakeADFun(data=data_gh, parameters=params_init,
                     map=map_init,
                     DLL="modD_BYM2_GH_v2", silent=TRUE)
opt_init = optim(obj_init$par, obj_init$fn, obj_init$gr,
                 method="BFGS", control=list(maxit=500))

# Step 2: full BYM2+GH
params_gh = list(
  log_tau=0,
  logit_phi=0,
  log_tauEps=as.numeric(opt_init$par["log_tauEps"]),
  alpha=as.numeric(opt_init$par["alpha"]),
  beta=as.numeric(opt_init$par[grep("^beta", names(opt_init$par))]),
  w_bym2Free=rep(0, nFree),
  u_bym2Free=rep(0, nFree)
)
obj_gh = MakeADFun(
  data=data_gh, parameters=params_gh,
  random=c("alpha", "beta", "w_bym2Free", "u_bym2Free"),
  DLL="modD_BYM2_GH_v2", silent=TRUE
)

t0 = proc.time()[3]
opt_gh = optim(obj_gh$par, obj_gh$fn, obj_gh$gr,
               method="BFGS", control=list(maxit=500))
time_gh = proc.time()[3] - t0

sd_gh = sdreport(obj_gh)
fix_gh = summary(sd_gh, "fixed")
ran_gh = summary(sd_gh, "random")

alpha_gh = ran_gh[rownames(ran_gh)=="alpha", "Estimate"][1]
alpha_gh_se = ran_gh[rownames(ran_gh)=="alpha", "Std. Error"][1]
beta_gh = ran_gh[rownames(ran_gh)=="beta", "Estimate"]
beta_gh_se = ran_gh[rownames(ran_gh)=="beta", "Std. Error"]
log_tau_gh = fix_gh[rownames(fix_gh)=="log_tau", "Estimate"][1]
log_tau_gh_se = fix_gh[rownames(fix_gh)=="log_tau", "Std. Error"][1]
logit_phi_gh = fix_gh[rownames(fix_gh)=="logit_phi", "Estimate"][1]
logit_phi_gh_se = fix_gh[rownames(fix_gh)=="logit_phi", "Std. Error"][1]
log_tauEps_gh = fix_gh[rownames(fix_gh)=="log_tauEps", "Estimate"][1]
log_tauEps_gh_se = fix_gh[rownames(fix_gh)=="log_tauEps", "Std. Error"][1]

sigmaBYM2_gh = exp(-0.5 * log_tau_gh)
sigmaBYM2_gh_se = 0.5 * sigmaBYM2_gh * log_tau_gh_se
phi_gh = 1/(1 + exp(-logit_phi_gh))
phi_gh_se = phi_gh * (1 - phi_gh) * logit_phi_gh_se
sigmaEps_gh = exp(-0.5 * log_tauEps_gh)
sigmaEps_gh_se = 0.5 * sigmaEps_gh * log_tauEps_gh_se

cat(sprintf("GH: conv=%d, NLL=%.4f, Time=%.1fs\n", opt_gh$convergence, opt_gh$value, time_gh))

# -------------------------
# Summary table
# -------------------------
cat("\n================================================================\n")
cat("DHS-only BYM2 (K=31/41): Laplace vs GH (Q=10)\n")
cat("================================================================\n")
cat(sprintf("%-25s %10s %8s %10s %8s %10s\n", "Parameter", "Laplace", "(SE)", "GH", "(SE)", "Truth"))
cat(sprintf("%-25s %10.4f %8.4f %10.4f %8.4f %10.4f\n", "alpha", alpha_lap, alpha_lap_se, alpha_gh, alpha_gh_se, trueAlpha))
for(j in 1:nBeta) {
  truthj = if(j <= length(trueBeta)) trueBeta[j] else NA
  cat(sprintf("%-25s %10.4f %8.4f %10.4f %8.4f %10.4f\n",
              paste0("beta[", covNames[j], "]"),
              beta_lap[j], beta_lap_se[j], beta_gh[j], beta_gh_se[j], truthj))
}
cat(sprintf("%-25s %10.4f %8.4f %10.4f %8.4f %10.4f\n",
            "sigmaBYM2", sigmaBYM2_lap, sigmaBYM2_lap_se, sigmaBYM2_gh, sigmaBYM2_gh_se, trueSigmaBYM2))
cat(sprintf("%-25s %10.4f %8.4f %10.4f %8.4f %10.4f\n",
            "phi", phi_lap, phi_lap_se, phi_gh, phi_gh_se, truePhi))
cat(sprintf("%-25s %10.4f %8.4f %10.4f %8.4f %10.4f\n",
            "sigmaEps", sigmaEps_lap, sigmaEps_lap_se, sigmaEps_gh, sigmaEps_gh_se, trueSigmaEps))
cat(sprintf("%-25s %10.4f %8s %10.4f\n", "NLL", opt_lap$value, "", opt_gh$value))
cat(sprintf("%-25s %10.1f %8s %10.1f\n", "Time (s)", time_lap, "", time_gh))

cat("\n--- Distance from truth (in SEs) ---\n")
cat(sprintf("%-25s %10s %10s\n", "Parameter", "Laplace", "GH"))
cat(sprintf("%-25s %10.2f %10.2f\n", "alpha",
            abs(alpha_lap - trueAlpha)/alpha_lap_se,
            abs(alpha_gh - trueAlpha)/alpha_gh_se))
for(j in 1:nBeta) {
  if(j <= length(trueBeta) && !is.na(trueBeta[j]) && trueBeta[j] != 0) {
    cat(sprintf("%-25s %10.2f %10.2f\n", paste0("beta[", covNames[j], "]"),
                abs(beta_lap[j] - trueBeta[j])/beta_lap_se[j],
                abs(beta_gh[j] - trueBeta[j])/beta_gh_se[j]))
  }
}
cat(sprintf("%-25s %10.2f %10.2f\n", "sigmaBYM2",
            abs(sigmaBYM2_lap - trueSigmaBYM2)/sigmaBYM2_lap_se,
            abs(sigmaBYM2_gh - trueSigmaBYM2)/sigmaBYM2_gh_se))
cat(sprintf("%-25s %10.2f %10.2f\n", "phi",
            abs(phi_lap - truePhi)/phi_lap_se,
            abs(phi_gh - truePhi)/phi_gh_se))
cat(sprintf("%-25s %10.2f %10.2f\n", "sigmaEps",
            abs(sigmaEps_lap - trueSigmaEps)/sigmaEps_lap_se,
            abs(sigmaEps_gh - trueSigmaEps)/sigmaEps_gh_se))

cat("\nDone.\n")
