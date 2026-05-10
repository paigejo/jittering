#!/usr/bin/env Rscript
# DHS-only BYM2 + GH, K=16/21, phi fixed at truth (0.8)

library(TMB)
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")

Qgh = 10
truePhi = 0.8
trueSigmaBYM2 = sqrt(0.5)
trueSigmaEps = sqrt(1.5)
phiLogitFixed = qlogis(truePhi)

# Data
load("savedOutput/simStudy1/simPopsSurveys_BYM2.RData")
ed = surveysDHS[[1]]
nameTab = rbind(c("N", "n"), c("Z", "y"))
for(i in 1:nrow(nameTab)) {
  fromN = nameTab[i,1]; toN = nameTab[i,2]
  if(!(toN %in% names(ed))) ed[[toN]] = ed[[fromN]]
}

isU = ed$urban
ysUrbDHS = ed$y[isU]; ysRurDHS = ed$y[!isU]
nsUrbDHS = ed$n[isU]; nsRurDHS = ed$n[!isU]

# K=16/21 integration points
intFile = "savedOutput/simStudy1/intPtsDHS_simStudy1_1_K16_21.RData"
if(file.exists(intFile)) {
  load(intFile)
} else {
  intPtsDHS = makeAllIntegrationPointsDHS(
    cbind(ed$east, ed$north), ed$urban,
    areaNames = ed$subarea,
    numPointsUrban = 16,
    numPointsRural = 21,
    JInnerUrban = 4, JOuterUrban = 0,
    JInnerRural = 4, JOuterRural = 1,
    popPrior = TRUE,
    adminMap = adm2Full,
    saveOutput = FALSE
  )
  save(intPtsDHS, file=intFile)
}

XUrb = as.matrix(intPtsDHS$covsUrb)
XRur = as.matrix(intPtsDHS$covsRur)
if("int" %in% colnames(XUrb)) XUrb = XUrb[, colnames(XUrb) != "int", drop=FALSE]
if("int" %in% colnames(XRur)) XRur = XRur[, colnames(XRur) != "int", drop=FALSE]

load("savedOutput/global/admFinalMat.RData")
bym2ArgsTMB = prepareBYM2argumentsForTMB(admFinalMat, u=0.5, alpha=2/3,
                                         constr=TRUE, scale.model=TRUE,
                                         matrixType="TsparseMatrix")
AUrbDHS = makeApointToArea(intPtsDHS$areasUrban, admFinal$NAME_FINAL)
ARurDHS = makeApointToArea(intPtsDHS$areasRural, admFinal$NAME_FINAL)
AUrbDHS = t(AUrbDHS); mode(AUrbDHS) = "numeric"
ARurDHS = t(ARurDHS); mode(ARurDHS) = "numeric"

alpha_pri = c(0, 100^2)
beta_pri = c(0, 10^2)
lambdaTau = getLambdaPCprec(u=1, alpha=0.1)
lambdaTauEps = getLambdaPCprec(u=1, alpha=0.1)

initUrbP = sum(ysUrbDHS)/sum(nsUrbDHS)
initRurP = sum(ysRurDHS)/sum(nsRurDHS)
initAlpha = log(initRurP/(1-initRurP))
initBeta1 = log(initUrbP/(1-initUrbP)) - initAlpha

unloadDynlibs()
if(!any(file.exists(paste0("code/modD_BYM2_GH_v2", c(".o", ".so", ".dll"))))) {
  compile("code/modD_BYM2_GH_v2.cpp", framework="TMBad", safebounds=FALSE)
}
dyn.load(dynlib("code/modD_BYM2_GH_v2"))

gh = fastGHQuad::gaussHermiteData(Qgh)
nAreas = ncol(bym2ArgsTMB$Q)
nFree = nAreas - 1
nBeta = ncol(XUrb)

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

# Step 1 FE+nug init with phi fixed and spatial fixed out
params_init = list(
  log_tau=0,
  logit_phi=phiLogitFixed,
  log_tauEps=0,
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

# Step 2 full BYM2+GH with phi fixed
params_full = list(
  log_tau=0,
  logit_phi=phiLogitFixed,
  log_tauEps=as.numeric(opt_init$par["log_tauEps"]),
  alpha=as.numeric(opt_init$par["alpha"]),
  beta=as.numeric(opt_init$par[grep("^beta", names(opt_init$par))]),
  w_bym2Free=rep(0, nFree),
  u_bym2Free=rep(0, nFree)
)
map_full = list(logit_phi=factor(NA))
obj_full = MakeADFun(
  data=data_gh, parameters=params_full,
  map=map_full,
  random=c("alpha", "beta", "w_bym2Free", "u_bym2Free"),
  DLL="modD_BYM2_GH_v2", silent=TRUE
)

start = proc.time()[3]
opt_full = optim(obj_full$par, obj_full$fn, obj_full$gr,
                 method="BFGS", control=list(maxit=500))
time_full = proc.time()[3] - start

sd_full = sdreport(obj_full)
fix = summary(sd_full, "fixed")
rn = rownames(fix)
log_tau = fix[rn=="log_tau", "Estimate"][1]
se_log_tau = fix[rn=="log_tau", "Std. Error"][1]
log_tauEps = fix[rn=="log_tauEps", "Estimate"][1]
se_log_tauEps = fix[rn=="log_tauEps", "Std. Error"][1]

sigmaBYM2 = exp(-0.5 * log_tau)
se_sigmaBYM2 = 0.5 * sigmaBYM2 * se_log_tau
sigmaEps = exp(-0.5 * log_tauEps)
se_sigmaEps = 0.5 * sigmaEps * se_log_tauEps

z_sigmaBYM2 = abs(sigmaBYM2 - trueSigmaBYM2) / se_sigmaBYM2

cat("\n============================================================\n")
cat("DHS BYM2+GH (K=16/21) with phi fixed to truth\n")
cat("============================================================\n")
cat(sprintf("Convergence: %d\n", opt_full$convergence))
cat(sprintf("NLL: %.4f\n", opt_full$value))
cat(sprintf("Time (s): %.1f\n", time_full))
cat(sprintf("phi fixed: %.4f\n", truePhi))
cat(sprintf("sigmaBYM2: %.4f (SE %.4f), truth %.4f, |diff|/SE=%.2f\n",
            sigmaBYM2, se_sigmaBYM2, trueSigmaBYM2, z_sigmaBYM2))
cat(sprintf("sigmaEps: %.4f (SE %.4f), truth %.4f\n",
            sigmaEps, se_sigmaEps, trueSigmaEps))

save(opt_full, sd_full, time_full,
  sigmaBYM2, se_sigmaBYM2,
  sigmaEps, se_sigmaEps,
  z_sigmaBYM2,
  file="savedOutput/fitDHS_BYM2_GH_K16_21_phiFixed.RData")
cat("Saved to savedOutput/fitDHS_BYM2_GH_K16_21_phiFixed.RData\n")
