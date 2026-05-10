#!/usr/bin/env Rscript
library(TMB)
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
source("code/modM_DMSep.R")

cat("step 1: makeInputsMDM\n")
inputsMDM = makeInputsMDM(ed, edMICS,
                          KMICS=100,
                          KDHSurb=16, JInnerUrban=4,
                          KDHSrur=21, JInnerRural=4, JOuterRural=1,
                          admMICS=admFinal, adm2DHS=adm2Full,
                          adm2AsCovariate=FALSE)
cat("step 1 done\n")

thisEnv = environment()
list2env(inputsMDM, envir=thisEnv)

out = load("savedOutput/global/admFinalMat.RData")
bym2ArgsTMB = prepareBYM2argumentsForTMB(admFinalMat, u=0.5, alpha=2/3,
                                         constr=TRUE, scale.model=TRUE, matrixType="TsparseMatrix")
cat("step 2: bym2 done\n")

lambdaTau = getLambdaPCprec(u=1, alpha=0.5)
lambdaTauEps = getLambdaPCprec(u=1, alpha=0.5)

gh = fastGHQuad::gaussHermiteData(10)

cat("step 3: build data list\n")
data_gh = list(
  y_iUrbanMICS=ysUrbMICS, y_iRuralMICS=ysRurMICS,
  n_iUrbanMICS=nsUrbMICS, n_iRuralMICS=nsRurMICS,
  areaidxlocUrbanMICS=as.integer(areaidxlocUrbanMICS),
  areaidxlocRuralMICS=as.integer(areaidxlocRuralMICS),
  X_betaUrbanMICS=intPtsMICS$XUrb,
  X_betaRuralMICS=intPtsMICS$XRur,
  wUrbanMICS=intPtsMICS$wUrban,
  wRuralMICS=intPtsMICS$wRural,
  y_iUrbanDHS=ysUrbDHS, y_iRuralDHS=ysRurDHS,
  n_iUrbanDHS=nsUrbDHS, n_iRuralDHS=nsRurDHS,
  areaidxlocUrbanDHS=as.integer(areaidxlocUrbanDHS),
  areaidxlocRuralDHS=as.integer(areaidxlocRuralDHS),
  X_betaUrbanDHS=intPtsDHS$covsUrb,
  X_betaRuralDHS=intPtsDHS$covsRur,
  wUrbanDHS=intPtsDHS$wUrban,
  wRuralDHS=intPtsDHS$wRural,
  Q_bym2=bym2ArgsTMB$Q,
  lchoose_urban_mics=lchoose(nsUrbMICS, ysUrbMICS),
  lchoose_rural_mics=lchoose(nsRurMICS, ysRurMICS),
  lchoose_urban_dhs=lchoose(nsUrbDHS, ysUrbDHS),
  lchoose_rural_dhs=lchoose(nsRurDHS, ysRurDHS),
  gh_nodes=gh$x, gh_weights=gh$w,
  alpha_pri=c(0, 100^2), beta_pri=c(0, sqrt(1000)),
  tr=bym2ArgsTMB$tr, gammaTildesm1=bym2ArgsTMB$gammaTildesm1,
  lambdaPhi=bym2ArgsTMB$lambda, lambdaTau=lambdaTau, lambdaTauEps=lambdaTauEps,
  options=0
)
cat("step 3 done\n")

nAreas = ncol(bym2ArgsTMB$Q)
nFree = nAreas - 1
nBeta = ncol(data_gh$X_betaUrbanMICS)

initUrbP = sum(c(ysUrbMICS, ysUrbDHS))/sum(c(nsUrbMICS, nsUrbDHS))
initRurP = sum(c(ysRurMICS, ysRurDHS))/sum(c(nsRurMICS, nsRurDHS))
initAlpha = logit(initRurP)
initBeta1 = logit(initUrbP) - initAlpha

params_init = list(
  log_tau=0, logit_phi=0, log_tauEps=0,
  alpha=initAlpha,
  beta=c(initBeta1, rep(0, nBeta-1)),
  w_bym2Free=rep(0, nFree),
  u_bym2Free=rep(0, nFree)
)

cat("step 4: unload/load dll\n")
unloadDynlibs()
dyn.load(dynlib("code/modMDM_BYM2_GH_v2"))
cat("step 4 done\n")

map_fe = list(
  w_bym2Free=factor(rep(NA, nFree)),
  u_bym2Free=factor(rep(NA, nFree)),
  log_tau=factor(NA),
  logit_phi=factor(NA)
)

cat("step 5: MakeADFun\n")
obj_fe = MakeADFun(data=data_gh, parameters=params_init, map=map_fe,
                   DLL="modMDM_BYM2_GH_v2", silent=FALSE)
cat("step 5 done\n")

cat("step 6: one fn eval\n")
val = obj_fe$fn(obj_fe$par)
cat("fn val:", val, "\n")
cat("done\n")
