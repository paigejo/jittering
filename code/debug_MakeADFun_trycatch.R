#!/usr/bin/env Rscript
library(TMB)
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
source("code/modM_DMSep.R")

cat('Building inputs via makeInputsMDM\n')
inputsMDM = makeInputsMDM(ed, edMICS,
                          KMICS=25,
                          KDHSurb=16, JInnerUrban=4,
                          KDHSrur=21, JInnerRural=4, JOuterRural=1,
                          admMICS=admFinal, adm2DHS=adm2Full,
                          adm2AsCovariate=FALSE)
thisEnv = environment(); list2env(inputsMDM, envir=thisEnv)

# print key object sizes
cat('nUrbMICS=', if(exists('nsUrbMICS')) length(nsUrbMICS) else 'MISSING', '\n')
cat('nRurMICS=', if(exists('nsRurMICS')) length(nsRurMICS) else 'MISSING', '\n')
cat('nUrbDHS=', if(exists('nsUrbDHS')) length(nsUrbDHS) else 'MISSING', '\n')
cat('nRurDHS=', if(exists('nsRurDHS')) length(nsRurDHS) else 'MISSING', '\n')
cat('areaidxlocUrbanMICS len=', if(exists('areaidxlocUrbanMICS')) length(areaidxlocUrbanMICS) else 'MISSING', '\n')
cat('areaidxlocRuralMICS len=', if(exists('areaidxlocRuralMICS')) length(areaidxlocRuralMICS) else 'MISSING', '\n')
cat('areaidxlocUrbanDHS len=', if(exists('areaidxlocUrbanDHS')) length(areaidxlocUrbanDHS) else 'MISSING', '\n')
cat('areaidxlocRuralDHS len=', if(exists('areaidxlocRuralDHS')) length(areaidxlocRuralDHS) else 'MISSING', '\n')
cat('X_betaUrbanMICS dims=', if(exists('intPtsMICS')) paste(dim(intPtsMICS$XUrb), collapse='x') else 'MISSING', '\n')
cat('X_betaRuralMICS dims=', if(exists('intPtsMICS')) paste(dim(intPtsMICS$XRur), collapse='x') else 'MISSING', '\n')

out = load('savedOutput/global/admFinalMat.RData')
bym2ArgsTMB = prepareBYM2argumentsForTMB(admFinalMat, u=0.5, alpha=2/3,
                                         constr=TRUE, scale.model=TRUE, matrixType='TsparseMatrix')

gh = fastGHQuad::gaussHermiteData(10)

data_gh = list(
  y_iUrbanMICS=ysUrbMICS, y_iRuralMICS=ysRurMICS,
  n_iUrbanMICS=nsUrbMICS, n_iRuralMICS=nsRurMICS,
  areaidxlocUrbanMICS=as.integer(areaidxlocUrbanMICS), areaidxlocRuralMICS=as.integer(areaidxlocRuralMICS),
  X_betaUrbanMICS=intPtsMICS$XUrb, X_betaRuralMICS=intPtsMICS$XRur,
  wUrbanMICS=intPtsMICS$wUrban, wRuralMICS=intPtsMICS$wRural,
  y_iUrbanDHS=ysUrbDHS, y_iRuralDHS=ysRurDHS,
  n_iUrbanDHS=nsUrbDHS, n_iRuralDHS=nsRurDHS,
  areaidxlocUrbanDHS=as.integer(areaidxlocUrbanDHS), areaidxlocRuralDHS=as.integer(areaidxlocRuralDHS),
  X_betaUrbanDHS=intPtsDHS$covsUrb, X_betaRuralDHS=intPtsDHS$covsRur,
  wUrbanDHS=intPtsDHS$wUrban, wRuralDHS=intPtsDHS$wRural,
  Q_bym2=bym2ArgsTMB$Q,
  lchoose_urban_mics=lchoose(nsUrbMICS, ysUrbMICS), lchoose_rural_mics=lchoose(nsRurMICS, ysRurMICS),
  lchoose_urban_dhs=lchoose(nsUrbDHS, ysUrbDHS), lchoose_rural_dhs=lchoose(nsRurDHS, ysRurDHS),
  gh_nodes=gh$x, gh_weights=gh$w,
  alpha_pri=c(0, 100^2), beta_pri=c(0, sqrt(1000)),
  tr=bym2ArgsTMB$tr, gammaTildesm1=bym2ArgsTMB$gammaTildesm1,
  lambdaPhi=bym2ArgsTMB$lambda, lambdaTau=1, lambdaTauEps=1,
  options=0
)

nAreas = ncol(bym2ArgsTMB$Q); nFree = nAreas - 1; nBeta = ncol(data_gh$X_betaUrbanMICS)
initUrbP = sum(c(ysUrbMICS, ysUrbDHS))/sum(c(nsUrbMICS, nsUrbDHS))
initRurP = sum(c(ysRurMICS, ysRurDHS))/sum(c(nsRurMICS, nsRurDHS))
initAlpha = logit(initRurP); initBeta1 = logit(initUrbP) - initAlpha
params_init = list(log_tau=0, logit_phi=0, log_tauEps=0, alpha=initAlpha,
                   beta=c(initBeta1, rep(0, nBeta-1)), w_bym2Free=rep(0, nFree), u_bym2Free=rep(0, nFree))

unloadDynlibs(); dyn.load(dynlib('code/modMDM_BYM2_GH_v2'))
map_fe = list(w_bym2Free=factor(rep(NA, nFree)), u_bym2Free=factor(rep(NA, nFree)), log_tau=factor(NA), logit_phi=factor(NA))

cat('About to call MakeADFun...\n')
res = tryCatch({
  obj_fe = MakeADFun(data=data_gh, parameters=params_init, map=map_fe, DLL='modMDM_BYM2_GH_v2', silent=TRUE)
  cat('MakeADFun finished OK\n')
  TRUE
}, error=function(e){
  cat('MAKEADFUN ERROR:\n')
  print(e)
  traceback()
  FALSE
}, warning=function(w){
  cat('MAKEADFUN WARNING:\n')
  print(w)
  invokeRestart('muffleWarning')
})

if(!isTRUE(res)){
  cat('MakeADFun failed; saving workspace to debug_MakeADFun_fail.RData\n')
  save.image('debug_MakeADFun_fail.RData')
}

cat('Script complete.\n')
