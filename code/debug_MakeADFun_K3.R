#!/usr/bin/env Rscript
library(TMB)
setwd('c:/Users/jpaige/git/jittering')
source('code/setup.R')
source('code/modM_DMSep.R')

cat('Building inputs KMICS=3\n')
inputsMDM = makeInputsMDM(ed, edMICS, KMICS=3,
                          KDHSurb=3, JInnerUrban=2,
                          KDHSrur=3, JInnerRural=2, JOuterRural=1,
                          admMICS=admFinal, adm2DHS=adm2Full,
                          adm2AsCovariate=FALSE)
thisEnv = environment(); list2env(inputsMDM, envir=thisEnv)

out = load('savedOutput/global/admFinalMat.RData')
bym2ArgsTMB = prepareBYM2argumentsForTMB(admFinalMat, u=0.5, alpha=2/3, constr=TRUE, scale.model=TRUE, matrixType='TsparseMatrix')
gh = fastGHQuad::gaussHermiteData(5)

# build minimal data_gh using smaller intPts objects (first few rows)
data_gh = list(y_iUrbanMICS=ysUrbMICS[1:10], y_iRuralMICS=ysRurMICS[1:10],
               n_iUrbanMICS=nsUrbMICS[1:10], n_iRuralMICS=nsRurMICS[1:10],
               areaidxlocUrbanMICS=as.integer(areaidxlocUrbanMICS[1:10]), areaidxlocRuralMICS=as.integer(areaidxlocRuralMICS[1:10]),
               X_betaUrbanMICS=matrix(0, nrow=3*10, ncol=ncol(intPtsMICS$XUrb)), X_betaRuralMICS=matrix(0, nrow=3*10, ncol=ncol(intPtsMICS$XRur)),
               wUrbanMICS=matrix(1, nrow=3*10, ncol=1), wRuralMICS=matrix(1, nrow=3*10, ncol=1),
               y_iUrbanDHS=ysUrbDHS[1:10], y_iRuralDHS=ysRurDHS[1:10],
               n_iUrbanDHS=nsUrbDHS[1:10], n_iRuralDHS=nsRurDHS[1:10],
               areaidxlocUrbanDHS=as.integer(areaidxlocUrbanDHS[1:10]), areaidxlocRuralDHS=as.integer(areaidxlocRuralDHS[1:10]),
               X_betaUrbanDHS=matrix(0, nrow=3*10, ncol=ncol(intPtsDHS$covsUrb)), X_betaRuralDHS=matrix(0, nrow=3*10, ncol=ncol(intPtsDHS$covsRur)),
               wUrbanDHS=matrix(1, nrow=3*10, ncol=1), wRuralDHS=matrix(1, nrow=3*10, ncol=1),
               Q_bym2=bym2ArgsTMB$Q,
               lchoose_urban_mics=lchoose(nsUrbMICS[1:10], ysUrbMICS[1:10]), lchoose_rural_mics=lchoose(nsRurMICS[1:10], ysRurMICS[1:10]),
               lchoose_urban_dhs=lchoose(nsUrbDHS[1:10], ysUrbDHS[1:10]), lchoose_rural_dhs=lchoose(nsRurDHS[1:10], ysRurDHS[1:10]),
               gh_nodes=gh$x, gh_weights=gh$w,
               alpha_pri=c(0,100^2), beta_pri=c(0,sqrt(1000)),
               tr=bym2ArgsTMB$tr, gammaTildesm1=bym2ArgsTMB$gammaTildesm1,
               lambdaPhi=bym2ArgsTMB$lambda, lambdaTau=1, lambdaTauEps=1, options=0)

nAreas = ncol(bym2ArgsTMB$Q); nFree = nAreas-1; nBeta = ncol(data_gh$X_betaUrbanMICS)
params_init = list(log_tau=0, logit_phi=0, log_tauEps=0, alpha=0, beta=rep(0,nBeta), w_bym2Free=rep(0,nFree), u_bym2Free=rep(0,nFree))

unloadDynlibs(); dyn.load(dynlib('code/modMDM_BYM2_GH_v2'))
map_fe = list(w_bym2Free=factor(rep(NA, nFree)), u_bym2Free=factor(rep(NA, nFree)), log_tau=factor(NA), logit_phi=factor(NA))

cat('Calling MakeADFun with tiny data...\n')
obj = tryCatch({MakeADFun(data=data_gh, parameters=params_init, map=map_fe, DLL='modMDM_BYM2_GH_v2', silent=TRUE); TRUE}, error=function(e){cat('ERROR:', e$message, '\n'); traceback(); FALSE})
cat('Result:', obj, '\n')
