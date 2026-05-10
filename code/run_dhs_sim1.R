#!/usr/bin/env Rscript
# Run DHS-only FE+nug model on simulated survey 1
library(TMB)
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")

simIndex <- 1
KMICS <- 100
KDHSu <- 16
KDHSr <- 21
Qgh <- 10
covs <- c("urban","access","elev","distRiversLakes","normPop")

# load simulated surveys and inputs
load('savedOutput/simStudy1/simPopsSurveys_BYM2.RData')
datDHS <- surveysDHS[[simIndex]]
inputsFile <- sprintf("savedOutput/simStudy1/inputs_sim%d_KMICS%03d_KDHS%02d_%02d.RData", simIndex, KMICS, KDHSu, KDHSr)
if(!file.exists(inputsFile)) stop('Inputs file not found: ', inputsFile)
load(inputsFile) # loads inputsMDM

# bym2 args
load("savedOutput/global/admFinalMat.RData")
bym2Args <- prepareBYM2argumentsForTMB(admFinalMat, u=0.5, alpha=2/3, constr=TRUE, scale.model=TRUE, matrixType="TsparseMatrix")

# prepare DHS-only data
intD <- inputsMDM$intPtsDHS
AUrbDHS <- inputsMDM$AUrbDHS
ARurDHS <- inputsMDM$ARurDHS

ysUrbDHS <- inputsMDM$ysUrbDHS
ysRurDHS <- inputsMDM$ysRurDHS
nsUrbDHS <- inputsMDM$nsUrbDHS
nsRurDHS <- inputsMDM$nsRurDHS

gh <- fastGHQuad::gaussHermiteData(Qgh)
## normalize DHS cov names to canonical names
fix_dhs_names <- function(mat) {
  if(is.null(mat)) return(mat)
  cn <- colnames(mat)
  if(!is.null(cn)) {
    cn <- gsub('^urb$', 'urban', cn)
    cn <- gsub('^pop$', 'normPop', cn)
    colnames(mat) <- cn
  }
  mat
}
intD$covsUrb <- fix_dhs_names(intD$covsUrb)
intD$covsRur <- fix_dhs_names(intD$covsRur)

data_dhs <- list(
  y_iUrbanDHS = ysUrbDHS, y_iRuralDHS = ysRurDHS,
  n_iUrbanDHS = nsUrbDHS, n_iRuralDHS = nsRurDHS,
  AprojUrbanDHS = AUrbDHS, AprojRuralDHS = ARurDHS,
  X_betaUrbanDHS = as.matrix(intD$covsUrb)[, covs],
  X_betaRuralDHS = as.matrix(intD$covsRur)[, covs],
  wUrbanDHS = intD$wUrban, wRuralDHS = intD$wRural,
  Q_bym2 = bym2Args$Q,
  lchoose_urban = lchoose(nsUrbDHS, ysUrbDHS),
  lchoose_rural = lchoose(nsRurDHS, ysRurDHS),
  gh_nodes = gh$x, gh_weights = gh$w,
  alpha_pri = c(0, 100^2), beta_pri = c(0, 10^2),
  tr = bym2Args$tr, gammaTildesm1 = bym2Args$gammaTildesm1,
  lambdaPhi = bym2Args$lambda, lambdaTau = getLambdaPCprec(1,0.1), lambdaTauEps = getLambdaPCprec(1,0.1),
  options=0
)

nFree <- ncol(bym2Args$Q) - 1
nBeta <- ncol(data_dhs$X_betaUrbanDHS)

initUrbP_dhs <- sum(ysUrbDHS)/sum(nsUrbDHS)
initRurP_dhs <- sum(ysRurDHS)/sum(nsRurDHS)
initAlpha_dhs <- log(initRurP_dhs/(1-initRurP_dhs))
initBeta1_dhs <- log(initUrbP_dhs/(1-initUrbP_dhs)) - initAlpha_dhs
params_dhs <- list(log_tau=0, logit_phi=0, log_tauEps=0,
                   alpha=initAlpha_dhs,
                   beta=c(initBeta1_dhs, rep(0, nBeta-1)),
                   w_bym2Free = rep(0, nFree), u_bym2Free = rep(0, nFree))
map_dhs <- list(w_bym2Free = factor(rep(NA, nFree)), u_bym2Free = factor(rep(NA, nFree)), log_tau=factor(NA), logit_phi=factor(NA))

# compile/load dll
dllName <- 'modD_BYM2_GH_v2'
if(!any(file.exists(paste0('code/', dllName, c('.o', '.so', '.dll'))))) {
  compile(paste0('code/', dllName, '.cpp'), framework='TMBad', safebounds=FALSE)
}
unloadDynlibs(); dyn.load(dynlib(paste0('code/', dllName)))

cat('Running DHS-only model...\n')
obj_dhs <- MakeADFun(data=data_dhs, parameters=params_dhs, map=map_dhs, DLL=dllName, silent=TRUE)

t0 <- proc.time()[3]
opt_dhs <- optim(obj_dhs$par, obj_dhs$fn, obj_dhs$gr, method='BFGS', control=list(maxit=2000, reltol=1e-6))
time_dhs <- proc.time()[3] - t0

cat('Computing Hessian and SEs...\n')
t0h <- proc.time()[3]
H_dhs <- optimHess(opt_dhs$par, obj_dhs$fn, obj_dhs$gr)
se_dhs <- sqrt(diag(solve(H_dhs)))
se_time <- proc.time()[3] - t0h

# extract
dhs_alpha <- as.numeric(opt_dhs$par['alpha']); dhs_alpha_se <- as.numeric(se_dhs['alpha'])
dhs_beta <- as.numeric(opt_dhs$par[grep('^beta', names(opt_dhs$par))]); dhs_beta_se <- as.numeric(se_dhs[grep('^beta', names(se_dhs))])

cat('\nDHS-only results:\n')
cat(sprintf('alpha = %.4f (%.4f)\n', dhs_alpha, dhs_alpha_se))
for(j in seq_along(dhs_beta)) {
  cat(sprintf('beta[%s] = %.4f (%.4f)\n', covs[j], dhs_beta[j], dhs_beta_se[j]))
}
cat(sprintf('opt_time=%.2f s; se_time=%.2f s; NLL=%.4f\n', time_dhs, se_time, opt_dhs$value))

cat('\nDone.\n')
