#!/usr/bin/env Rscript
# Run only the combined MICS+DHS FE+nug model on simulated survey 1
library(TMB)
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
source("code/makeInputsTMB.R")

simIndex <- 1
KMICS <- 100
KDHSu <- 16
KDHSr <- 21
Qgh <- 10
covs <- c("urban","access","elev","distRiversLakes","normPop")

inputsFile <- sprintf("savedOutput/simStudy1/inputs_sim%d_KMICS%03d_KDHS%02d_%02d.RData", simIndex, KMICS, KDHSu, KDHSr)
if(!file.exists(inputsFile)) stop("Inputs file not found: ", inputsFile)
load(inputsFile) # loads inputsMDM
if(!exists("inputsMDM")) stop("inputsMDM not found in inputs file")

# fix DHS cov names if needed
fix_dhs_names <- function(mat) {
  if(is.null(mat)) return(mat)
  cn <- colnames(mat)
  if(!is.null(cn)) {
    cn <- gsub("^urb$", "urban", cn)
    cn <- gsub("^pop$", "normPop", cn)
    colnames(mat) <- cn
  }
  mat
}
inputsMDM$intPtsDHS$covsUrb <- fix_dhs_names(inputsMDM$intPtsDHS$covsUrb)
inputsMDM$intPtsDHS$covsRur <- fix_dhs_names(inputsMDM$intPtsDHS$covsRur)

# BYM2 args
load("savedOutput/global/admFinalMat.RData")
bym2Args <- prepareBYM2argumentsForTMB(admFinalMat, u=0.5, alpha=2/3, constr=TRUE, scale.model=TRUE, matrixType="TsparseMatrix")

# prepare data for combined model
gh <- fastGHQuad::gaussHermiteData(Qgh)
intM <- inputsMDM$intPtsMICS
intD <- inputsMDM$intPtsDHS

data_comb <- list(
  y_iUrbanMICS = inputsMDM$ysUrbMICS, y_iRuralMICS = inputsMDM$ysRurMICS,
  n_iUrbanMICS = inputsMDM$nsUrbMICS, n_iRuralMICS = inputsMDM$nsRurMICS,
  areaidxlocUrbanMICS = as.integer(inputsMDM$areaidxlocUrbanMICS),
  areaidxlocRuralMICS = as.integer(inputsMDM$areaidxlocRuralMICS),
  X_betaUrbanMICS = as.matrix(intM$XUrb)[, covs],
  X_betaRuralMICS = as.matrix(intM$XRur)[, covs],
  wUrbanMICS = intM$wUrban, wRuralMICS = intM$wRural,

  y_iUrbanDHS = inputsMDM$ysUrbDHS, y_iRuralDHS = inputsMDM$ysRurDHS,
  n_iUrbanDHS = inputsMDM$nsUrbDHS, n_iRuralDHS = inputsMDM$nsRurDHS,
  areaidxlocUrbanDHS = as.integer(inputsMDM$areaidxlocUrbanDHS),
  areaidxlocRuralDHS = as.integer(inputsMDM$areaidxlocRuralDHS),
  X_betaUrbanDHS = as.matrix(intD$covsUrb)[, covs],
  X_betaRuralDHS = as.matrix(intD$covsRur)[, covs],
  wUrbanDHS = intD$wUrban, wRuralDHS = intD$wRural,

  Q_bym2 = bym2Args$Q,
  lchoose_urban_mics = lchoose(inputsMDM$nsUrbMICS, inputsMDM$ysUrbMICS),
  lchoose_rural_mics = lchoose(inputsMDM$nsRurMICS, inputsMDM$ysRurMICS),
  lchoose_urban_dhs = lchoose(inputsMDM$nsUrbDHS, inputsMDM$ysUrbDHS),
  lchoose_rural_dhs = lchoose(inputsMDM$nsRurDHS, inputsMDM$ysRurDHS),
  gh_nodes = gh$x, gh_weights = gh$w,
  alpha_pri = c(0, 100^2), beta_pri = c(0, 10^2),
  tr = bym2Args$tr, gammaTildesm1 = bym2Args$gammaTildesm1,
  lambdaPhi = bym2Args$lambda, lambdaTau = getLambdaPCprec(1, 0.1), lambdaTauEps = getLambdaPCprec(1,0.1),
  options = 0
)

# initial params
initUrbP <- sum(c(inputsMDM$ysUrbMICS, inputsMDM$ysUrbDHS))/sum(c(inputsMDM$nsUrbMICS, inputsMDM$nsUrbDHS))
initRurP <- sum(c(inputsMDM$ysRurMICS, inputsMDM$ysRurDHS))/sum(c(inputsMDM$nsRurMICS, inputsMDM$nsRurDHS))
initAlpha <- logit(initRurP)
initBeta1 <- logit(initUrbP) - initAlpha
nBeta <- length(covs)
nFree <- ncol(bym2Args$Q) - 1
params_comb <- list(log_tau=0, logit_phi=0, log_tauEps=0,
                    alpha=initAlpha,
                    beta=c(initBeta1, rep(0, nBeta-1)),
                    w_bym2Free = rep(0, nFree), u_bym2Free = rep(0, nFree))
map_comb <- list(w_bym2Free = factor(rep(NA, nFree)), u_bym2Free = factor(rep(NA, nFree)), log_tau=factor(NA), logit_phi=factor(NA))

# compile if needed and run
dll <- 'modMDM_BYM2_GH_v2'
if(!any(file.exists(paste0('code/', dll, c('.o','.so','.dll'))))) {
  compile(paste0('code/', dll, '.cpp'), framework='TMBad', safebounds=FALSE)
}
unloadDynlibs(); dyn.load(dynlib(paste0('code/', dll)))

cat('Running combined model...\n')
obj_comb <- MakeADFun(data=data_comb, parameters=params_comb, map=map_comb, DLL=dll, silent=TRUE)

t0 <- proc.time()[3]
opt_comb <- optim(obj_comb$par, obj_comb$fn, obj_comb$gr, method='BFGS', control=list(maxit=2000, reltol=1e-6))
opt_time <- proc.time()[3] - t0

cat('Computing Hessian and SEs...\n')
t0h <- proc.time()[3]
H_comb <- optimHess(opt_comb$par, obj_comb$fn, obj_comb$gr)
se_comb <- sqrt(diag(solve(H_comb)))
se_time <- proc.time()[3] - t0h

# extract
alpha_comb <- opt_comb$par['alpha']; alpha_comb_se <- se_comb['alpha']
beta_comb <- opt_comb$par[grep('^beta', names(opt_comb$par))]; beta_comb_se <- se_comb[grep('^beta', names(se_comb))]

cat('\nCombined model results:\n')
cat(sprintf('alpha = %.4f (%.4f)\n', alpha_comb, alpha_comb_se))
for(i in seq_along(covs)) {
  nm <- covs[i]
  bnm <- beta_comb[i]; bse <- beta_comb_se[i]
  cat(sprintf('beta[%s] = %.4f (%.4f)\n', nm, bnm, bse))
}
cat(sprintf('opt_time=%.2f s; se_time=%.2f s; NLL=%.4f\n', opt_time, se_time, opt_comb$value))

cat('\nDone.\n')
