#!/usr/bin/env Rscript
# Run DHS-only BYM2 (GH) and IID+Nugget (TMB) and compare results
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
source("code/makeInputsTMB.R")
source("code/modFED.R")

Qgh <- 10
canonNames <- c("urban", "access", "elev", "distRiversLakes", "normPop")

# Run BYM2 GH for DHS K=16/21 via fitFED
cat("\n=== Running DHS FE+nug GH via fitFED (K=16/21) ===\n")
res_bym2_fit <- fitFED(datDHS=ed, KDHSurb=16, KDHSrur=21, Qgh=Qgh,
                       fixedEffectsOnly=TRUE, getSDs=TRUE, verbose=FALSE)

pe_bym2 <- summary(res_bym2_fit$TMBsd, "fixed")
res_bym2 <- list(
  alpha = pe_bym2["alpha","Estimate"], alpha_se = pe_bym2["alpha","Std. Error"],
  beta = pe_bym2[grep("^beta", rownames(pe_bym2)), "Estimate"],
  beta_se = pe_bym2[grep("^beta", rownames(pe_bym2)), "Std. Error"],
  covNames = res_bym2_fit$covNames,
  sigEps = exp(-0.5 * pe_bym2["log_tauEps","Estimate"]),
  NLL = res_bym2_fit$opt$objective,
  time = res_bym2_fit$totalTime
)
cat("BYM2 (DHS K=16/21) run complete.\n")

# ---- IID + Nugget model on DHS-only (using modM_MIIDonly.cpp) ----
cat("\n=== Running IID + nugget TMB model on DHS-only ===\n")
inputsMDM <- makeInputsMDM(datDHS=ed, datMICS=edMICS, KMICS=100, KDHSurb=16, KDHSrur=21)
intD <- inputsMDM$intPtsDHS

XUrb <- as.matrix(intD$covsUrb)
XRur <- as.matrix(intD$covsRur)
if(!is.null(colnames(XUrb))) {
  cn <- colnames(XUrb)
  cn <- gsub("^urb$", "urban", cn)
  cn <- gsub("^pop$", "normPop", cn)
  colnames(XUrb) <- cn
}
if(!is.null(colnames(XRur))) {
  cn2 <- colnames(XRur)
  cn2 <- gsub("^urb$", "urban", cn2)
  cn2 <- gsub("^pop$", "normPop", cn2)
  colnames(XRur) <- cn2
}
if("int" %in% colnames(XUrb)) XUrb <- XUrb[, colnames(XUrb) != "int", drop=FALSE]
if("int" %in% colnames(XRur)) XRur <- XRur[, colnames(XRur) != "int", drop=FALSE]

# canonical covariate order
covNamesIID <- colnames(XUrb)
cat("IID covariates:", paste(covNamesIID, collapse=", "), "\n")

ysU <- inputsMDM$ysUrbDHS; ysR <- inputsMDM$ysRurDHS
nsU <- inputsMDM$nsUrbDHS; nsR <- inputsMDM$nsRurDHS
areaU <- inputsMDM$areaidxlocUrbanDHS; areaR <- inputsMDM$areaidxlocRuralDHS

nAreas <- max(c(areaU, areaR)) + 1

# priors
lambdaTau <- getLambdaPCprec(1, 0.1)
lambdaTauEps <- getLambdaPCprec(1, 0.1)

# initial values
initUrbP <- sum(ysU)/sum(nsU)
initRurP <- sum(ysR)/sum(nsR)
initAlpha <- log(initRurP/(1-initRurP))
initBeta1 <- log(initUrbP/(1-initUrbP)) - initAlpha

nBeta <- ncol(XUrb)
params_iid <- list(
  log_tau=0, log_tauEps=0,
  alpha=initAlpha,
  beta=c(initBeta1, rep(0, nBeta-1)),
  u_spatial=rep(0, nAreas),
  nuggetUrbMICS=rep(0, length(ysU)),
  nuggetRurMICS=rep(0, length(ysR))
)

data_iid <- list(
  y_iUrbanMICS=ysU, y_iRuralMICS=ysR,
  n_iUrbanMICS=nsU, n_iRuralMICS=nsR,
  areaidxlocUrbanMICS=areaU, areaidxlocRuralMICS=areaR,
  X_betaUrbanMICS=XUrb, X_betaRuralMICS=XRur,
  wUrbanMICS=intD$wUrban, wRuralMICS=intD$wRural,
  nAreas=nAreas,
  alpha_pri=c(0,100^2), beta_pri=c(0,10^2),
  lambdaTau=lambdaTau, lambdaTauEps=lambdaTauEps
)

# compile/load IID TMB model
unloadDynlibs()
dllIID <- "modM_MIIDonly"
if(!any(file.exists(paste0("code/", dllIID, c(".o", ".so", ".dll"))))) {
  cat("Compiling code/", dllIID, "...\n")
  compile(paste0("code/", dllIID, ".cpp"), framework="TMBad", safebounds=FALSE)
}

dyn.load(dynlib(paste0("code/", dllIID)))
TMB::config(tmbad.sparse_hessian_compress = 1)

rand_effs <- c('alpha', 'beta', 'u_spatial', 'nuggetUrbMICS', 'nuggetRurMICS')
obj_iid <- MakeADFun(data=data_iid, parameters=params_iid, random=rand_effs, hessian=TRUE, DLL=dllIID)

cat("Optimizing IID model (this may take a bit)...\n")
t0 <- proc.time()[3]
opt_iid <- optim(obj_iid$par, obj_iid$fn, obj_iid$gr, method="BFGS", control=list(maxit=1000, reltol=1e-6))
time_iid <- proc.time()[3] - t0

H_iid <- optimHess(opt_iid$par, obj_iid$fn, obj_iid$gr)
se_iid <- sqrt(diag(solve(H_iid)))
names(se_iid) <- names(opt_iid$par)

iid_alpha <- as.numeric(opt_iid$par["alpha"]); iid_alpha_se <- as.numeric(se_iid["alpha"])
iid_beta <- as.numeric(opt_iid$par[grep("^beta", names(opt_iid$par))])
iid_beta_se <- as.numeric(se_iid[grep("^beta", names(se_iid))])
iid_logTau <- as.numeric(opt_iid$par["log_tau"]) ; iid_sigU <- exp(-0.5 * iid_logTau)
iid_logTauEps <- as.numeric(opt_iid$par["log_tauEps"]) ; iid_sigEps <- exp(-0.5 * iid_logTauEps)
iid_NLL <- opt_iid$value

# Print comparison table
cat("\n\n=== Comparison: BYM2 (GH) vs IID+Nugget (DHS K=16/21) ===\n")
cat(sprintf("%-25s %25s %25s\n", "Parameter", "BYM2 (DHS)", "IID + Nug (DHS)"))
cat(paste0(rep("-", 80), "\n"))
cat(sprintf("%-25s %8.4f (%6.4f) %8.4f (%6.4f)\n", "alpha",
            res_bym2$alpha, res_bym2$alpha_se, iid_alpha, iid_alpha_se))

for(cn in canonNames) {
  # BYM2 index
  dI <- match(cn, res_bym2$covNames)
  dTxt <- if(!is.na(dI)) sprintf("%8.4f (%6.4f)", res_bym2$beta[dI], res_bym2$beta_se[dI]) else sprintf("%16s", "—")
  # IID index
  iI <- match(cn, covNamesIID)
  iTxt <- if(!is.na(iI)) sprintf("%8.4f (%6.4f)", iid_beta[iI], iid_beta_se[iI]) else sprintf("%16s", "—")
  cat(sprintf("%-25s %s %s\n", paste0("beta[", cn, "]"), dTxt, iTxt))
}
cat(sprintf("%-25s %20.4f %20.4f\n", "sigmaEps (BYM2)", res_bym2$sigEps, iid_sigEps))
cat(sprintf("%-25s %20.4f %20.4f\n", "sigmaU (IID)", NA, iid_sigU))
cat(sprintf("%-25s %20.4f %20.4f\n", "NLL", res_bym2$NLL, iid_NLL))

# Save results
save(res_bym2, list=ls(pattern="^iid_"), file="savedOutput/global/DHS_IID_vs_BYM2_16_21.RData")
cat("\nResults saved to savedOutput/global/DHS_IID_vs_BYM2_16_21.RData\n")

cat("\nDone.\n")
