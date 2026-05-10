#!/usr/bin/env Rscript
# Fit DHS FE+nugget GH model at TRUE cluster coordinates (K=1)
# to test whether biased estimates are due to integration points or something else.
# Compares: true-coords (K=1) vs integration points (K=11/16)

library(TMB)
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")

# ── Load BYM2 simulated data ─────────────────────────────────────
load("savedOutput/simStudy1/simPopsSurveys_BYM2.RData")
ed = surveysDHS[[1]]
nameTab = rbind(c("N", "n"), c("Z", "y"))
for(i in 1:nrow(nameTab)) {
  fromN = nameTab[i,1]; toN = nameTab[i,2]
  if(!(toN %in% names(ed))) ed[[toN]] = ed[[fromN]]
}

isU = ed$urban
nU = sum(isU); nR = sum(!isU)
cat("DHS clusters:", nrow(ed), "(urban:", nU, "rural:", nR, ")\n")

# ── Extract covariates at TRUE cluster coordinates ────────────────
cat("Extracting covariates at true cluster (lon, lat)...\n")
trueCoords = cbind(ed$lon, ed$lat)
covMat = getDesignMatPopNorm(trueCoords, useThreshPopMat=TRUE)
cat("covMat columns:", colnames(covMat), "\n")
cat("covMat dim:", nrow(covMat), "x", ncol(covMat), "\n")

# Build X_beta matrices (nObs x nBeta, K=1 so no stacking needed)
# covMat from getDesignMatPopNorm has: int, pop, urb, access, elev, distRiversLakes, urbanicity
# intPtsDHS$covsUrb (before removing int) has: int, urb, access, elev, distRiversLakes, pop
# Match the intPtsDHS column order: urb, access, elev, distRiversLakes, pop
keepCols = c("urb", "access", "elev", "distRiversLakes", "pop")
covMat_sel = covMat[, keepCols, drop=FALSE]
cat("Selected covariate columns:", colnames(covMat_sel), "\n")
cat("Any NAs?:", any(is.na(covMat_sel)), "\n")
if(any(is.na(covMat_sel))) {
  cat("NA counts per column:", colSums(is.na(covMat_sel)), "\n")
  # Fill NAs with 0 (normalized covariates, 0 = mean)
  covMat_sel[is.na(covMat_sel)] = 0
}
XUrb_true = as.matrix(covMat_sel[isU, , drop=FALSE])
XRur_true = as.matrix(covMat_sel[!isU, , drop=FALSE])
cat("XUrb_true dim:", nrow(XUrb_true), "x", ncol(XUrb_true), "\n")
cat("XRur_true dim:", nrow(XRur_true), "x", ncol(XRur_true), "\n")

covNames = colnames(XUrb_true)
nBeta = ncol(XUrb_true)
cat("nBeta:", nBeta, "covNames:", covNames, "\n")

# Weights: K=1, weight=1
wUrb_true = matrix(1, nrow=nU, ncol=1)
wRur_true = matrix(1, nrow=nR, ncol=1)

# Obs data
ysUrbDHS = ed$y[isU]; ysRurDHS = ed$y[!isU]
nsUrbDHS = ed$n[isU]; nsRurDHS = ed$n[!isU]

# Aproj: nObs x nArea, each obs maps to its true area
# Construct from ed$area (admin1) mapped to admFinal
load("savedOutput/global/admFinalMat.RData")
bym2ArgsTMB = prepareBYM2argumentsForTMB(admFinalMat, u=0.5, alpha=2/3,
                                         constr=TRUE, scale.model=TRUE,
                                         matrixType="TsparseMatrix")
nAreas = ncol(bym2ArgsTMB$Q)
nFree = nAreas - 1

# Area assignment: ed$subarea -> admFinal$NAME_FINAL
areaNames = admFinal$NAME_FINAL
AUrbDHS_true = makeApointToArea(ed$subarea[isU], areaNames)
ARurDHS_true = makeApointToArea(ed$subarea[!isU], areaNames)
AUrbDHS_true = t(AUrbDHS_true); mode(AUrbDHS_true) = "numeric"
ARurDHS_true = t(ARurDHS_true); mode(ARurDHS_true) = "numeric"
cat("AUrb dim:", dim(AUrbDHS_true), "ARur dim:", dim(ARurDHS_true), "\n")

# ── Priors ────────────────────────────────────────────────────────
alpha_pri = c(0, 100^2)
beta_pri = c(0, 10^2)
lambdaTau = getLambdaPCprec(u=1, alpha=0.1)
lambdaTauEps = getLambdaPCprec(u=1, alpha=0.1)

# ── True values ──────────────────────────────────────────────────
trueAlpha = -1.25
trueBeta = c(1.00, 0, 0, 0, 0.5)
trueSigmaEps = sqrt(1.5)

# ── Initial values ───────────────────────────────────────────────
initUrbP = sum(ysUrbDHS)/sum(nsUrbDHS)
initRurP = sum(ysRurDHS)/sum(nsRurDHS)
initAlpha = log(initRurP/(1-initRurP))
initBeta1 = log(initUrbP/(1-initUrbP)) - initAlpha

# ── Compile ──────────────────────────────────────────────────────
unloadDynlibs()
cat("Compiling modD_BYM2_GH_v2.cpp...\n")
compile("code/modD_BYM2_GH_v2.cpp", framework="TMBad", safebounds=FALSE)
dyn.load(dynlib("code/modD_BYM2_GH_v2"))

# ── Function to run FE+nug at given Q, given data ────────────────
run_FEnug = function(Qval, XUrb, XRur, wUrb, wRur, AUrb, ARur, label) {
  gh = fastGHQuad::gaussHermiteData(Qval)
  
  data_gh = list(
    y_iUrbanDHS=ysUrbDHS, y_iRuralDHS=ysRurDHS,
    n_iUrbanDHS=nsUrbDHS, n_iRuralDHS=nsRurDHS,
    AprojUrbanDHS=AUrb, AprojRuralDHS=ARur,
    X_betaUrbanDHS=XUrb, X_betaRuralDHS=XRur,
    wUrbanDHS=wUrb, wRuralDHS=wRur,
    Q_bym2=bym2ArgsTMB$Q,
    lchoose_urban=lchoose(nsUrbDHS, ysUrbDHS),
    lchoose_rural=lchoose(nsRurDHS, ysRurDHS),
    gh_nodes=gh$x, gh_weights=gh$w,
    alpha_pri=alpha_pri, beta_pri=beta_pri,
    tr=bym2ArgsTMB$tr, gammaTildesm1=bym2ArgsTMB$gammaTildesm1,
    lambdaPhi=bym2ArgsTMB$lambda, lambdaTau=lambdaTau, lambdaTauEps=lambdaTauEps,
    options=0
  )
  
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
  
  obj = MakeADFun(data=data_gh, parameters=params_init,
                  map=map_init,
                  DLL='modD_BYM2_GH_v2', silent=TRUE)
  
  # Check initial value
  fn0 = obj$fn(obj$par)
  cat(sprintf("Initial NLL: %.6f\n", fn0))
  if(!is.finite(fn0)) {
    cat("WARNING: Initial NLL is not finite! Checking params...\n")
    cat("obj$par:", obj$par, "\n")
    return(NULL)
  }
  
  t0 = proc.time()[3]
  opt = optim(obj$par, obj$fn, obj$gr, method="BFGS",
              control=list(trace=1, REPORT=1, maxit=500))
  elapsed = proc.time()[3] - t0
  
  pe = opt$par
  H = optimHess(pe, obj$fn, obj$gr)
  se = sqrt(diag(solve(H)))
  names(se) = names(pe)
  
  list(opt=opt, pe=pe, se=se, time=elapsed, Q=Qval, label=label)
}

# ── Also load integration point data for comparison ──────────────
load("savedOutput/simStudy1/intPtsDHS_simStudy1_1.RData")
intPtsDHS$covsUrb = intPtsDHS$covsUrb[,-1]  # remove intercept
intPtsDHS$covsRur = intPtsDHS$covsRur[,-1]

AUrbDHS_int = makeApointToArea(intPtsDHS$areasUrban, areaNames)
ARurDHS_int = makeApointToArea(intPtsDHS$areasRural, areaNames)
AUrbDHS_int = t(AUrbDHS_int); mode(AUrbDHS_int) = "numeric"
ARurDHS_int = t(ARurDHS_int); mode(ARurDHS_int) = "numeric"

# ── Run: TRUE coordinates (K=1, Q=25) ────────────────────────────
cat("\n============================================================\n")
cat("FE+nug GH Q=25: TRUE COORDINATES (K=1)\n")
cat("============================================================\n")
res_true = run_FEnug(25, XUrb_true, XRur_true, wUrb_true, wRur_true,
                     AUrbDHS_true, ARurDHS_true, "TrueCoords")

cat(sprintf("\nTrueCoords: convergence=%d, NLL=%.6f, Time=%.1f s\n",
            res_true$opt$convergence, res_true$opt$value, res_true$time))

# ── Run: Integration points (K=11/16, Q=25) ──────────────────────
cat("\n============================================================\n")
cat("FE+nug GH Q=25: INTEGRATION POINTS (K=11/16)\n")
cat("============================================================\n")
res_int = run_FEnug(25, intPtsDHS$covsUrb, intPtsDHS$covsRur,
                    intPtsDHS$wUrban, intPtsDHS$wRural,
                    AUrbDHS_int, ARurDHS_int, "IntPts")

cat(sprintf("\nIntPts: convergence=%d, NLL=%.6f, Time=%.1f s\n",
            res_int$opt$convergence, res_int$opt$value, res_int$time))

# ── Comparison table ─────────────────────────────────────────────
extract = function(res) {
  pe = res$pe; se = res$se
  alpha = as.numeric(pe["alpha"])
  beta = as.numeric(pe[grep("^beta", names(pe))])
  sigmaEps = exp(-0.5 * as.numeric(pe["log_tauEps"]))
  se_alpha = as.numeric(se["alpha"])
  se_beta = as.numeric(se[grep("^beta", names(se))])
  se_logTauEps = as.numeric(se["log_tauEps"])
  se_sigmaEps = 0.5 * sigmaEps * se_logTauEps
  list(alpha=alpha, beta=beta, sigmaEps=sigmaEps,
       se_alpha=se_alpha, se_beta=se_beta, se_sigmaEps=se_sigmaEps,
       NLL=res$opt$value, time=res$time)
}

eT = extract(res_true)
eI = extract(res_int)

cat("\n============================================================\n")
cat("COMPARISON: True Coords vs Integration Points (FE+nug, Q=25)\n")
cat("============================================================\n")
cat(sprintf("%-25s %10s %8s %10s %8s %10s\n", "Parameter", "TrueCoord", "(SE)", "IntPts", "(SE)", "Truth"))
cat(sprintf("%-25s %10.4f %8.4f %10.4f %8.4f %10.4f\n",
            "alpha", eT$alpha, eT$se_alpha, eI$alpha, eI$se_alpha, trueAlpha))
for(j in 1:nBeta) {
  cat(sprintf("%-25s %10.4f %8.4f %10.4f %8.4f %10.4f\n",
              paste0("beta[",covNames[j],"]"),
              eT$beta[j], eT$se_beta[j], eI$beta[j], eI$se_beta[j], trueBeta[j]))
}
cat(sprintf("%-25s %10.4f %8.4f %10.4f %8.4f %10.4f\n",
            "sigmaEps", eT$sigmaEps, eT$se_sigmaEps, eI$sigmaEps, eI$se_sigmaEps, trueSigmaEps))
cat(sprintf("%-25s %10.4f %8s %10.4f\n", "NLL", eT$NLL, "", eI$NLL))
cat(sprintf("%-25s %10.1f %8s %10.1f\n", "Time (s)", eT$time, "", eI$time))

# Show how many SEs from truth
cat("\n--- Distance from truth (in SEs) ---\n")
cat(sprintf("%-25s %10s %10s\n", "Parameter", "TrueCoord", "IntPts"))
cat(sprintf("%-25s %10.2f %10.2f\n", "alpha",
            abs(eT$alpha - trueAlpha)/eT$se_alpha,
            abs(eI$alpha - trueAlpha)/eI$se_alpha))
for(j in 1:nBeta) {
  if(trueBeta[j] != 0) {
    cat(sprintf("%-25s %10.2f %10.2f\n", paste0("beta[",covNames[j],"]"),
                abs(eT$beta[j] - trueBeta[j])/eT$se_beta[j],
                abs(eI$beta[j] - trueBeta[j])/eI$se_beta[j]))
  }
}
cat(sprintf("%-25s %10.2f %10.2f\n", "sigmaEps",
            abs(eT$sigmaEps - trueSigmaEps)/eT$se_sigmaEps,
            abs(eI$sigmaEps - trueSigmaEps)/eI$se_sigmaEps))

cat("\nDone.\n")
