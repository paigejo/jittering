#!/usr/bin/env Rscript
# Compare FE+nug GH at Q=10:
#   (A) True coordinates (K=1)
#   (B) Original integration points (K=11 urban, K=16 rural)
#   (C) More integration points (K=31 urban, K=41 rural)
#       10 pts/ring, 3 inner + 1 outer (rural only)

library(TMB)
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")

Qval = 10

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

# ── True values ──────────────────────────────────────────────────
trueAlpha = -1.25
trueBeta = c(1.00, 0, 0, 0, 0.5)
trueSigmaEps = sqrt(1.5)
covNames = c("urb", "access", "elev", "distRiversLakes", "pop")

# ── BYM2 setup (for mapped-out spatial effect) ───────────────────
load("savedOutput/global/admFinalMat.RData")
bym2ArgsTMB = prepareBYM2argumentsForTMB(admFinalMat, u=0.5, alpha=2/3,
                                         constr=TRUE, scale.model=TRUE,
                                         matrixType="TsparseMatrix")
nAreas = ncol(bym2ArgsTMB$Q)
nFree = nAreas - 1
areaNames = admFinal$NAME_FINAL

# ── Obs data ─────────────────────────────────────────────────────
ysUrbDHS = ed$y[isU]; ysRurDHS = ed$y[!isU]
nsUrbDHS = ed$n[isU]; nsRurDHS = ed$n[!isU]
nBeta = 5

# ── Initial values ───────────────────────────────────────────────
initUrbP = sum(ysUrbDHS)/sum(nsUrbDHS)
initRurP = sum(ysRurDHS)/sum(nsRurDHS)
initAlpha = log(initRurP/(1-initRurP))
initBeta1 = log(initUrbP/(1-initUrbP)) - initAlpha

# ── Priors ────────────────────────────────────────────────────────
alpha_pri = c(0, 100^2)
beta_pri = c(0, 10^2)
lambdaTau = getLambdaPCprec(u=1, alpha=0.1)
lambdaTauEps = getLambdaPCprec(u=1, alpha=0.1)

# ── Compile ──────────────────────────────────────────────────────
unloadDynlibs()
cat("Compiling modD_BYM2_GH_v2.cpp...\n")
compile("code/modD_BYM2_GH_v2.cpp", framework="TMBad", safebounds=FALSE)
dyn.load(dynlib("code/modD_BYM2_GH_v2"))

# ── FE+nug fitting function ──────────────────────────────────────
run_FEnug = function(Q, XUrb, XRur, wUrb, wRur, AUrb, ARur, label) {
  gh = fastGHQuad::gaussHermiteData(Q)
  
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
  
  fn0 = obj$fn(obj$par)
  cat(sprintf("[%s] Initial NLL: %.6f\n", label, fn0))
  if(!is.finite(fn0)) {
    cat("WARNING: Initial NLL is not finite!\n")
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
  
  list(opt=opt, pe=pe, se=se, time=elapsed, Q=Q, label=label)
}

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

# ══════════════════════════════════════════════════════════════════
# (A) TRUE COORDINATES (K=1)
# ══════════════════════════════════════════════════════════════════
cat("\n=== Preparing TRUE COORDINATES (K=1) ===\n")
trueCoords = cbind(ed$lon, ed$lat)
covMat = getDesignMatPopNorm(trueCoords, useThreshPopMat=TRUE)
keepCols = c("urb", "access", "elev", "distRiversLakes", "pop")
covMat_sel = covMat[, keepCols, drop=FALSE]
covMat_sel[is.na(covMat_sel)] = 0
XUrb_true = as.matrix(covMat_sel[isU, , drop=FALSE])
XRur_true = as.matrix(covMat_sel[!isU, , drop=FALSE])
wUrb_true = matrix(1, nrow=nU, ncol=1)
wRur_true = matrix(1, nrow=nR, ncol=1)
AUrb_true = makeApointToArea(ed$subarea[isU], areaNames)
ARur_true = makeApointToArea(ed$subarea[!isU], areaNames)
AUrb_true = t(AUrb_true); mode(AUrb_true) = "numeric"
ARur_true = t(ARur_true); mode(ARur_true) = "numeric"

res_true = run_FEnug(Qval, XUrb_true, XRur_true, wUrb_true, wRur_true,
                     AUrb_true, ARur_true, "TrueCoords_K1")
cat(sprintf("\nTrueCoords: conv=%d, NLL=%.4f, Time=%.1fs\n",
            res_true$opt$convergence, res_true$opt$value, res_true$time))

# ══════════════════════════════════════════════════════════════════
# (B) ORIGINAL INTEGRATION POINTS (K=11/16)
# ══════════════════════════════════════════════════════════════════
cat("\n=== Loading ORIGINAL integration points (K=11/16) ===\n")
load("savedOutput/simStudy1/intPtsDHS_simStudy1_1.RData")
intOrig = intPtsDHS
intOrig$covsUrb = intOrig$covsUrb[,-1]
intOrig$covsRur = intOrig$covsRur[,-1]
AUrb_orig = makeApointToArea(intOrig$areasUrban, areaNames)
ARur_orig = makeApointToArea(intOrig$areasRural, areaNames)
AUrb_orig = t(AUrb_orig); mode(AUrb_orig) = "numeric"
ARur_orig = t(ARur_orig); mode(ARur_orig) = "numeric"
cat("Original: wUrban", dim(intOrig$wUrban), "wRural", dim(intOrig$wRural), "\n")

res_orig = run_FEnug(Qval, intOrig$covsUrb, intOrig$covsRur,
                     intOrig$wUrban, intOrig$wRural,
                     AUrb_orig, ARur_orig, "IntPts_K11_16")
cat(sprintf("\nOriginal: conv=%d, NLL=%.4f, Time=%.1fs\n",
            res_orig$opt$convergence, res_orig$opt$value, res_orig$time))

# ══════════════════════════════════════════════════════════════════
# (C) MORE INTEGRATION POINTS (K=16 urban, K=21 rural)
#     5 pts/ring, JInner=4 (center+3 inner), JOuter=0 urban / 1 rural
# ══════════════════════════════════════════════════════════════════
cat("\n=== Generating MORE integration points (K=16/21) ===\n")
cat("Urban: 1 + 3*5 = 16, JInner=4, JOuter=0\n")
cat("Rural: 1 + 3*5 + 1*5 = 21, JInner=4, JOuter=1\n")

t0_gen = proc.time()[3]
intMore = makeAllIntegrationPointsDHS(
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
t_gen = proc.time()[3] - t0_gen
cat(sprintf("Integration point generation took %.1f s\n", t_gen))
cat("New: wUrban", dim(intMore$wUrban), "wRural", dim(intMore$wRural), "\n")
cat("New: covsUrb", dim(intMore$covsUrb), "covsRur", dim(intMore$covsRur), "\n")

# Save the new integration points for reuse
intPtsDHS = intMore
save(intPtsDHS, file="savedOutput/simStudy1/intPtsDHS_simStudy1_1_K16_21.RData")
cat("Saved to savedOutput/simStudy1/intPtsDHS_simStudy1_1_K16_21.RData\n")

intMore$covsUrb = intMore$covsUrb[,-1]
intMore$covsRur = intMore$covsRur[,-1]
AUrb_more = makeApointToArea(intMore$areasUrban, areaNames)
ARur_more = makeApointToArea(intMore$areasRural, areaNames)
AUrb_more = t(AUrb_more); mode(AUrb_more) = "numeric"
ARur_more = t(ARur_more); mode(ARur_more) = "numeric"

res_more = run_FEnug(Qval, intMore$covsUrb, intMore$covsRur,
                     intMore$wUrban, intMore$wRural,
                     AUrb_more, ARur_more, "IntPts_K16_21")
cat(sprintf("\nMore IntPts: conv=%d, NLL=%.4f, Time=%.1fs\n",
            res_more$opt$convergence, res_more$opt$value, res_more$time))

# ══════════════════════════════════════════════════════════════════
# COMPARISON TABLE
# ══════════════════════════════════════════════════════════════════
eT = extract(res_true)
eO = extract(res_orig)
eM = extract(res_more)

cat("\n================================================================\n")
cat(sprintf("FE+nug GH Q=%d: True Coords vs K=11/16 vs K=16/21\n", Qval))
cat("================================================================\n")
cat(sprintf("%-20s %9s %7s %9s %7s %9s %7s %9s\n",
            "Parameter", "TrueCrd", "(SE)", "K=11/16", "(SE)", "K=16/21", "(SE)", "Truth"))

cat(sprintf("%-20s %9.4f %7.4f %9.4f %7.4f %9.4f %7.4f %9.4f\n",
            "alpha", eT$alpha, eT$se_alpha, eO$alpha, eO$se_alpha, eM$alpha, eM$se_alpha, trueAlpha))
for(j in 1:nBeta) {
  cat(sprintf("%-20s %9.4f %7.4f %9.4f %7.4f %9.4f %7.4f %9.4f\n",
              paste0("beta[",covNames[j],"]"),
              eT$beta[j], eT$se_beta[j], eO$beta[j], eO$se_beta[j],
              eM$beta[j], eM$se_beta[j], trueBeta[j]))
}
cat(sprintf("%-20s %9.4f %7.4f %9.4f %7.4f %9.4f %7.4f %9.4f\n",
            "sigmaEps", eT$sigmaEps, eT$se_sigmaEps, eO$sigmaEps, eO$se_sigmaEps,
            eM$sigmaEps, eM$se_sigmaEps, trueSigmaEps))
cat(sprintf("%-20s %9.4f %7s %9.4f %7s %9.4f\n", "NLL", eT$NLL, "", eO$NLL, "", eM$NLL))
cat(sprintf("%-20s %9.1f %7s %9.1f %7s %9.1f\n", "Time (s)", eT$time, "", eO$time, "", eM$time))

cat("\n--- Distance from truth (in SEs) ---\n")
cat(sprintf("%-20s %9s %9s %9s\n", "Parameter", "TrueCrd", "K=11/16", "K=16/21"))
cat(sprintf("%-20s %9.2f %9.2f %9.2f\n", "alpha",
            abs(eT$alpha - trueAlpha)/eT$se_alpha,
            abs(eO$alpha - trueAlpha)/eO$se_alpha,
            abs(eM$alpha - trueAlpha)/eM$se_alpha))
for(j in 1:nBeta) {
  if(trueBeta[j] != 0) {
    cat(sprintf("%-20s %9.2f %9.2f %9.2f\n", paste0("beta[",covNames[j],"]"),
                abs(eT$beta[j] - trueBeta[j])/eT$se_beta[j],
                abs(eO$beta[j] - trueBeta[j])/eO$se_beta[j],
                abs(eM$beta[j] - trueBeta[j])/eM$se_beta[j]))
  }
}
cat(sprintf("%-20s %9.2f %9.2f %9.2f\n", "sigmaEps",
            abs(eT$sigmaEps - trueSigmaEps)/eT$se_sigmaEps,
            abs(eO$sigmaEps - trueSigmaEps)/eO$se_sigmaEps,
            abs(eM$sigmaEps - trueSigmaEps)/eM$se_sigmaEps))

cat("\nDone.\n")
