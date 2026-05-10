#!/usr/bin/env Rscript
# DHS-only FE+nugget GH (Q=10), all 5 covariates, K=16/21 integration points
# Uses real DHS data (ed), modD_BYM2_GH_v2.cpp with spatial mapped out

library(TMB)
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")

Qgh = 10

# ── Data ──────────────────────────────────────────────────────────
nameTab = rbind(c("N", "n"), c("Z", "y"))
for(i in 1:nrow(nameTab)) {
  fromN = nameTab[i,1]; toN = nameTab[i,2]
  if(!(toN %in% names(ed))) ed[[toN]] = ed[[fromN]]
}

isU = ed$urban
nU = sum(isU); nR = sum(!isU)
cat("DHS clusters:", nrow(ed), "(urban:", nU, "rural:", nR, ")\n")

# ── Integration points K=16/21 ───────────────────────────────────
cat("Generating K=16/21 integration points for real DHS data...\n"); flush.console()
t0_int = proc.time()[3]
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
time_int = proc.time()[3] - t0_int
cat(sprintf("Integration points generated in %.1fs\n", time_int))

# ── Design matrices (all covariates) ─────────────────────────────
XUrb = as.matrix(intPtsDHS$covsUrb)
XRur = as.matrix(intPtsDHS$covsRur)
if("int" %in% colnames(XUrb)) XUrb = XUrb[, colnames(XUrb) != "int", drop=FALSE]
if("int" %in% colnames(XRur)) XRur = XRur[, colnames(XRur) != "int", drop=FALSE]
nBeta = ncol(XUrb)
covNames = colnames(XUrb)
cat("Covariates:", paste(covNames, collapse=", "), "\n")

ysUrbDHS = ed$y[isU]; ysRurDHS = ed$y[!isU]
nsUrbDHS = ed$n[isU]; nsRurDHS = ed$n[!isU]

# ── Area mappings ─────────────────────────────────────────────────
load("savedOutput/global/admFinalMat.RData")
bym2ArgsTMB = prepareBYM2argumentsForTMB(admFinalMat, u=0.5, alpha=2/3,
                                         constr=TRUE, scale.model=TRUE,
                                         matrixType="TsparseMatrix")
AUrbDHS = makeApointToArea(intPtsDHS$areasUrban, admFinal$NAME_FINAL)
ARurDHS = makeApointToArea(intPtsDHS$areasRural, admFinal$NAME_FINAL)
AUrbDHS = t(AUrbDHS); mode(AUrbDHS) = "numeric"
ARurDHS = t(ARurDHS); mode(ARurDHS) = "numeric"

# ── Priors ────────────────────────────────────────────────────────
alpha_pri = c(0, 100^2)
beta_pri = c(0, 10^2)
lambdaTau = getLambdaPCprec(u=1, alpha=0.1)
lambdaTauEps = getLambdaPCprec(u=1, alpha=0.1)

# ── Initial values ────────────────────────────────────────────────
initUrbP = sum(ysUrbDHS)/sum(nsUrbDHS)
initRurP = sum(ysRurDHS)/sum(nsRurDHS)
initAlpha = log(initRurP/(1-initRurP))
initBeta1 = log(initUrbP/(1-initUrbP)) - initAlpha

# ════════════════════════════════════════════════════════════════════
# GH FE-only Q=10, all covariates
# ════════════════════════════════════════════════════════════════════
cat("\n============================================================\n")
cat("  DHS-only FE+nug GH Q=10, K=16/21 — ALL covariates\n")
cat("============================================================\n\n"); flush.console()

unloadDynlibs()
if(!any(file.exists(paste0("code/modD_BYM2_GH_v2", c(".o", ".so", ".dll"))))) {
  cat("Compiling code/modD_BYM2_GH_v2.cpp...\n")
  compile("code/modD_BYM2_GH_v2.cpp", framework="TMBad", safebounds=FALSE)
}
dyn.load(dynlib("code/modD_BYM2_GH_v2"))

gh = fastGHQuad::gaussHermiteData(Qgh)
nFree = ncol(bym2ArgsTMB$Q) - 1

params_gh = list(
  log_tau=0, logit_phi=0, log_tauEps=0,
  alpha=initAlpha,
  beta=c(initBeta1, rep(0, nBeta-1)),
  w_bym2Free=rep(0, nFree),
  u_bym2Free=rep(0, nFree)
)
# Map out spatial effects and BYM2 hyperparams → FE+nugget only
map_gh = list(
  w_bym2Free=factor(rep(NA, nFree)),
  u_bym2Free=factor(rep(NA, nFree)),
  log_tau=factor(NA),
  logit_phi=factor(NA)
)

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

obj_gh = MakeADFun(data=data_gh, parameters=params_gh,
                   map=map_gh,
                   DLL="modD_BYM2_GH_v2", silent=TRUE)

t0 = proc.time()[3]
opt_gh = optim(obj_gh$par, obj_gh$fn, obj_gh$gr,
               method="BFGS", control=list(maxit=1000, reltol=1e-6))
time_gh = proc.time()[3] - t0

# Hessian-based SEs
H = optimHess(opt_gh$par, obj_gh$fn, obj_gh$gr)
se_gh = sqrt(diag(solve(H)))
names(se_gh) = names(opt_gh$par)

alpha_gh = as.numeric(opt_gh$par["alpha"])
alpha_gh_se = as.numeric(se_gh["alpha"])
beta_gh = as.numeric(opt_gh$par[grep("^beta", names(opt_gh$par))])
beta_gh_se = as.numeric(se_gh[grep("^beta", names(se_gh))])
log_tauEps_gh = as.numeric(opt_gh$par["log_tauEps"])
log_tauEps_gh_se = as.numeric(se_gh["log_tauEps"])
sigmaEps_gh = exp(-0.5 * log_tauEps_gh)

cat(sprintf("Convergence: %d | NLL: %.4f | Time: %.1fs\n",
            opt_gh$convergence, opt_gh$value, time_gh))

# ═════════════════════════════════════════════════════════════════
# Results
# ═════════════════════════════════════════════════════════════════
cat("\n============================================================\n")
cat("  DHS FE+nug GH Q=10, K=16/21, ALL covariates — Results\n")
cat("============================================================\n\n")

cat(sprintf("%-25s %10s %8s\n", "Parameter", "Estimate", "(SE)"))
cat(sprintf("%-25s %10s %8s\n", "---", "---", "---"))
cat(sprintf("%-25s %10.4f %8.4f\n", "alpha", alpha_gh, alpha_gh_se))
for(j in 1:nBeta) {
  cat(sprintf("%-25s %10.4f %8.4f\n",
              paste0("beta[", covNames[j], "]"),
              beta_gh[j], beta_gh_se[j]))
}
cat(sprintf("%-25s %10.4f\n", "sigmaEps", sigmaEps_gh))
cat(sprintf("%-25s %10.4f\n", "NLL", opt_gh$value))
cat(sprintf("%-25s %10.1f\n", "Time (s)", time_gh))

# ═════════════════════════════════════════════════════════════════
# Comparison with MICS all-cov results (from previous run)
# ═════════════════════════════════════════════════════════════════
cat("\n============================================================\n")
cat("  COMPARISON: DHS vs MICS (all 5 covariates, GH Q=10)\n")
cat("============================================================\n\n")

# MICS all-cov results (from test_FEnug_GH_allcov.R)
mics_alpha = -2.3706; mics_alpha_se = 0.1634
mics_beta = c(1.0241, -0.3645, -0.7734, 0.5032, 0.7918)
mics_beta_se = c(0.1853, 0.1469, 0.0683, 0.0580, 0.2575)
# DHS integration points use "urb"/"pop"; MICS used "urban"/"normPop"
mics_covNames = c("urban", "access", "elev", "distRiversLakes", "normPop")
# mapping: DHS name -> MICS name
dhsToMics = c(urb="urban", pop="normPop")
mics_sigmaEps = 1.6377
mics_NLL = 3280.6077

cat(sprintf("%-25s %12s %12s\n", "Parameter", "DHS K=16/21", "MICS K=100"))
cat(sprintf("%-25s %12s %12s\n", "---", "---", "---"))
cat(sprintf("%-25s %8.4f (%5.4f) %8.4f (%5.4f)\n", "alpha",
            alpha_gh, alpha_gh_se, mics_alpha, mics_alpha_se))
for(j in 1:nBeta) {
  lookupName = ifelse(covNames[j] %in% names(dhsToMics), dhsToMics[covNames[j]], covNames[j])
  matchJ = match(lookupName, mics_covNames)
  if(!is.na(matchJ)) {
    cat(sprintf("%-25s %8.4f (%5.4f) %8.4f (%5.4f)\n",
                paste0("beta[", covNames[j], "]"),
                beta_gh[j], beta_gh_se[j],
                mics_beta[matchJ], mics_beta_se[matchJ]))
  } else {
    cat(sprintf("%-25s %8.4f (%5.4f) %12s\n",
                paste0("beta[", covNames[j], "]"),
                beta_gh[j], beta_gh_se[j], "—"))
  }
}
cat(sprintf("%-25s %12.4f %12.4f\n", "sigmaEps", sigmaEps_gh, mics_sigmaEps))
cat(sprintf("%-25s %12.4f %12.4f\n", "NLL", opt_gh$value, mics_NLL))
cat(sprintf("%-25s %12.1f %12s\n", "Time (s)", time_gh, "8.0"))

cat("\nDone.\n")
