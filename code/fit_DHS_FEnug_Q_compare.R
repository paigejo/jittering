#!/usr/bin/env Rscript
# Compare FE+nugget-only (no spatial) GH model at Q=10 vs Q=25
# Uses modD_BYM2_GH_v2.cpp with BYM2/spatial parameters mapped out

library(TMB)
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")

# ── Load BYM2 simulated data ─────────────────────────────────────
load("savedOutput/simStudy1/simPopsSurveys_BYM2.RData")
ed = surveysDHS[[1]]
load("savedOutput/simStudy1/intPtsDHS_simStudy1_1.RData")

nameTab = rbind(c("N", "n"), c("Z", "y"))
for(i in 1:nrow(nameTab)) {
  fromN = nameTab[i,1]; toN = nameTab[i,2]
  if(!(toN %in% names(ed))) ed[[toN]] = ed[[fromN]]
}

# ── Build shared data objects ─────────────────────────────────────
ysUrbDHS = ed$y[ed$urban];  ysRurDHS = ed$y[!ed$urban]
nsUrbDHS = ed$n[ed$urban];  nsRurDHS = ed$n[!ed$urban]
AUrbDHS = makeApointToArea(intPtsDHS$areasUrban, admFinal$NAME_FINAL)
ARurDHS = makeApointToArea(intPtsDHS$areasRural, admFinal$NAME_FINAL)
AUrbDHS = t(AUrbDHS); mode(AUrbDHS) = "numeric"
ARurDHS = t(ARurDHS); mode(ARurDHS) = "numeric"
intPtsDHS$covsUrb = intPtsDHS$covsUrb[,-1]
intPtsDHS$covsRur = intPtsDHS$covsRur[,-1]
nBeta = ncol(intPtsDHS$covsUrb)
covNames = colnames(intPtsDHS$covsUrb)

# ── BYM2 setup (needed for Q_bym2 etc even though mapped out) ────
load("savedOutput/global/admFinalMat.RData")
bym2ArgsTMB = prepareBYM2argumentsForTMB(admFinalMat, u=0.5, alpha=2/3,
                                         constr=TRUE, scale.model=TRUE,
                                         matrixType="TsparseMatrix")
nAreas = ncol(bym2ArgsTMB$Q)
nFree = nAreas - 1

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

# ── Compile once ─────────────────────────────────────────────────
unloadDynlibs()
cat("Compiling modD_BYM2_GH_v2.cpp...\n")
compile("code/modD_BYM2_GH_v2.cpp", framework="TMBad", safebounds=FALSE)
dyn.load(dynlib("code/modD_BYM2_GH_v2"))

# ── Function to run FE+nug at given Q ────────────────────────────
run_FEnug = function(Qval) {
  gh = fastGHQuad::gaussHermiteData(Qval)
  
  data_gh = list(
    y_iUrbanDHS=ysUrbDHS, y_iRuralDHS=ysRurDHS,
    n_iUrbanDHS=nsUrbDHS, n_iRuralDHS=nsRurDHS,
    AprojUrbanDHS=AUrbDHS, AprojRuralDHS=ARurDHS,
    X_betaUrbanDHS=intPtsDHS$covsUrb, X_betaRuralDHS=intPtsDHS$covsRur,
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
  
  t0 = proc.time()[3]
  opt = optim(obj$par, obj$fn, obj$gr, method="BFGS",
              control=list(trace=1, REPORT=1, maxit=500))
  elapsed = proc.time()[3] - t0
  
  pe = opt$par
  # Compute SEs from Hessian
  H = optimHess(pe, obj$fn, obj$gr)
  se = sqrt(diag(solve(H)))
  names(se) = names(pe)
  list(opt=opt, pe=pe, se=se, time=elapsed, Q=Qval)
}

# ── Run Q=10 ─────────────────────────────────────────────────────
cat("\n============================================================\n")
cat("FE + nugget(GH) Q=10\n")
cat("============================================================\n")
res10 = run_FEnug(10)

cat(sprintf("\nQ=10: convergence=%d, NLL=%.6f, Time=%.1f s\n",
            res10$opt$convergence, res10$opt$value, res10$time))

# ── Run Q=25 ─────────────────────────────────────────────────────
cat("\n============================================================\n")
cat("FE + nugget(GH) Q=25\n")
cat("============================================================\n")
res25 = run_FEnug(25)

cat(sprintf("\nQ=25: convergence=%d, NLL=%.6f, Time=%.1f s\n",
            res25$opt$convergence, res25$opt$value, res25$time))

# ── Comparison table ─────────────────────────────────────────────
extract = function(res) {
  pe = res$pe; se = res$se
  alpha = as.numeric(pe["alpha"])
  beta = as.numeric(pe[grep("^beta", names(pe))])
  sigmaEps = exp(-0.5 * as.numeric(pe["log_tauEps"]))
  se_alpha = as.numeric(se["alpha"])
  se_beta = as.numeric(se[grep("^beta", names(se))])
  se_logTauEps = as.numeric(se["log_tauEps"])
  # Delta method: sigmaEps = exp(-0.5*log_tauEps), d/d(log_tauEps) = -0.5*sigmaEps
  se_sigmaEps = 0.5 * sigmaEps * se_logTauEps
  list(alpha=alpha, beta=beta, sigmaEps=sigmaEps,
       se_alpha=se_alpha, se_beta=se_beta, se_sigmaEps=se_sigmaEps,
       NLL=res$opt$value, time=res$time)
}

e10 = extract(res10)
e25 = extract(res25)

cat("\n============================================================\n")
cat("COMPARISON: FE+nug Q=10 vs Q=25\n")
cat("============================================================\n")
cat(sprintf("%-25s %10s %8s %10s %8s %10s\n", "Parameter", "Q=10", "(SE)", "Q=25", "(SE)", "Truth"))
cat(sprintf("%-25s %10.4f %8.4f %10.4f %8.4f %10.4f\n",
            "alpha", e10$alpha, e10$se_alpha, e25$alpha, e25$se_alpha, trueAlpha))
for(j in 1:nBeta) {
  cat(sprintf("%-25s %10.4f %8.4f %10.4f %8.4f %10.4f\n",
              paste0("beta[",covNames[j],"]"),
              e10$beta[j], e10$se_beta[j], e25$beta[j], e25$se_beta[j], trueBeta[j]))
}
cat(sprintf("%-25s %10.4f %8.4f %10.4f %8.4f %10.4f\n",
            "sigmaEps", e10$sigmaEps, e10$se_sigmaEps, e25$sigmaEps, e25$se_sigmaEps, trueSigmaEps))
cat(sprintf("%-25s %10.4f %8s %10.4f\n", "NLL", e10$NLL, "", e25$NLL))
cat(sprintf("%-25s %10.1f %8s %10.1f\n", "Time (s)", e10$time, "", e25$time))
cat(sprintf("\n|NLL difference|: %.8f\n", abs(e10$NLL - e25$NLL)))
cat("Done.\n")
