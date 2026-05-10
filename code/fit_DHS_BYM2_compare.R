#!/usr/bin/env Rscript
# Compare DHS-only models on BYM2-simulated data:
#   (A) Existing Laplace BYM2 (modBYM2JitterDHS.cpp) with nuggets as random effects
#   (B) New BYM2+GH (modD_BYM2_GH_v2.cpp) with nuggets marginalized via GH Q=10

library(TMB)
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")

# ── Load BYM2 simulated data ─────────────────────────────────────
load("savedOutput/simStudy1/simPopsSurveys_BYM2.RData")
ed = surveysDHS[[1]]
load("savedOutput/simStudy1/intPtsDHS_simStudy1_1.RData")

# Rename columns
nameTab = rbind(c("N", "n"), c("Z", "y"))
for(i in 1:nrow(nameTab)) {
  fromN = nameTab[i,1]; toN = nameTab[i,2]
  if(!(toN %in% names(ed))) ed[[toN]] = ed[[fromN]]
}

cat("DHS clusters:", nrow(ed), "(urban:", sum(ed$urban), "rural:", sum(!ed$urban), ")\n")

# ── Build shared data objects ─────────────────────────────────────
ysUrbDHS = ed$y[ed$urban]
ysRurDHS = ed$y[!ed$urban]
nsUrbDHS = ed$n[ed$urban]
nsRurDHS = ed$n[!ed$urban]

# Build Aproj: nObs x nArea (one row per obs, not per integration point)
AUrbDHS = makeApointToArea(intPtsDHS$areasUrban, admFinal$NAME_FINAL)
ARurDHS = makeApointToArea(intPtsDHS$areasRural, admFinal$NAME_FINAL)
AUrbDHS = t(AUrbDHS); mode(AUrbDHS) = "numeric"
ARurDHS = t(ARurDHS); mode(ARurDHS) = "numeric"

# Area index (0-based) for each obs
areaidxlocUrbanDHS = apply(AUrbDHS, 1, function(x) match(1, x)) - 1L
areaidxlocRuralDHS = apply(ARurDHS, 1, function(x) match(1, x)) - 1L

# Covariates: remove intercept column
intPtsDHS$covsUrb = intPtsDHS$covsUrb[,-1] # remove 'int'
intPtsDHS$covsRur = intPtsDHS$covsRur[,-1]

cat("DHS covariate columns:", colnames(intPtsDHS$covsUrb), "\n")
nBeta = ncol(intPtsDHS$covsUrb)
cat("Number of covariates:", nBeta, "\n")

# ── BYM2 setup ───────────────────────────────────────────────────
load("savedOutput/global/admFinalMat.RData")
bym2ArgsTMB = prepareBYM2argumentsForTMB(admFinalMat, u=0.5, alpha=2/3,
                                         constr=TRUE, scale.model=TRUE,
                                         matrixType="TsparseMatrix")
nAreas = ncol(bym2ArgsTMB$Q)
cat("nAreas:", nAreas, "\n")

# ── Priors ───────────────────────────────────────────────────────
# dnorm(x, mean, sd) — edDHS.R uses alpha_pri=c(0, 100^2), beta_pri=c(0, 10^2)
alpha_pri = c(0, 100^2)
beta_pri = c(0, 10^2)
lambdaTau = getLambdaPCprec(u=1, alpha=0.1)
lambdaTauEps = getLambdaPCprec(u=1, alpha=0.1)

# ── True values ──────────────────────────────────────────────────
trueAlpha = -1.25
trueBeta = c(1.00, 0, 0, 0, 0.5) # urb, access, elev, distRiversLakes, normPop
trueSigmaBYM2 = sqrt(0.5)
truePhi = 0.8
trueSigmaEps = sqrt(1.5)

# ── Initial values ───────────────────────────────────────────────
initUrbP = sum(ysUrbDHS)/sum(nsUrbDHS)
initRurP = sum(ysRurDHS)/sum(nsRurDHS)
initAlpha = log(initRurP/(1-initRurP))
initBeta1 = log(initUrbP/(1-initUrbP)) - initAlpha

# ══════════════════════════════════════════════════════════════════
# (A) Laplace BYM2 (existing modBYM2JitterDHS.cpp)
# ══════════════════════════════════════════════════════════════════
cat("\n============================================================\n")
cat("(A) DHS Laplace BYM2 (nuggets as random effects)\n")
cat("============================================================\n")

data_laplace = list(
  y_iUrbanDHS=ysUrbDHS, y_iRuralDHS=ysRurDHS,
  n_iUrbanDHS=nsUrbDHS, n_iRuralDHS=nsRurDHS,
  AprojUrbanDHS=AUrbDHS, AprojRuralDHS=ARurDHS,
  X_betaUrbanDHS=intPtsDHS$covsUrb, X_betaRuralDHS=intPtsDHS$covsRur,
  wUrbanDHS=intPtsDHS$wUrban, wRuralDHS=intPtsDHS$wRural,
  V_bym2=bym2ArgsTMB$V, Q_bym2=bym2ArgsTMB$Q,
  alpha_pri=alpha_pri, beta_pri=beta_pri,
  tr=bym2ArgsTMB$tr, gammaTildesm1=bym2ArgsTMB$gammaTildesm1,
  lambdaPhi=bym2ArgsTMB$lambda, lambdaTau=lambdaTau, lambdaTauEps=lambdaTauEps,
  options=0
)

params_laplace = list(
  alpha=initAlpha,
  beta=c(initBeta1, rep(0, nBeta-1)),
  log_tau=0, logit_phi=0, log_tauEps=0,
  Epsilon_bym2=rep(0, nAreas),
  nuggetUrbDHS=rep(0, length(ysUrbDHS)),
  nuggetRurDHS=rep(0, length(ysRurDHS))
)
rand_laplace = c('alpha', 'beta', 'Epsilon_bym2', 'nuggetUrbDHS', 'nuggetRurDHS')

unloadDynlibs()
cat("Compiling modBYM2JitterDHS.cpp...\n")
compile("code/modBYM2JitterDHS.cpp", framework="TMBad", safebounds=FALSE)
dyn.load(dynlib("code/modBYM2JitterDHS"))

cat("Building Laplace AD tape...\n"); flush.console()
obj_lap = MakeADFun(data=data_laplace, parameters=params_laplace,
                    random=rand_laplace, hessian=TRUE,
                    DLL='modBYM2JitterDHS', silent=TRUE)
TMB::config(tmbad.sparse_hessian_compress = 1)

cat("Outer params:", names(obj_lap$par), "\n")
cat("Random effects:", length(obj_lap$env$random), "params\n")
cat("Starting BFGS...\n\n"); flush.console()

# Check gradient at starting point
grad0 = obj_lap$gr(obj_lap$par)
cat("Gradient at start:", grad0, "\n"); flush.console()

t0 = proc.time()[3]
opt_lap = optim(obj_lap$par, obj_lap$fn, obj_lap$gr, method="BFGS",
                control=list(trace=1, REPORT=1, maxit=500, reltol=1e-12))
time_lap = proc.time()[3] - t0

cat(sprintf("\nLaplace convergence: %d, NLL: %.6f, Time: %.1f s\n",
            opt_lap$convergence, opt_lap$value, time_lap))
cat("Final par:", opt_lap$par, "\n")
cat("Final gradient:", obj_lap$gr(opt_lap$par), "\n"); flush.console()

cat("Computing sdreport...\n"); flush.console()
sd_lap = sdreport(obj_lap)
est_fix_lap = summary(sd_lap, "fixed")
est_ran_lap = summary(sd_lap, "random")

# alpha/beta are now inner (random) parameters
pe_r = est_ran_lap[,"Estimate"]; nms_r = rownames(est_ran_lap)
alpha_lap = pe_r[nms_r=="alpha"][1]
beta_lap = pe_r[nms_r=="beta"]
pe_f = est_fix_lap[,"Estimate"]; nms_f = rownames(est_fix_lap)

cat("\n--- Laplace DHS BYM2 ---\n")
cat(sprintf("  alpha:      %7.4f  truth: %7.4f\n", alpha_lap, trueAlpha))
covNames = colnames(intPtsDHS$covsUrb)
for(j in 1:length(beta_lap)) {
  cat(sprintf("  beta[%s]: %7.4f  truth: %7.4f\n", covNames[j], beta_lap[j], trueBeta[j]))
}
cat(sprintf("  sigmaBYM2:  %7.4f  truth: %7.4f\n",
            exp(-0.5*pe_f[nms_f=="log_tau"]), trueSigmaBYM2))
cat(sprintf("  phi:        %7.4f  truth: %7.4f\n",
            1/(1+exp(-pe_f[nms_f=="logit_phi"])), truePhi))
cat(sprintf("  sigmaEps:   %7.4f  truth: %7.4f\n",
            exp(-0.5*pe_f[nms_f=="log_tauEps"]), trueSigmaEps))

# ══════════════════════════════════════════════════════════════════
# (B) BYM2 + GH Q=10 (modD_BYM2_GH_v2.cpp)
#     Step 1: FE+nug init (map out spatial, pure optim)
#     Step 2: Full BYM2+GH warm-started from Step 1
# ══════════════════════════════════════════════════════════════════
cat("\n============================================================\n")
cat("(B) DHS BYM2+GH Q=10 (nuggets marginalized)\n")
cat("============================================================\n")

Q = 10
gh = fastGHQuad::gaussHermiteData(Q)
nFree = nAreas - 1

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

unloadDynlibs()
try(dyn.unload(dynlib("code/modBYM2JitterDHS")), silent=TRUE)
try(dyn.unload(dynlib("code/modD_BYM2_GH_v2")), silent=TRUE)
cat("Compiling modD_BYM2_GH_v2.cpp...\n")
compile("code/modD_BYM2_GH_v2.cpp", framework="TMBad", safebounds=FALSE)
dyn.load(dynlib("code/modD_BYM2_GH_v2"))

# ── Step 1: FE + nugget init (spatial mapped out, pure optimization) ──
cat("\n--- Step 1: FE + nugget(GH) initialization ---\n"); flush.console()

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

obj_init = MakeADFun(data=data_gh, parameters=params_init,
                     map=map_init,
                     DLL='modD_BYM2_GH_v2', silent=TRUE)

cat("Init params:", names(obj_init$par), "\n")
cat("Number of params:", length(obj_init$par), "\n"); flush.console()

t0_init = proc.time()[3]
opt_init = optim(obj_init$par, obj_init$fn, obj_init$gr, method="BFGS",
                 control=list(trace=1, REPORT=1, maxit=500))
time_init = proc.time()[3] - t0_init

pe_init = opt_init$par
cat(sprintf("\nFE+nug init: convergence=%d, NLL=%.4f, Time=%.1f s\n",
            opt_init$convergence, opt_init$value, time_init))
beta_init_vals = as.numeric(pe_init[grep("^beta", names(pe_init))])
cat(sprintf("  alpha:     %7.4f\n", pe_init["alpha"]))
for(j in 1:nBeta) {
  cat(sprintf("  beta[%s]: %7.4f\n", covNames[j], beta_init_vals[j]))
}
cat(sprintf("  sigmaEps:  %7.4f\n", exp(-0.5*pe_init["log_tauEps"])))
flush.console()

# ── Step 2: Full BYM2+GH warm-started from init ──
cat("\n--- Step 2: Full BYM2+GH (warm-started) ---\n"); flush.console()

params_gh = list(
  log_tau=0, logit_phi=0,
  log_tauEps=as.numeric(pe_init["log_tauEps"]),
  alpha=as.numeric(pe_init["alpha"]),
  beta=as.numeric(pe_init[grep("^beta", names(pe_init))]),
  w_bym2Free=rep(0, nFree),
  u_bym2Free=rep(0, nFree)
)
rand_gh = c('alpha', 'beta', 'w_bym2Free', 'u_bym2Free')

obj_gh = MakeADFun(data=data_gh, parameters=params_gh,
                   random=rand_gh,
                   DLL='modD_BYM2_GH_v2', silent=TRUE)
TMB::config(tmbad.sparse_hessian_compress = 1)

cat("Outer params:", names(obj_gh$par), "\n")
cat("Random effects:", length(obj_gh$env$random), "params\n")
cat("Starting BFGS...\n\n"); flush.console()

t0 = proc.time()[3]
opt_gh = optim(obj_gh$par, obj_gh$fn, obj_gh$gr, method="BFGS",
               control=list(trace=1, REPORT=1, maxit=500))
time_gh = proc.time()[3] - t0

cat(sprintf("\nGH convergence: %d, NLL: %.4f, Time: %.1f s\n",
            opt_gh$convergence, opt_gh$value, time_gh))

cat("Computing sdreport...\n"); flush.console()
sd_gh = sdreport(obj_gh)
est_fix_gh = summary(sd_gh, "fixed")
est_ran_gh = summary(sd_gh, "random")

pe_r2 = est_ran_gh[,"Estimate"]; nms_r2 = rownames(est_ran_gh)
alpha_gh = pe_r2[nms_r2=="alpha"][1]
beta_gh = pe_r2[nms_r2=="beta"]
pe_f2 = est_fix_gh[,"Estimate"]; nms_f2 = rownames(est_fix_gh)

cat("\n--- GH DHS BYM2 ---\n")
cat(sprintf("  alpha:      %7.4f  truth: %7.4f\n", alpha_gh, trueAlpha))
for(j in 1:length(beta_gh)) {
  cat(sprintf("  beta[%s]: %7.4f  truth: %7.4f\n", covNames[j], beta_gh[j], trueBeta[j]))
}
cat(sprintf("  sigmaBYM2:  %7.4f  truth: %7.4f\n",
            exp(-0.5*pe_f2[nms_f2=="log_tau"]), trueSigmaBYM2))
cat(sprintf("  phi:        %7.4f  truth: %7.4f\n",
            1/(1+exp(-pe_f2[nms_f2=="logit_phi"])), truePhi))
cat(sprintf("  sigmaEps:   %7.4f  truth: %7.4f\n",
            exp(-0.5*pe_f2[nms_f2=="log_tauEps"]), trueSigmaEps))

# ══════════════════════════════════════════════════════════════════
# COMPARISON
# ══════════════════════════════════════════════════════════════════
# FE+nug init results
alpha_init = as.numeric(pe_init["alpha"])
beta_init = as.numeric(pe_init[grep("^beta", names(pe_init))])
sigmaEps_init = exp(-0.5*as.numeric(pe_init["log_tauEps"]))

cat("\n============================================================\n")
cat("COMPARISON: Laplace vs FE+nug init vs GH (DHS-only BYM2)\n")
cat("============================================================\n")
cat(sprintf("%-20s %10s %10s %10s %10s\n", "Parameter", "Laplace", "FE+nug", "GH", "Truth"))
cat(sprintf("%-20s %10.4f %10.4f %10.4f %10.4f\n", "alpha", alpha_lap, alpha_init, alpha_gh, trueAlpha))
for(j in 1:length(beta_lap)) {
  cat(sprintf("%-20s %10.4f %10.4f %10.4f %10.4f\n",
              paste0("beta[",covNames[j],"]"), beta_lap[j], beta_init[j], beta_gh[j], trueBeta[j]))
}
cat(sprintf("%-20s %10.4f %10s %10.4f %10.4f\n", "sigmaBYM2",
            exp(-0.5*pe_f[nms_f=="log_tau"]), "---", exp(-0.5*pe_f2[nms_f2=="log_tau"]), trueSigmaBYM2))
cat(sprintf("%-20s %10.4f %10s %10.4f %10.4f\n", "phi",
            1/(1+exp(-pe_f[nms_f=="logit_phi"])), "---", 1/(1+exp(-pe_f2[nms_f2=="logit_phi"])), truePhi))
cat(sprintf("%-20s %10.4f %10.4f %10.4f %10.4f\n", "sigmaEps",
            exp(-0.5*pe_f[nms_f=="log_tauEps"]), sigmaEps_init, exp(-0.5*pe_f2[nms_f2=="log_tauEps"]), trueSigmaEps))
cat(sprintf("\n%-20s %10.1f %10.1f %10.1f\n", "Time (s)", time_lap, time_init, time_gh))

# ── Save ──────────────────────────────────────────────────────────
results_DHS = list(
  laplace=list(opt=opt_lap, sd=sd_lap, time=time_lap),
  gh_init=list(opt=opt_init, time=time_init),
  gh=list(opt=opt_gh, sd=sd_gh, time=time_gh, Q=Q)
)
save(results_DHS, file="savedOutput/fitDHS_BYM2_comparison.RData")
cat("\nSaved to savedOutput/fitDHS_BYM2_comparison.RData\n")
cat("Done.\n")
