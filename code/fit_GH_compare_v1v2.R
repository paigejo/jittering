#!/usr/bin/env Rscript
# Compare v1 (triple-loop) vs v2 (matrix-ops) GH templates
# Config A only: FE + nugget(GH) Q=10, no random effects

library(TMB)
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
source("code/modM_MIIDonly.R")

# ── Load BYM2 simulated data ────────────────────────────────────
load("savedOutput/simStudy1/simPopsSurveys_BYM2.RData")
edMICS = surveysMICS[[1]]

nameTab = rbind(c("N", "ns"), c("N", "n"), c("Z", "ys"), c("Z", "y"))
for(i in 1:nrow(nameTab)) {
  fromN = nameTab[i,1]; toN = nameTab[i,2]
  if(!(toN %in% names(edMICS))) edMICS[[toN]] = edMICS[[fromN]]
  if(!(toN %in% names(ed))) ed[[toN]] = ed[[fromN]]
}
if(!("Stratum" %in% names(edMICS)))
  edMICS$Stratum = adm2ToStratumMICS(edMICS$subarea)

inputsMDM = makeInputsMDM(datDHS=ed, datMICS=edMICS,
                          KMICS=100, admMICS=admFinal, adm2DHS=adm2Full)

allCovNames = colnames(inputsMDM$intPtsMICS$XUrb)
keepIdx = which(allCovNames %in% c("urban", "normPop"))
inputsMDM$intPtsMICS$XUrb = inputsMDM$intPtsMICS$XUrb[, keepIdx, drop=FALSE]
inputsMDM$intPtsMICS$XRur = inputsMDM$intPtsMICS$XRur[, keepIdx, drop=FALSE]
cat("X columns:", colnames(inputsMDM$intPtsMICS$XUrb), "\n")

thisEnv = environment()
list2env(inputsMDM, envir=thisEnv)

lambdaTau = getLambdaPCprec(u=1, alpha=0.1)
lambdaTauEps = getLambdaPCprec(u=1, alpha=0.1)
nAreas = max(c(areaidxlocUrbanMICS, areaidxlocRuralMICS)) + 1

Q = 10
gh = fastGHQuad::gaussHermiteData(Q)

# ── Initial values ───────────────────────────────────────────────
initUrbP = sum(ysUrbMICS)/sum(nsUrbMICS)
initRurP = sum(ysRurMICS)/sum(nsRurMICS)
initAlpha = log(initRurP/(1-initRurP))
initBeta1 = log(initUrbP/(1-initUrbP)) - initAlpha

# ── Shared data for v1 ──────────────────────────────────────────
data_v1 = list(
  y_iUrbanMICS=ysUrbMICS, y_iRuralMICS=ysRurMICS,
  n_iUrbanMICS=nsUrbMICS, n_iRuralMICS=nsRurMICS,
  areaidxlocUrbanMICS=areaidxlocUrbanMICS,
  areaidxlocRuralMICS=areaidxlocRuralMICS,
  X_betaUrbanMICS=intPtsMICS$XUrb,
  X_betaRuralMICS=intPtsMICS$XRur,
  wUrbanMICS=intPtsMICS$wUrban,
  wRuralMICS=intPtsMICS$wRural,
  nAreas=nAreas,
  gh_nodes=gh$x, gh_weights=gh$w,
  alpha_pri=c(0, 100^2), beta_pri=c(0, sqrt(1000)),
  lambdaTau=lambdaTau, lambdaTauEps=lambdaTauEps
)

# ── Data for v2: same + lchoose ──────────────────────────────────
data_v2 = c(data_v1, list(
  lchoose_urban = lchoose(nsUrbMICS, ysUrbMICS),
  lchoose_rural = lchoose(nsRurMICS, ysRurMICS)
))

params_A = list(
  log_tau=0, log_tauEps=0,
  alpha=initAlpha, beta_urban=initBeta1, beta_normPop=0,
  u_spatial=rep(0, nAreas)
)
mapA = list(u_spatial=factor(rep(NA, nAreas)), log_tau=factor(NA))

# ══════════════════════════════════════════════════════════════════
# V1: Original triple-loop
# ══════════════════════════════════════════════════════════════════
cat("\n=== Compiling v1 (original) ===\n")
unloadDynlibs()
compile("code/modM_MIIDonly_GH.cpp", framework="TMBad", safebounds=FALSE)
dyn.load(dynlib("code/modM_MIIDonly_GH"))

obj_v1 = MakeADFun(data=data_v1, parameters=params_A,
                    map=mapA, DLL='modM_MIIDonly_GH', silent=TRUE)

cat("V1 initial NLL:", obj_v1$fn(obj_v1$par), "\n")

t0 = proc.time()[3]
opt_v1 = optim(obj_v1$par, obj_v1$fn, obj_v1$gr, method="BFGS",
               control=list(trace=1, REPORT=1, maxit=500))
time_v1 = proc.time()[3] - t0

cat(sprintf("\nV1 time: %.2f s, NLL: %.4f, convergence: %d\n",
            time_v1, opt_v1$value, opt_v1$convergence))

sd_v1 = sdreport(obj_v1)
est_v1 = summary(sd_v1, "fixed")
cat("V1 estimates:\n"); print(est_v1)

# ══════════════════════════════════════════════════════════════════
# V2: Matrix-ops optimized
# ══════════════════════════════════════════════════════════════════
cat("\n=== Compiling v2 (optimized) ===\n")
unloadDynlibs()
compile("code/modM_MIIDonly_GH_v2.cpp", framework="TMBad", safebounds=FALSE)
dyn.load(dynlib("code/modM_MIIDonly_GH_v2"))

obj_v2 = MakeADFun(data=data_v2, parameters=params_A,
                    map=mapA, DLL='modM_MIIDonly_GH_v2', silent=TRUE)

cat("V2 initial NLL:", obj_v2$fn(obj_v2$par), "\n")

t0 = proc.time()[3]
opt_v2 = optim(obj_v2$par, obj_v2$fn, obj_v2$gr, method="BFGS",
               control=list(trace=1, REPORT=1, maxit=500))
time_v2 = proc.time()[3] - t0

cat(sprintf("\nV2 time: %.2f s, NLL: %.4f, convergence: %d\n",
            time_v2, opt_v2$value, opt_v2$convergence))

sd_v2 = sdreport(obj_v2)
est_v2 = summary(sd_v2, "fixed")
cat("V2 estimates:\n"); print(est_v2)

# ══════════════════════════════════════════════════════════════════
# Comparison
# ══════════════════════════════════════════════════════════════════
cat("\n============================================================\n")
cat("COMPARISON: V1 (loop) vs V2 (matrix)\n")
cat("============================================================\n")

# Extract estimates
p1 = opt_v1$par; p2 = opt_v2$par
cat(sprintf("  log_tauEps: v1=%.6f  v2=%.6f  diff=%.2e\n",
            p1["log_tauEps"], p2["log_tauEps"],
            abs(p1["log_tauEps"]-p2["log_tauEps"])))
cat(sprintf("  alpha:      v1=%.6f  v2=%.6f  diff=%.2e\n",
            p1["alpha"], p2["alpha"],
            abs(p1["alpha"]-p2["alpha"])))
cat(sprintf("  beta_urban: v1=%.6f  v2=%.6f  diff=%.2e\n",
            p1["beta_urban"], p2["beta_urban"],
            abs(p1["beta_urban"]-p2["beta_urban"])))
cat(sprintf("  beta_normPop: v1=%.6f  v2=%.6f  diff=%.2e\n",
            p1["beta_normPop"], p2["beta_normPop"],
            abs(p1["beta_normPop"]-p2["beta_normPop"])))

cat(sprintf("\n  NLL:  v1=%.6f  v2=%.6f  diff=%.2e\n",
            opt_v1$value, opt_v2$value,
            abs(opt_v1$value - opt_v2$value)))
cat(sprintf("  Time: v1=%.2f s  v2=%.2f s  speedup=%.2fx\n",
            time_v1, time_v2, time_v1/time_v2))

cat("\nDone.\n")
