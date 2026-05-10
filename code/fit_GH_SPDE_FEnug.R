#!/usr/bin/env Rscript
# FE + nugget(GH) with Q=10 on SPDE-simulated population (sim 1)
# Uses optimized v2 template

library(TMB)
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
source("code/modM_MIIDonly.R")

# ── Load SPDE simulated data (not BYM2) ──────────────────────────
cat("Loading SPDE simulation...\n")
load("savedOutput/simStudy1/simPopsSurveys.RData")
ed = surveysDHS[[1]]
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

# Subset to urban + normPop (same as BYM2 FE+nug runs)
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
cat(sprintf("Using %d Gauss-Hermite quadrature nodes\n", Q))
cat(sprintf("nAreas: %d, nUrbObs: %d, nRurObs: %d\n",
            nAreas, length(ysUrbMICS), length(ysRurMICS)))

# ── True values (same DGP as BYM2 sim) ──────────────────────────
trueAlpha = -1.25; trueUrban = 1.00; trueNormPop = 0.50
trueSigmaEps = sqrt(1.5)

# ── Initial values ───────────────────────────────────────────────
initUrbP = sum(ysUrbMICS)/sum(nsUrbMICS)
initRurP = sum(ysRurMICS)/sum(nsRurMICS)
initAlpha = log(initRurP/(1-initRurP))
initBeta1 = log(initUrbP/(1-initUrbP)) - initAlpha

# ── TMB Data (v2: includes lchoose) ─────────────────────────────
data_gh = list(
  y_iUrbanMICS=ysUrbMICS, y_iRuralMICS=ysRurMICS,
  n_iUrbanMICS=nsUrbMICS, n_iRuralMICS=nsRurMICS,
  areaidxlocUrbanMICS=areaidxlocUrbanMICS,
  areaidxlocRuralMICS=areaidxlocRuralMICS,
  X_betaUrbanMICS=intPtsMICS$XUrb,
  X_betaRuralMICS=intPtsMICS$XRur,
  wUrbanMICS=intPtsMICS$wUrban,
  wRuralMICS=intPtsMICS$wRural,
  lchoose_urban=lchoose(nsUrbMICS, ysUrbMICS),
  lchoose_rural=lchoose(nsRurMICS, ysRurMICS),
  gh_nodes=gh$x, gh_weights=gh$w,
  alpha_pri=c(0, 100^2), beta_pri=c(0, sqrt(1000)),
  lambdaTau=lambdaTau, lambdaTauEps=lambdaTauEps
)

# Clean up
rm(surveysDHS, surveysMICS, inputsMDM); gc()

# ── Compile v2 ───────────────────────────────────────────────────
unloadDynlibs()
cat("\nCompiling modM_MIIDonly_GH_v2.cpp...\n")
compile("code/modM_MIIDonly_GH_v2.cpp", framework="TMBad", safebounds=FALSE)
dyn.load(dynlib("code/modM_MIIDonly_GH_v2"))

# ── Config A: FE+nug, Q=10, BFGS ────────────────────────────────
cat("\n============================================================\n")
cat("SPDE sim 1: FE+NUGGET(GH Q=10 v2), BFGS\n")
cat("============================================================\n")

params_A = list(
  log_tau=0, log_tauEps=0,
  alpha=initAlpha, beta_urban=initBeta1, beta_normPop=0,
  u_spatial=rep(0, nAreas)
)
mapA = list(u_spatial=factor(rep(NA, nAreas)), log_tau=factor(NA))

obj_A = MakeADFun(data=data_gh, parameters=params_A,
                  map=mapA, DLL='modM_MIIDonly_GH_v2', silent=TRUE)

cat("Parameters to optimize:", names(obj_A$par), "\n")
cat("Number of params:", length(obj_A$par), "\n\n"); flush.console()

t0 = proc.time()[3]
opt_A = optim(obj_A$par, obj_A$fn, obj_A$gr, method="BFGS",
              control=list(trace=1, REPORT=1, maxit=500))
time_A = proc.time()[3] - t0

cat(sprintf("\nConvergence: %d\n", opt_A$convergence))
cat(sprintf("NLL: %.4f\n", opt_A$value))
cat(sprintf("Time: %.1f s\n", time_A))

sd_A = sdreport(obj_A)
est_A = summary(sd_A, "fixed")
cat("\nFixed parameter estimates:\n"); print(est_A)

pe = est_A[,"Estimate"]; se = est_A[,"Std. Error"]
nms = rownames(est_A)
sigEps_est = exp(-0.5*pe[nms=="log_tauEps"])

cat(sprintf("\n--- SPDE sim1: FE+nugget(GH Q=10 v2) ---\n"))
cat(sprintf("  alpha:       %7.4f (SE %.4f)  truth: %7.4f\n",
            pe[nms=="alpha"], se[nms=="alpha"], trueAlpha))
cat(sprintf("  beta_urban:  %7.4f (SE %.4f)  truth: %7.4f\n",
            pe[nms=="beta_urban"], se[nms=="beta_urban"], trueUrban))
cat(sprintf("  beta_normPop:%7.4f (SE %.4f)  truth: %7.4f\n",
            pe[nms=="beta_normPop"], se[nms=="beta_normPop"], trueNormPop))
cat(sprintf("  sigmaEps:    %7.4f             truth: %7.4f\n",
            sigEps_est, trueSigmaEps))
cat(sprintf("  Time: %.1f s\n", time_A))

# ── Save ─────────────────────────────────────────────────────────
result_spde_FEnug = list(opt=opt_A, sd=sd_A, time=time_A, Q=Q)
save(result_spde_FEnug, file="savedOutput/fitGH_SPDE_FEnug_Q10.RData")
cat("\nSaved to savedOutput/fitGH_SPDE_FEnug_Q10.RData\n")
cat("Done.\n")
