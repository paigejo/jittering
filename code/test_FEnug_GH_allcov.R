#!/usr/bin/env Rscript
# Compare FE+nugget GH (Q=10, KMICS=100) on real edMICS:
#   (A) covariates = c("urban", "normPop")
#   (B) all covariates (urban, access, elev, distRiversLakes, normPop)

library(TMB)
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
source("code/modM_FEnug.R")

# ── Build shared inputs once (KMICS=100) ─────────────────────────
cat("Building integration point inputs (KMICS=100)...\n"); flush.console()
inputsMDM = makeInputsMDM(datDHS=ed, datMICS=edMICS,
                          KMICS=100, admMICS=admFinal, adm2DHS=adm2Full)
cat("Available covariates:", colnames(inputsMDM$intPtsMICS$XUrb), "\n")

# ═══════ (A) 2 covariates: urban + normPop ═══════════════════════
cat("\n============================================================\n")
cat("  (A) FE+nugget GH Q=10 — covariates: urban, normPop\n")
cat("============================================================\n\n"); flush.console()

fitA = fitMFEM(datDHS=ed, datMICS=edMICS, inputsMDM=inputsMDM,
               covariates=c("urban", "normPop"),
               KMICS=100, Qgh=10,
               fixedEffectsOnly=TRUE, getSDs=TRUE, verbose=FALSE)

cat(sprintf("Convergence: %d | NLL: %.4f | Time: %.1fs + %.1fs\n",
            fitA$opt$convergence, fitA$opt$objective, fitA$totalTime, fitA$sdTime))
if(!is.null(fitA$TMBsd)) {
  estA = summary(fitA$TMBsd, "fixed")
  cat("\nFixed parameter estimates (A):\n"); print(estA)
}

# ═══════ (B) All covariates ══════════════════════════════════════
cat("\n============================================================\n")
cat("  (B) FE+nugget GH Q=10 — ALL covariates\n")
cat("============================================================\n\n"); flush.console()

fitB = fitMFEM(datDHS=ed, datMICS=edMICS, inputsMDM=inputsMDM,
               covariates=NULL,   # NULL = use all covariates
               KMICS=100, Qgh=10,
               fixedEffectsOnly=TRUE, getSDs=TRUE, verbose=FALSE)

cat(sprintf("Convergence: %d | NLL: %.4f | Time: %.1fs + %.1fs\n",
            fitB$opt$convergence, fitB$opt$objective, fitB$totalTime, fitB$sdTime))
if(!is.null(fitB$TMBsd)) {
  estB = summary(fitB$TMBsd, "fixed")
  cat("\nFixed parameter estimates (B):\n"); print(estB)
}

# ═══════ Comparison table ════════════════════════════════════════
cat("\n\n============================================================\n")
cat("  COMPARISON: 2-covariate vs all-covariate (real edMICS)\n")
cat("============================================================\n\n")

if(!is.null(fitA$TMBsd) && !is.null(fitB$TMBsd)) {
  peA = summary(fitA$TMBsd, "fixed"); peB = summary(fitB$TMBsd, "fixed")

  # sigmaEps from log_tauEps
  sigA = exp(-0.5 * peA["log_tauEps", "Estimate"])
  sigB = exp(-0.5 * peB["log_tauEps", "Estimate"])

  cat(sprintf("%-20s %12s %12s\n", "Parameter", "2-cov (A)", "All-cov (B)"))
  cat(sprintf("%-20s %12s %12s\n", "---", "---", "---"))

  # alpha
  cat(sprintf("%-20s %8.4f (%5.4f) %8.4f (%5.4f)\n", "alpha",
              peA["alpha","Estimate"], peA["alpha","Std. Error"],
              peB["alpha","Estimate"], peB["alpha","Std. Error"]))

  # beta — Model A has 2, Model B has all
  betaIdxA = which(rownames(peA) == "beta")
  betaIdxB = which(rownames(peB) == "beta")

  covNamesA = fitA$covNames
  covNamesB = fitB$covNames

  for(j in seq_along(betaIdxA)) {
    cat(sprintf("%-20s %8.4f (%5.4f)",
                paste0("beta[", covNamesA[j], "]"),
                peA[betaIdxA[j], "Estimate"], peA[betaIdxA[j], "Std. Error"]))
    # find matching covariate in B
    matchJ = match(covNamesA[j], covNamesB)
    if(!is.na(matchJ)) {
      cat(sprintf(" %8.4f (%5.4f)", peB[betaIdxB[matchJ], "Estimate"], peB[betaIdxB[matchJ], "Std. Error"]))
    } else {
      cat(sprintf(" %12s", "—"))
    }
    cat("\n")
  }
  # extra covariates only in B
  extraJ = setdiff(seq_along(covNamesB), match(covNamesA, covNamesB))
  for(j in extraJ) {
    cat(sprintf("%-20s %12s %8.4f (%5.4f)\n",
                paste0("beta[", covNamesB[j], "]"), "—",
                peB[betaIdxB[j], "Estimate"], peB[betaIdxB[j], "Std. Error"]))
  }

  cat(sprintf("%-20s %12.4f %12.4f\n", "sigmaEps", sigA, sigB))
  cat(sprintf("%-20s %12.4f %12.4f\n", "NLL", fitA$opt$objective, fitB$opt$objective))
  cat(sprintf("%-20s %12.1f %12.1f\n", "Time (s)", fitA$totalTime, fitB$totalTime))
}

cat("\nDone.\n")
