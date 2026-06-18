# Phi diagnostic on M_M (BYM2 generative, sim 1): is the low phi-hat coming
# from the PC prior, or does the (Laplace) likelihood itself prefer low phi?
#
#   Fit A: default PC prior on phi  (fitMM defaults: u=0.5, alpha=1/3)
#   Fit B: uniform prior on phi     (uniform_phi_prior=TRUE)
#
# For each: phi-hat, sdreport SD on logit_phi (+ delta-method SD on phi), and
# an NLL slice over a dense phi grid (other params at that fit's MLE).
# NOTE obj$fn includes the prior, so Fit B's slice is likelihood + flat-prior
# Jacobian only — that's the "pure likelihood" view on the phi scale.
#
# Interpretation:
#   - phiB ~ 0.8 or B-slice flat from 0.3-0.9  -> prior was the distortion
#   - phiB still ~ 0.3 with sharp B-slice      -> likelihood genuinely prefers
#     low phi (spatial confounding with covariates, or a deeper model bug)
#   - 2*(NLL_B(0.8) - NLL_B(phiB)) > 3.84      -> likelihood "rejects" truth
#     at ~5% (1 dof) — quantifies how strongly

source("code/setup.R")
options(warn = 1)

cat(sprintf("\n=== %s | phi diagnostic: M_M, BYM2 sim 1 ===\n", format(Sys.time())))

simEnv <- new.env()
simulateSurveys("bym2", nsim = 1, regenerate = FALSE, envir = simEnv)
micsEnv <- new.env()
load("savedOutput/global/intPtsMICS_100.RData", envir = micsEnv)

ip <- buildInputsForSim(1, simEnv$surveysDHS, simEnv$surveysMICS,
                        micsEnv$intPtsMICS, "bym2",
                        KMICS = 100, KDHSu = 16, KDHSr = 21)

phiGrid <- c(0.05, 0.10, 0.20, 0.26, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95)

report <- function(res, label) {
    lp   <- res$opt$par["logit_phi"]
    phiH <- plogis(lp)
    sdLP <- tryCatch(sqrt(diag(res$TMBsd$cov.fixed))["logit_phi"],
                     error = function(e) NA)
    sdPhi <- sdLP * phiH * (1 - phiH)        # delta method
    cat(sprintf("\n[%s] phi-hat = %.4f   sd(logit_phi) = %.4f   sd(phi) ~= %.4f\n",
                label, phiH, sdLP, sdPhi))
    cat(sprintf("[%s] NLL slice over phi grid (other params at this fit's MLE):\n", label))
    prof <- nllProfileInPhi(res, phiVals = phiGrid)
    print(prof, row.names = FALSE)
    # likelihood-ratio-style evidence against truth (treating slice as profile)
    nll08  <- prof$NLL[prof$phi == 0.80]
    nllMin <- min(prof$NLL)
    cat(sprintf("[%s] 2*(NLL(0.8) - min NLL) = %.2f   (>3.84 = 'rejects' 0.8 at 5%%, 1 dof)\n",
                label, 2 * (nll08 - nllMin)))
    invisible(prof)
}

cat(sprintf("\n=== %s | Fit A: PC prior on phi ===\n", format(Sys.time())))
tA <- proc.time()[3]
fitA <- fitMM(datDHS = ip$datDHS, datMICS = ip$datMICS, inputsMDM = ip$inputsMDM,
              KMICS = 100, Qgh = 10, getSDs = TRUE, verbose = FALSE)
cat(sprintf("fit A time: %.1f min\n", (proc.time()[3] - tA)/60))
profA <- report(fitA, "A: PC prior")

cat(sprintf("\n=== %s | Fit B: uniform prior on phi ===\n", format(Sys.time())))
tB <- proc.time()[3]
fitB <- fitMM(datDHS = ip$datDHS, datMICS = ip$datMICS, inputsMDM = ip$inputsMDM,
              KMICS = 100, Qgh = 10, getSDs = TRUE, verbose = FALSE,
              uniform_phi_prior = TRUE)
cat(sprintf("fit B time: %.1f min\n", (proc.time()[3] - tB)/60))
profB <- report(fitB, "B: uniform prior")

# Quantify the PC prior's pull: its log-density along the grid
cat("\n=== PC prior log-density along the grid (u=0.5, alpha=1/3) ===\n")
out <- load("savedOutput/global/admFinalMat.RData")
bArgs <- prepareBYM2argumentsForTMB(admFinalMat, u = 0.5, alpha = 1/3,
                                    constr = TRUE, scale.model = TRUE,
                                    matrixType = "TsparseMatrix")
priLD <- sapply(phiGrid, function(p)
    dBYM2phiPC(p, lambda = bArgs$lambda,
               gammaTildesm1 = bArgs$gammaTildesm1, tr = bArgs$tr, doLog = TRUE))
print(data.frame(phi = phiGrid, priorLogDens = round(priLD, 3),
                 priorPullVsPhi08 = round(priLD - priLD[phiGrid == 0.80], 3)),
      row.names = FALSE)

cat("\n=== Summary ===\n")
cat(sprintf("phi truth        : 0.80\n"))
cat(sprintf("phi-hat PC prior : %.4f\n", plogis(fitA$opt$par["logit_phi"])))
cat(sprintf("phi-hat uniform  : %.4f\n", plogis(fitB$opt$par["logit_phi"])))
cat(sprintf("\n=== %s | done ===\n", format(Sys.time())))
