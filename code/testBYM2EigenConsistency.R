# Does simPopBYM2 produce IDENTICAL BYM2 effects when given the bym2ArgsTMB
# struct vs when it computes Q + eigen.spam internally? With the
# pass-through wired up, both paths should produce bit-identical output
# given the same seed.
#
# Also verifies: with the new sum-to-zero zeroing of the constant-mode
# eigenvariance, sum(Epsilon_bym2) over all areas is exactly 0 (machine
# precision), matching the TMB template's sum(w_bym2)=0 constraint.

source("code/setup.R")

cat(sprintf("\n=== %s | BYM2 eigen-consistency test ===\n", format(Sys.time())))

# Build bym2ArgsTMB once (same call as modM_DMSep.R uses)
out <- load("savedOutput/global/admFinalMat.RData")
bym2ArgsTMB <- prepareBYM2argumentsForTMB(admFinalMat, u = 0.5, alpha = 2/3,
                                          constr = TRUE, scale.model = TRUE,
                                          matrixType = "TsparseMatrix")
cat(sprintf("bym2ArgsTMB: nrow(Q) = %d, length(gammaTildesm1) = %d\n",
            nrow(bym2ArgsTMB$Q), length(bym2ArgsTMB$gammaTildesm1)))

# Minimal inputs for simPopBYM2: we only care about the BYM2 effect sampling,
# not the survey machinery. Pass nsim=1 with a fixed seed each call.
# Make a small popMat tied to admFinalMat's areas.
nAreas    <- nrow(admFinalMat)
mockPopMat <- data.frame(
    east     = seq(1, nAreas) * 10,
    north    = seq(1, nAreas) * 10,
    urban    = rep(c(0, 1), length.out = nAreas),
    stratumMICS = seq_len(nAreas)
)

# Pull the BYM2-effect-sampling block out of simPopBYM2 into a tiny helper
# so we can compare deterministically across the two paths. Mirrors lines
# 700-735 of simData.R (the part we care about).
sampleBYM2 <- function(bym2ArgsTMB_arg = NULL, seed = 42L,
                        sigmaBYM2 = sqrt(0.5), phi = 0.8) {
    set.seed(seed)
    if (is.null(bym2ArgsTMB_arg)) {
        Q             <- makeQBesag(admFinalMat, constr = TRUE,
                                    scale.model = TRUE, matrixType = "spam")
        eigQ          <- eigen.spam(Q, symmetric = TRUE)
        gammas        <- eigQ$values
        V             <- eigQ$vectors
        tol           <- 1e-8
        gammaTildes   <- 1 / gammas
        gammaTildes[abs(gammas) < tol] <- 0
        gammaTildesm1 <- gammaTildes - 1
        path          <- "built-from-scratch"
    } else {
        Q             <- bym2ArgsTMB_arg$Q
        V             <- bym2ArgsTMB_arg$V
        gammaTildesm1 <- bym2ArgsTMB_arg$gammaTildesm1
        gammaTildes   <- gammaTildesm1 + 1
        tol           <- 1e-8
        gammas        <- ifelse(abs(gammaTildes) > tol, 1 / gammaTildes, 0)
        path          <- "from-bym2ArgsTMB"
    }
    nGraphAreas <- nrow(Q)
    tau         <- 1 / sigmaBYM2^2
    eigenVars   <- (1/tau) * (1 + phi * gammaTildesm1)
    eigenVars[eigenVars < 0] <- 0
    zeroModeIdx <- which(abs(gammas) < tol)
    if (length(zeroModeIdx) > 0) eigenVars[zeroModeIdx] <- 0
    eigenSDs    <- sqrt(eigenVars)
    z           <- rnorm(nGraphAreas)
    Eps         <- as.numeric(V %*% (eigenSDs * z))
    list(path = path, Eps = Eps, sumEps = sum(Eps),
         maxAbsEps = max(abs(Eps)), nZeroMode = length(zeroModeIdx))
}

A <- sampleBYM2(bym2ArgsTMB_arg = NULL,        seed = 42L)
B <- sampleBYM2(bym2ArgsTMB_arg = bym2ArgsTMB, seed = 42L)

cat(sprintf("\n[A] %s : sum(Eps) = %.3e   max|Eps| = %.4f   nZeroModes = %d\n",
            A$path, A$sumEps, A$maxAbsEps, A$nZeroMode))
cat(sprintf("[B] %s : sum(Eps) = %.3e   max|Eps| = %.4f   nZeroModes = %d\n",
            B$path, B$sumEps, B$maxAbsEps, B$nZeroMode))

cat(sprintf("\nmax |A$Eps - B$Eps| = %.3e\n", max(abs(A$Eps - B$Eps))))
cat(sprintf("identical(A$Eps, B$Eps) : %s\n", identical(A$Eps, B$Eps)))

# Sum-to-zero check: with the new zeroing of the constant mode, both paths
# should give sum(Eps) at machine epsilon (~1e-14 to 1e-15 due to floating
# rounding from the V*z multiplication, NOT 1e-2 like before the fix).
threshold <- 1e-9
cat(sprintf("\nSum-to-zero check (threshold %g):\n", threshold))
cat(sprintf("  A: |sum(Eps)| = %.3e  %s\n", abs(A$sumEps),
            ifelse(abs(A$sumEps) < threshold, "PASS", "FAIL")))
cat(sprintf("  B: |sum(Eps)| = %.3e  %s\n", abs(B$sumEps),
            ifelse(abs(B$sumEps) < threshold, "PASS", "FAIL")))

# Reproducibility check: same seed, same eigenvars, same Eps?
if (identical(A$Eps, B$Eps)) {
    cat("\nPASS: paths produce bit-identical BYM2 effects.\n")
} else if (max(abs(A$Eps - B$Eps)) < 1e-12) {
    cat("\nPASS (within 1e-12): paths produce numerically-identical BYM2 effects.\n")
} else {
    cat("\nFAIL: paths differ. This indicates Q / V / gammaTildesm1 are not the same.\n")
    cat(sprintf("Max diff = %.3e\n", max(abs(A$Eps - B$Eps))))
}

cat(sprintf("\n=== %s | done ===\n", format(Sys.time())))
