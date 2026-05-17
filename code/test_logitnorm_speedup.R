# Accuracy / speed comparison for E[expit(mu + sigma Z)] approximations,
# Z ~ N(0,1).  Three approximators are compared against numerical integration:
#
#   - GH-Q       : Gauss-Hermite quadrature, fully vectorized over pixels & draws
#   - bilin2D    : pretrained 2D lookup table + bilinear interpolation
#   - spline1D   : the current logitNormMeanSplineApprox path (per-draw monotone
#                  spline; what scoreSpatial was using before the patch).
#
# Reference: numerical integration via integrate(), per (mu, sigma) pair.
#
# Test 1: small grid of (mu, sigma) values — accuracy.
# Test 2: realistic-scale workload — speed at nGrid * nDraw shape.

source("setup.R")
options(error=traceback)

# ============================================================================
# 1) Reference implementation: brute-force numerical integration
# ============================================================================
refMean <- function(mu, sigma) {
    if(sigma == 0) return(plogis(mu))
    integrate(function(z) plogis(mu + sigma*z) * dnorm(z),
              lower = -10, upper = 10, abs.tol = 1e-12)$value
}

# ============================================================================
# 2) Gauss-Hermite quadrature (vectorized)
# ============================================================================
# E[expit(mu + sigma Z)] = (1/sqrt(pi)) * sum_q w_q * expit(mu + sqrt(2) sigma u_q)
makeGHfn <- function(Q = 30) {
    gh <- fastGHQuad::gaussHermiteData(Q)
    nodes <- gh$x; weights <- gh$w
    # Vectorized: mus and sigmas are matrices of the same shape (nGrid x nDraw).
    # We assume sigma is constant within each column of `mus` — i.e., sigmas is
    # a length-ncol(mus) vector (one sigma per draw). This matches the predGrid
    # use case and lets us avoid per-pixel sigma broadcasting.
    function(mus, sigmasPerCol) {
        out <- matrix(0, nrow = nrow(mus), ncol = ncol(mus))
        sqrt2sigma <- sqrt(2) * sigmasPerCol           # length ncol(mus)
        for(q in seq_along(nodes)) {
            shifts <- nodes[q] * sqrt2sigma
            out <- out + weights[q] * plogis(sweep(mus, 2, shifts, "+"))
        }
        out / sqrt(pi)
    }
}

# ============================================================================
# 3) Pretrained 2D bilinear lookup
# ============================================================================
# Build lookup table f(mu, sigma) = E[expit(mu + sigma Z)] on a regular grid.
# Symmetry isn't exploited to keep the code minimal; saturation handled by snap-
# ping out-of-range inputs to the grid boundary (where expit is essentially 0/1).
buildBilinGrid <- function(muRange = c(-15, 15), nMu = 200,
                           sigmaRange = c(0.01, 5), nSigma = 100) {
    muG    <- seq(muRange[1],    muRange[2],    length.out = nMu)
    sigmaG <- seq(sigmaRange[1], sigmaRange[2], length.out = nSigma)
    G      <- matrix(0, nrow = nMu, ncol = nSigma)
    for(j in seq_along(sigmaG)) {
        s <- sigmaG[j]
        for(i in seq_along(muG)) G[i, j] <- refMean(muG[i], s)
    }
    list(muGrid = muG, sigmaGrid = sigmaG, G = G)
}

# Vectorized bilinear interpolation. `mus`, `sigmas` same shape (typically
# matrices nGrid x nDraw) or broadcastable. Returns shape of `mus`.
bilinLookup <- function(grid, mus, sigmas) {
    muG <- grid$muGrid;  sG <- grid$sigmaGrid
    nMu <- length(muG);  nSig <- length(sG)
    Gmat <- grid$G

    iM <- pmin(pmax(findInterval(mus,    muG), 1L), nMu  - 1L)
    iS <- pmin(pmax(findInterval(sigmas, sG ), 1L), nSig - 1L)

    mu0 <- muG[iM];  mu1 <- muG[iM + 1L]
    s0  <- sG [iS];  s1  <- sG [iS + 1L]
    tM  <- (mus    - mu0) / (mu1 - mu0)
    tS  <- (sigmas - s0 ) / (s1  - s0 )

    flat00 <- iM        + (iS - 1L) * nMu
    flat10 <- (iM + 1L) + (iS - 1L) * nMu
    flat01 <- iM        +  iS       * nMu
    flat11 <- (iM + 1L) +  iS       * nMu

    res <- Gmat[flat00] * (1 - tM) * (1 - tS) +
           Gmat[flat10] *  tM      * (1 - tS) +
           Gmat[flat01] * (1 - tM) *  tS      +
           Gmat[flat11] *  tM      *  tS
    if(is.matrix(mus)) dim(res) <- dim(mus)
    res
}

# ============================================================================
# 4) Existing per-draw spline (what scoreSpatial used previously)
# ============================================================================
spline1DPerDraw <- function(mus, sigmasPerCol) {
    nDraw <- ncol(mus)
    out <- matrix(0, nrow = nrow(mus), ncol = nDraw)
    for(s in seq_len(nDraw)) {
        out[, s] <- logitNormMeanSplineApprox(mus[, s], sigmasPerCol[s])$vals
    }
    out
}

# ============================================================================
# Test 1: ACCURACY on a small grid of (mu, sigma) values
# ============================================================================
cat("\n========== Test 1: accuracy ==========\n")
testMu    <- seq(-8, 8, by = 0.5)
testSigma <- c(0.1, 0.5, 1.0, 1.225, 2.0, 3.0)
testPts   <- expand.grid(mu = testMu, sigma = testSigma)
testPts$ref <- mapply(refMean, testPts$mu, testPts$sigma)

# GH
ghFn <- makeGHfn(Q = 30)
testPts$gh <- mapply(function(mu, sg) ghFn(matrix(mu), sg), testPts$mu, testPts$sigma)

# 2D bilinear (build once)
cat("Building 2D lookup grid (200 x 100) ... ")
t0 <- proc.time()[3]
bgrid <- buildBilinGrid()
cat(sprintf("built in %.1f sec\n", proc.time()[3] - t0))
testPts$bilin <- bilinLookup(bgrid, testPts$mu, testPts$sigma)

cat("Max abs error vs reference:\n")
cat(sprintf("  GH Q=30          : %.3e\n", max(abs(testPts$gh    - testPts$ref))))
cat(sprintf("  bilin 2D lookup  : %.3e\n", max(abs(testPts$bilin - testPts$ref))))

# Also tabulate worst (mu, sigma) for each approximator
worst <- function(col) {
    i <- which.max(abs(testPts[[col]] - testPts$ref))
    sprintf("worst at mu=%.2f sigma=%.2f, est=%.6f ref=%.6f (err=%.2e)",
            testPts$mu[i], testPts$sigma[i], testPts[[col]][i],
            testPts$ref[i], testPts[[col]][i] - testPts$ref[i])
}
cat("  GH:    ", worst("gh"),    "\n")
cat("  bilin: ", worst("bilin"), "\n")

# ============================================================================
# Test 2: SPEED at realistic predGrid shapes
# ============================================================================
cat("\n========== Test 2: speed ==========\n")
nGrid <- 36564L      # nrow(popMatNGAThresh)
nDraw <- 5000L
set.seed(42)
bigMu        <- matrix(rnorm(nGrid * nDraw, mean = 0, sd = 2), nrow = nGrid)
sigmaPerCol  <- runif(nDraw, 0.7, 1.6)

# (a) 2D bilinear, broken into build vs run -----------------------------------
cat("2D bilinear  build (200 mu x 100 sigma = 20000 integrations) ... ")
t0 <- proc.time()[3]
bgrid <- buildBilinGrid()
buildSec <- proc.time()[3] - t0
cat(sprintf("%.2f sec\n", buildSec))

cat("2D bilinear   run (one vectorized lookup, ", nGrid, "x", nDraw, ") ... ", sep = "")
sigmaMat <- matrix(sigmaPerCol, nrow = nGrid, ncol = nDraw, byrow = TRUE)
t0 <- proc.time()[3]
res_bilin <- bilinLookup(bgrid, bigMu, sigmaMat)
cat(sprintf("%.2f sec  (build amortizes if reused)\n", proc.time()[3] - t0))

# (b) GH Q=30 -----------------------------------------------------------------
cat("GH Q=30 vectorized over Q nodes,      ", nGrid, "x", nDraw, " ... ", sep = "")
t0 <- proc.time()[3]
res_gh <- ghFn(bigMu, sigmaPerCol)
cat(sprintf("%.2f sec\n", proc.time()[3] - t0))

# (c) Per-draw 1D spline (hand-rolled for-loop) at several nDraw values to
#     verify that the cost is genuinely linear-in-nDraw and to give the user
#     a real (not extrapolated) number at nDraw=1000.
cat("\nPer-draw 1D spline (hand-rolled for loop), linearity check:\n")
for(nD in c(100L, 500L, 1000L)) {
    t0 <- proc.time()[3]
    res_sp1d <- spline1DPerDraw(bigMu[, 1:nD, drop = FALSE], sigmaPerCol[1:nD])
    el <- proc.time()[3] - t0
    cat(sprintf("  nDraw=%4d : %6.2f sec   (per-call %.3f sec)\n",
                nD, el, el/nD))
}

# (d) logitNormMeanGrouped — the function predGrid's main branch actually uses.
#     Takes a (1+nGrid) x nDraw matrix (sigma row stacked on top of mus).
cat("\nlogitNormMeanGrouped(splineApprox=TRUE) — actual main-branch fn:\n")
for(nD in c(100L, 500L, 1000L)) {
    mat_in <- rbind(sigmaPerCol[1:nD], bigMu[, 1:nD, drop = FALSE])
    t0 <- proc.time()[3]
    res_grp <- logitNormMeanGrouped(mat_in, logisticApprox = FALSE, splineApprox = TRUE)
    el <- proc.time()[3] - t0
    cat(sprintf("  nDraw=%4d : %6.2f sec   (per-call %.3f sec)\n",
                nD, el, el/nD))
}

# Sanity check across methods at full scale
maxDiff <- max(abs(res_gh - res_bilin))
cat(sprintf("\nMax abs diff bilin vs GH at full scale: %.3e\n", maxDiff))

cat("\nDone.\n")
