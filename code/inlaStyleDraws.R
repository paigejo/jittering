# ============================================================================
# INLA-style smooth hyperparameter integration for BYM2 GH templates.
#
# Implements the same scheme used by INLA's inla.hyperpar.sample():
#   1. Mode + Hessian + eigen-rotation of cov.fixed restricted to hypers.
#   2. Multi-step axial walk along each eigen-axis with step δ_z (default 1.0)
#      until ΔNLL exceeds δ_π (default 2.5), capped at maxAxialSteps.
#   3. At each axial point, hold the three hypers fixed (via `map`) and re-fit
#      the inner block (alpha, beta, w_bym2Free, u_bym2Free) under Laplace.
#      Cache the inner mode μ_k and inner precision Q_k (TMB's spHess).
#   4. Per-axis skew-normal scales (σⱼ⁺, σⱼ⁻) by least-squares fit of
#      ΔNLL_i = z_i² / (2σ²) across the axial points on each side.
#   5. Sample hypers via product of split-Gaussians:
#        direction ~ Bernoulli(σ⁺ / (σ⁺+σ⁻)); z = ±|N(0, σ_direction²)|.
#   6. For each draw, pick the closest CCD point (in z-space) and sample
#      inner ~ N(μ_k, Q_k⁻¹) from the cached Cholesky of Q_k.
#
# Matches INLA's R code for steps 1, 4, 5; differs from INLA's C engine in
# CCD point set (no Cartesian product) and improve.marginals post-step
# (skipped here).
# ============================================================================

# -- internal: draw n values from a two-piece (split) Gaussian with σ⁺=sp, σ⁻=sn
.rSplitGaussian <- function(n, sp, sn) {
    p_pos  <- sp / (sp + sn)
    sign_v <- ifelse(runif(n) < p_pos, +1, -1)
    sd_v   <- ifelse(sign_v > 0, sp, sn)
    sign_v * abs(rnorm(n, sd = sd_v))
}

# -- internal: evaluate (nll, μ, Q) at one CCD point REUSING a single shared
# TMB object instead of rebuilding the AD tape. walkObj has the hyperparameters
# as free OUTER parameters and (alpha, beta, w_bym2Free, u_bym2Free) as random,
# so walkObj$fn(theta) returns the Laplace-marginal NLL at hyper = theta after
# internally optimizing the inner random effects — equivalent to rebuilding a
# fresh fixed-hyper MakeADFun per point, but without re-taping (~24 s/point).
# Across a ~19-point walk this turns ~19 tape builds into one (validated
# bit-equivalent to the per-point-rebuild draws).
#
# Two subtleties this function handles, both confirmed by testing:
#  (1) read last.par (the MOST RECENT evaluation = [theta, inner_mode(theta)]),
#      NOT last.par.best (which tracks the lowest-NLL theta across all calls).
#  (2) DEEP-COPY Q after spHess: spHess reuses one internal value buffer per
#      object, so on a shared object every returned Q aliases and the next call
#      overwrites it — without the copy all points collapse onto the last Q.
.evalCCDreuse <- function(walkObj, theta_k, innerStart = NULL) {
    randIdx <- walkObj$env$random
    # Match the rebuild path's warm-start: every CCD point optimizes the
    # random effects starting from the SAME inner mode (the centre mode),
    # not sequentially from the previous point. The inner problem has a flat
    # (unidentified) spatial ridge, so the optimizer's stopping point depends
    # on its start; equalizing the start makes reuse reproduce the rebuild
    # modes (and hence draws) to optimizer tolerance. TMB warm-starts the
    # inner solve from last.par[random], so we overwrite that block here.
    if(!is.null(innerStart)) {
        lp <- walkObj$env$last.par
        lp[randIdx] <- innerStart
        walkObj$env$last.par <- lp
        if(!is.null(walkObj$env$last.par.best)) {
            lpb <- walkObj$env$last.par.best
            lpb[randIdx] <- innerStart
            walkObj$env$last.par.best <- lpb
        }
    }
    theta <- as.numeric(theta_k[names(walkObj$par)])  # order to walkObj$par
    nll   <- as.numeric(walkObj$fn(theta))
    fullPar <- walkObj$env$last.par                   # [theta, inner_mode(theta)]
    mu      <- fullPar[randIdx]                        # named inner mode
    Q       <- tryCatch(walkObj$env$spHess(fullPar, random = TRUE),
                        error = function(e) NULL)
    # CRITICAL: spHess returns a matrix that reuses ONE internal value buffer
    # per object. Because we reuse a single walkObj across all CCD points, the
    # next spHess call OVERWRITES this matrix's values in place — so every Q we
    # store would alias and collapse onto the LAST point's Hessian (confirmed:
    # a stored reference's values changed after the next spHess call, and ended
    # up equal to the next point's Q). Deep-copy here so each point keeps its
    # own precision. `+ 0` materialises an independent sparse copy. The rebuild
    # path is immune (fresh object per point, spHess called once), which is why
    # it was correct. Without this copy, reuse silently produces wrong draws.
    if(!is.null(Q)) Q <- Q + 0
    list(nll = nll, mu = mu, Q = Q)
}

# -- internal: fit per-side σ via least-squares: ΔNLL_i = z_i² / (2σ²)
.fitSigmaLS <- function(zList, dNLLList) {
    # zList, dNLLList are paired numeric vectors (same length)
    keep <- which(dNLLList > 0)
    if(length(keep) == 0) return(1)        # degenerate; pretend symmetric Gaussian
    x <- zList[keep]^2
    y <- dNLLList[keep]
    # Closed-form OLS for y = β x:  β̂ = Σxy / Σx².  Then σ = √(1/(2β̂)).
    beta_hat <- sum(x*y) / sum(x^2)
    if(beta_hat <= 0) return(1)
    1 / sqrt(2 * beta_hat)
}

# ============================================================================
# Single dispatch entry point. Use this everywhere (predGrid, scoring, tests)
# so the two sampling strategies live in exactly one place.
#
#   posteriorDraws(res, NDRAWS, useInla, ...)
#     useInla = "auto"  (default)  ->  try joint-Gaussian first; if Cholesky
#                                      of jointPrecision fails (e.g. pdHess
#                                      = FALSE), automatically fall back to
#                                      inlaStyleDraws().
#     useInla = TRUE   / "yes"     ->  force inlaStyleDraws().
#     useInla = FALSE  / "no"      ->  force the joint-Gaussian draws.
#                                      Errors propagate if the Cholesky fails.
#
# Both branches return a matrix with rownames = parameter names and
# ncol = NDRAWS. Auto-fallback is the recommended default: the Gaussian path
# is ~25× faster and indistinguishable when it works, and the INLA path is
# robust when it doesn't.
# ============================================================================
# Per-name agreement check between draw-implied SDs and sdreport's reported
# SDs. Two independent TMB code paths: if their per-name SDs agree, the
# draw matrix is correctly labelled (and per-name disagreement localizes any
# permutation bug). Returns the diagnostic table invisibly and warns on any
# large mismatch (default tolerance scales as relSE / sqrt(NDRAWS)).
.checkDrawLabels <- function(draws, res, relTol = 8) {
    SD <- res$TMBsd
    if(!inherits(SD, "sdreport")) return(invisible(NULL))
    # Path B: per-name SEs via summary(SD) — keyed by user-declared names.
    summTab <- tryCatch(summary(SD, "report"), error = function(e) NULL)
    if(is.null(summTab) || nrow(summTab) == 0)
        summTab <- tryCatch(summary(SD), error = function(e) NULL)
    if(is.null(summTab) || nrow(summTab) == 0) return(invisible(NULL))
    sePaths <- summTab[, "Std. Error"]
    seNames <- rownames(summTab)
    # Path A: marginal SDs from the draws (keyed by rownames(draws), the
    # names we attached when sampling).
    drawSDs   <- apply(draws, 1, sd)
    drawNames <- rownames(draws)
    # Per-name compare on the intersection.
    common <- intersect(drawNames, seNames)
    if(length(common) == 0) {
        warning(".checkDrawLabels: no shared names between draws and sdreport summary")
        return(invisible(NULL))
    }
    tab <- data.frame(
        name      = common,
        sd_draws  = drawSDs[common],
        sd_report = sePaths[common],
        diff      = drawSDs[common] - sePaths[common],
        relDiff   = (drawSDs[common] - sePaths[common]) /
                    pmax(abs(sePaths[common]), 1e-12),
        stringsAsFactors = FALSE
    )
    # Expected MC SE on sd of N draws is sd / sqrt(2*(N-1)). Flag when
    # observed |relDiff| > relTol * 1/sqrt(2*(N-1)).
    NDRAWS <- ncol(draws)
    relTolEff <- relTol / sqrt(2 * max(NDRAWS - 1, 1))
    bad <- which(abs(tab$relDiff) > relTolEff)
    if(length(bad) > 0) {
        warning(sprintf(
            ".checkDrawLabels: %d/%d parameters disagree beyond %.3f*MC-SE:\n%s",
            length(bad), nrow(tab), relTol,
            paste0("  [", tab$name[bad], "] draws=", round(tab$sd_draws[bad], 5),
                   " report=", round(tab$sd_report[bad], 5),
                   " relDiff=", round(tab$relDiff[bad], 3), collapse = "\n")))
    }
    invisible(tab)
}

posteriorDraws <- function(res, NDRAWS = 1000, useInla = "auto", ...) {
    mode <- if(identical(useInla, "auto"))           "auto"
            else if(isTRUE(useInla)  || identical(useInla, "yes")) "inla"
            else if(isFALSE(useInla) || identical(useInla, "no"))  "gauss"
            else stop("posteriorDraws: useInla must be one of \"auto\", TRUE/\"yes\", or FALSE/\"no\"")

    if(mode == "inla") {
        out <- inlaStyleDraws(res, NDRAWS = NDRAWS, ...)
        .checkDrawLabels(out, res)
        return(out)
    }

    if(mode == "gauss") {
        out <- .gaussianPosteriorDraws(res, NDRAWS = NDRAWS,
                                       requireJoint = TRUE)
        .checkDrawLabels(out, res)
        return(out)
    }

    # auto: try Gaussian (with strict jointPrecision Cholesky); on failure
    # fall back to INLA-style and emit a one-line note so the choice is
    # visible in scoring logs.
    out <- tryCatch(
        .gaussianPosteriorDraws(res, NDRAWS = NDRAWS, requireJoint = TRUE),
        error = function(e) {
            message("[posteriorDraws] Gaussian draws failed (", conditionMessage(e),
                    "); falling back to INLA-style.")
            NULL
        })
    if(is.null(out)) out <- inlaStyleDraws(res, NDRAWS = NDRAWS, ...)
    .checkDrawLabels(out, res)
    out
}

# Joint-Gaussian draws. When `requireJoint = TRUE`, errors out if the
# jointPrecision Cholesky fails — the auto branch of posteriorDraws relies on
# that to detect non-PSD and trigger the INLA fallback. When FALSE, falls
# back to cov.fixed-only (the historical predGrid path; loses random-effect
# uncertainty but never crashes outright).
.gaussianPosteriorDraws <- function(res, NDRAWS = 1000, requireJoint = TRUE) {
    SD <- res$TMBsd
    if(!inherits(SD, "sdreport"))
        stop(".gaussianPosteriorDraws: res$TMBsd must be a valid sdreport object")

    if(!is.null(SD$jointPrecision)) {
        L <- tryCatch(Matrix::Cholesky(SD$jointPrecision, LDL = FALSE, super = NA),
                      error = function(e) e)
        if(!inherits(L, "error")) {
            mode  <- c(SD$par.fixed, SD$par.random)
            # Defensive: TMB's sdreport sometimes reorders parameters internally
            # such that the row/col labels of jointPrecision don't match
            # names(c(par.fixed, par.random)). If we draw and assume the
            # canonical order without checking, we'd get mis-labelled draws
            # (right marginals attached to the wrong parameter). Verify here.
            jpNames <- rownames(SD$jointPrecision)
            if(!is.null(jpNames)) {
                if(length(jpNames) != length(mode) || !all(jpNames == names(mode))) {
                    # Reorder Q to match c(par.fixed, par.random) before drawing.
                    perm  <- match(names(mode), jpNames)
                    if(anyNA(perm))
                        stop(".gaussianPosteriorDraws: names(c(par.fixed, par.random)) ",
                             "not all found in rownames(jointPrecision); cannot reorder ",
                             "safely. Missing: ",
                             paste(names(mode)[is.na(perm)], collapse = ", "))
                    warning(".gaussianPosteriorDraws: jointPrecision label order ",
                            "differs from c(par.fixed, par.random); reordering. ",
                            "First mismatch at idx ",
                            which(jpNames != names(mode))[1])
                    Qr <- SD$jointPrecision[perm, perm]
                    L  <- tryCatch(Matrix::Cholesky(Qr, LDL = FALSE, super = NA),
                                   error = function(e) e)
                    if(inherits(L, "error"))
                        stop("Cholesky of reordered jointPrecision failed: ",
                             conditionMessage(L))
                }
            }
            z     <- matrix(rnorm(length(mode) * NDRAWS), nrow = length(mode))
            shift <- Matrix::solve(L, z, system = "Lt")
            shift <- Matrix::solve(L, shift, system = "Pt")
            out   <- as.matrix(shift) + mode
            rownames(out) <- names(mode)
            return(out)
        }
        if(requireJoint)
            stop("jointPrecision Cholesky failed: ", conditionMessage(L))
    } else if(requireJoint) {
        stop("res$TMBsd has no jointPrecision (set getJointPrecision=TRUE in sdreport)")
    }
    # cov.fixed-only fallback (predGrid's historical path)
    mF <- SD$par.fixed
    Lc <- chol(SD$cov.fixed)
    z  <- matrix(rnorm(length(mF) * NDRAWS), nrow = length(mF))
    fe <- mF + t(Lc) %*% z
    rownames(fe) <- names(mF)
    fe
}

# ============================================================================
# INLA-style hyperparameter integration (called via posteriorDraws(useInla=T))
# ============================================================================
inlaStyleDraws <- function(res, NDRAWS = 1000,
                            deltaZ = 1.0, deltaPi = 2.5, maxAxialSteps = 4,
                            hyperNames = NULL,
                            verbose = FALSE) {
    SD <- res$TMBsd
    if(!inherits(SD, "sdreport"))
        stop("res$TMBsd must be a valid sdreport object")
    obj            <- res$TMBobj
    dataList       <- obj$env$data
    paramsTemplate <- obj$env$parameters
    DLL            <- obj$env$DLL

    parFixed   <- SD$par.fixed
    parRandom  <- SD$par.random
    covFixed   <- SD$cov.fixed
    parNames   <- names(parFixed)
    # The hyperparameters are the outer (fixed) params that are NOT fixed effects.
    # Default to whatever the fit actually has (works for FE-only {log_tauEps},
    # BYM2 {log_tau,logit_phi,log_tauEps}, commensurate {+log_sigma_comm}, etc.)
    # rather than assuming a fixed BYM2 hyper set.
    if(is.null(hyperNames))
        hyperNames <- setdiff(unique(parNames), c("alpha","alpha_M","beta","beta_M"))
    hyperIdx   <- which(parNames %in% hyperNames)
    if(length(hyperIdx) == 0)
        stop("inlaStyleDraws: no hyperparameters found in par.fixed (have: ",
             paste(unique(parNames), collapse = ", "), ")")
    missingHN <- setdiff(hyperNames, parNames)
    if(length(missingHN))
        stop("Could not find all hyperNames in par.fixed: missing ",
             paste(missingHN, collapse = ", "))
    innerFEIdx <- setdiff(seq_along(parNames), hyperIdx)
    feNames    <- parNames[innerFEIdx]

    theta_mode <- parFixed[hyperIdx]
    fe_mode    <- parFixed[innerFEIdx]
    re_mode    <- parRandom

    # Eigen-rotation on hyper cov
    Sigma_hh <- covFixed[hyperIdx, hyperIdx, drop = FALSE]
    eig      <- eigen(Sigma_hh, symmetric = TRUE)
    if(any(eig$values <= 0)) {
        warning("hyper covariance non-PD; clamping to |eigenvalues|")
        eig$values <- abs(eig$values)
    }
    V        <- eig$vectors
    eig_sd   <- sqrt(eig$values)
    p        <- length(hyperIdx)

    z_to_theta <- function(z)
        setNames(as.numeric(theta_mode + V %*% (eig_sd * z)), hyperNames)

    # Build ONE shared TMB object for the whole walk: hyperparameters as free
    # outer params (no map), inner effects random. walkObj$fn(theta) gives the
    # Laplace-marginal NLL at hyper=theta. Reused at every CCD point so the AD
    # tape is built once, not ~19 times (see .evalCCDreuse). Inner params start
    # from paramsTemplate (the fit's MLE values) and warm-start sequentially.
    # Reuse the ORIGINAL fitted object as the walk object. It already has exactly
    # the right random-effect set, map, and a built AD tape, so walkObj$fn(theta)
    # gives the Laplace-marginal NLL at hyper=theta. The previous code rebuilt with
    # a HARDCODED BYM2 random set (and no map), which dropped real REs (alpha_M,
    # beta_M) and/or reintroduced mapped-out ones (a spurious spatial field for
    # FE-only fits) -> corrupted inner mode -> wrong conditional means.
    walkObj <- obj
    # .evalCCDreuse warm-starts by overwriting walkObj$env$last.par(.best) at each
    # config, which mutates the fitted object's inner state. Save and restore it on
    # exit so res$TMBobj is left at its MLE (a later predGrid/report reads last.par).
    .savedLP  <- walkObj$env$last.par
    .savedLPB <- if(!is.null(walkObj$env$last.par.best)) walkObj$env$last.par.best else NULL
    on.exit({
        walkObj$env$last.par <- .savedLP
        if(!is.null(.savedLPB)) walkObj$env$last.par.best <- .savedLPB
    }, add = TRUE)

    # ── Step 1: centre evaluation ─────────────────────────────────────────
    cat(sprintf("[inla] axial walk (reused tape): deltaZ=%.1f deltaPi=%.1f maxSteps=%d\n",
                deltaZ, deltaPi, maxAxialSteps))
    t0 <- proc.time()[3]
    z_centre <- rep(0, p)
    centre   <- .evalCCDreuse(walkObj, z_to_theta(z_centre))
    cat(sprintf("  centre   NLL=%.4f  (%.1f s)\n", centre$nll, proc.time()[3]-t0))
    # Centre inner mode — used to warm-start every axial point identically
    # (matches the rebuild path; see .evalCCDreuse).
    centreInner <- centre$mu

    # ── Step 2: multi-step axial walk along each eigen-axis ────────────────
    ccd <- list(centre = list(z = z_centre, nll = centre$nll,
                              mu = centre$mu, Q = centre$Q))
    axialZ_plus  <- vector("list", p); axialDN_plus  <- vector("list", p)
    axialZ_minus <- vector("list", p); axialDN_minus <- vector("list", p)

    # Walk one direction (sign = +1 or -1) along eigen-axis j; return a list
    # of (tag, z, ev, dN, abs_distance) entries for points that produced a
    # usable precision matrix. Stops on ΔNLL > δ_π, maxAxialSteps reached, or
    # the first step that yields non-finite NLL / NULL Q.
    walkAxis <- function(j, sign) {
        entries <- list()
        for(s in 1:maxAxialSteps) {
            z <- rep(0, p); z[j] <- sign * s * deltaZ
            t0 <- proc.time()[3]
            ev <- tryCatch(.evalCCDreuse(walkObj, z_to_theta(z),
                                         innerStart = centreInner),
                           error = function(e) list(nll = NaN, mu = NULL, Q = NULL,
                                                    err = conditionMessage(e)))
            dN <- ev$nll - centre$nll
            tag <- sprintf("a%d_%s%d", j, if(sign > 0) "p" else "m", s)
            elapsed <- proc.time()[3] - t0
            if(!is.finite(dN) || is.null(ev$Q)) {
                cat(sprintf("  %s  dz=%+.1f  dNLL=%s  (%.1f s) [skipped: %s]\n",
                            tag, sign*s*deltaZ,
                            if(is.finite(dN)) sprintf("%.4f", dN) else "NaN",
                            elapsed,
                            if(!is.null(ev$err)) ev$err
                            else if(!is.finite(dN)) "non-finite NLL" else "no precision"))
                break
            }
            cat(sprintf("  %s  dz=%+.1f  dNLL=%.4f  (%.1f s)\n",
                        tag, sign*s*deltaZ, dN, elapsed))
            entries[[length(entries) + 1]] <- list(tag = tag, z = z, ev = ev,
                                                   dN = dN, absStep = s * deltaZ)
            if(dN > deltaPi) break
        }
        entries
    }

    for(j in seq_len(p)) {
        for(side in c(+1, -1)) {
            entries <- walkAxis(j, side)
            for(en in entries) {
                ccd[[en$tag]] <- list(z = en$z, nll = en$ev$nll,
                                      mu = en$ev$mu, Q = en$ev$Q)
                if(side > 0) {
                    axialZ_plus[[j]]  <- c(axialZ_plus[[j]],  en$absStep)
                    axialDN_plus[[j]] <- c(axialDN_plus[[j]], en$dN)
                } else {
                    axialZ_minus[[j]]  <- c(axialZ_minus[[j]],  en$absStep)
                    axialDN_minus[[j]] <- c(axialDN_minus[[j]], en$dN)
                }
            }
        }
    }

    # ── Step 3: per-axis skew scales via least-squares ─────────────────────
    sigma_pos <- numeric(p); sigma_neg <- numeric(p)
    for(j in seq_len(p)) {
        sigma_pos[j] <- .fitSigmaLS(axialZ_plus[[j]],  axialDN_plus[[j]])
        sigma_neg[j] <- .fitSigmaLS(axialZ_minus[[j]], axialDN_minus[[j]])
    }
    cat("[inla] per-axis skew-normal scales (eigen-coordinates):\n")
    for(j in seq_len(p))
        cat(sprintf("  axis %d:  σ⁺=%.3f (from %d pts)  σ⁻=%.3f (from %d pts)  asym=%.2f\n",
                    j, sigma_pos[j], length(axialZ_plus[[j]]),
                       sigma_neg[j], length(axialZ_minus[[j]]),
                       sigma_pos[j] / sigma_neg[j]))

    # ── Step 4: sample z* from product of split Gaussians ──────────────────
    z_draws <- matrix(0, nrow = p, ncol = NDRAWS)
    for(j in seq_len(p))
        z_draws[j, ] <- .rSplitGaussian(NDRAWS, sigma_pos[j], sigma_neg[j])

    # Transform z → θ
    theta_draws <- theta_mode + V %*% (eig_sd * z_draws)
    rownames(theta_draws) <- hyperNames

    # ── Step 5: per-draw closest-CCD-point inner sample ────────────────────
    ccd_z   <- do.call(cbind, lapply(ccd, `[[`, "z"))
    ccd_names <- names(ccd)
    nInner    <- length(ccd[[1]]$mu)

    # Closest CCD for each draw via squared-distance argmin (needs only the
    # small z-coordinates, not the precisions).
    dists <- matrix(0, nrow = length(ccd), ncol = NDRAWS)
    for(k in seq_along(ccd)) {
        diff <- z_draws - ccd_z[, k]
        dists[k, ] <- colSums(diff^2)
    }
    closest <- max.col(-t(dists))    # argmin per draw

    # Sample inner ~ N(mu_k, Q_k^{-1}) via L = Cholesky(Q_k);  x = mu_k + L^{-T} z
    #
    # MEMORY: factor each Q_k lazily INSIDE this loop and discard it before the
    # next iteration, so only ONE Cholesky factor is ever alive — instead of
    # pre-building and holding a factor for all ~25 CCD points at once (the
    # accumulation that drove the multi-tens-of-GB spikes). We keep supernodal
    # (super = TRUE, the faster BLAS-3 factorization); the problem was never
    # supernodal per factor, only holding 25 of them. We also only factor the
    # CCD points that actually own draws (unique(closest)), which is often
    # fewer than all of them, and free each precision right after use.
    inner_draws <- matrix(0, nrow = nInner, ncol = NDRAWS)
    rownames(inner_draws) <- names(centre$mu)

    # A CCD point whose inner evaluation failed carries a NULL/invalid Q (and mu)
    # yet can still be the nearest point for some draws — Cholesky(NULL) then
    # throws "'n' must be a non-negative integer", killing the whole draw set
    # (this broke the auto-fallback for non-PD leave-half-out fits, e.g.
    # M_D/Delta). For each such point, redirect its draws to the NEAREST VALID
    # CCD point in eigen-coordinates; ties broken by highest density (lowest nll).
    validQ <- vapply(ccd, function(cc)
        !(is.null(cc$Q) || is.null(dim(cc$Q)) || nrow(cc$Q) < 1L || any(!is.finite(cc$mu))),
        logical(1))
    if(!any(validQ))
        stop("inlaStyleDraws: no CCD config point has a valid inner precision; ",
             "fit is degenerate (non-PD inner Hessian everywhere) — cannot draw.")
    ccd_nll <- vapply(ccd, function(cc) if(is.null(cc$nll)) Inf else cc$nll, numeric(1))
    .srcOf <- function(k) {
        if(validQ[k]) return(k)
        cand <- which(validQ)
        d2   <- colSums((ccd_z[, cand, drop = FALSE] - ccd_z[, k])^2)
        best <- cand[d2 == min(d2)]              # nearest valid point(s)
        best[which.min(ccd_nll[best])]           # tie-break: highest density (lowest nll)
    }
    usedK   <- unique(closest)
    srcForK <- vapply(usedK, .srcOf, integer(1))
    for(i in which(usedK != srcForK))
        warning(sprintf("inlaStyleDraws: CCD point %s has invalid inner precision; ",
                        "using nearest valid point %s", ccd_names[usedK[i]], ccd_names[srcForK[i]]))

    # Factor each unique SOURCE config once (memory: one factor alive at a time),
    # serving all draws routed to it.
    for(src in unique(srcForK)) {
        ksHere <- usedK[srcForK == src]
        idx    <- which(closest %in% ksHere)
        nDr    <- length(idx)
        Q      <- ccd[[src]]$Q
        muK    <- ccd[[src]]$mu
        L   <- tryCatch(Matrix::Cholesky(Q, super = TRUE), error = function(e) NULL)
        if(is.null(L)) {
            warning(sprintf("Q at CCD point %s not PD; perturbing", ccd_names[src]))
            L <- Matrix::Cholesky(Q + Matrix::Diagonal(nrow(Q), 1e-6), super = TRUE)
        }
        z   <- matrix(rnorm(nInner * nDr), nrow = nInner)
        # Reverse permutation+triangular solve. Pattern from predGrid::rmvnorm_prec
        z   <- Matrix::solve(L, z, system = "Lt")
        z   <- Matrix::solve(L, z, system = "Pt")
        inner_draws[, idx] <- muK + as.matrix(z)
        ccd[[src]]$Q <- NULL        # free this precision
        rm(Q, L, z); gc(FALSE)      # and its factor before the next point
    }

    draws <- rbind(theta_draws, inner_draws)
    attr(draws, "sigma_pos") <- sigma_pos
    attr(draws, "sigma_neg") <- sigma_neg
    attr(draws, "ccd_nll")   <- sapply(ccd, `[[`, "nll")
    attr(draws, "n_ccd")     <- length(ccd)
    attr(draws, "closest_counts") <- as.numeric(table(factor(closest, levels = seq_along(ccd))))
    draws
}
