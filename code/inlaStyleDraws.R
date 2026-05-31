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

# -- internal: rebuild MakeADFun with the three hypers held fixed at fixedHyper
#    and the four inner blocks declared random. Uses the data/parameters from
#    the original obj.
.makeFixedHyperObj <- function(dataList, paramsTemplate, fixedHyper, DLL,
                               innerWarm = NULL, verbose = FALSE) {
    pars <- paramsTemplate
    if(!is.null(innerWarm)) for(nm in names(innerWarm)) pars[[nm]] <- innerWarm[[nm]]
    pars$log_tau    <- as.numeric(fixedHyper["log_tau"])
    pars$logit_phi  <- as.numeric(fixedHyper["logit_phi"])
    pars$log_tauEps <- as.numeric(fixedHyper["log_tauEps"])
    mapList <- list(log_tau    = factor(NA),
                    logit_phi  = factor(NA),
                    log_tauEps = factor(NA))
    MakeADFun(data = dataList, parameters = pars,
              random = c("alpha", "beta", "w_bym2Free", "u_bym2Free"),
              map    = mapList,
              DLL    = DLL, silent = !verbose)
}

# -- internal: evaluate NLL and cache (μ, Q) at one CCD point.
.evalCCD <- function(dataList, paramsTemplate, theta_k, DLL,
                     innerWarm = NULL, verbose = FALSE) {
    obj_k <- .makeFixedHyperObj(dataList, paramsTemplate, theta_k, DLL,
                                innerWarm = innerWarm, verbose = verbose)
    nll <- as.numeric(obj_k$fn(numeric()))
    mu  <- obj_k$env$last.par.best
    # Hessian of joint NLL w.r.t. random parameters at the inner mode.
    # This IS the precision of the Gaussian approximation of p(x | θ_k, y).
    Q   <- tryCatch(obj_k$env$spHess(mu, random = TRUE),
                    error = function(e) NULL)
    if(is.null(Q)) {
        # Fallback: use sdreport's jointPrecision
        SD <- tryCatch(TMB::sdreport(obj_k, getJointPrecision = TRUE),
                       error = function(e) NULL)
        if(!is.null(SD)) Q <- SD$jointPrecision
    }
    list(nll = nll, mu = mu, Q = Q, obj = obj_k)
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
posteriorDraws <- function(res, NDRAWS = 1000, useInla = "auto", ...) {
    mode <- if(identical(useInla, "auto"))           "auto"
            else if(isTRUE(useInla)  || identical(useInla, "yes")) "inla"
            else if(isFALSE(useInla) || identical(useInla, "no"))  "gauss"
            else stop("posteriorDraws: useInla must be one of \"auto\", TRUE/\"yes\", or FALSE/\"no\"")

    if(mode == "inla")
        return(inlaStyleDraws(res, NDRAWS = NDRAWS, ...))

    if(mode == "gauss")
        return(.gaussianPosteriorDraws(res, NDRAWS = NDRAWS,
                                       requireJoint = TRUE))

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
    if(!is.null(out)) return(out)
    inlaStyleDraws(res, NDRAWS = NDRAWS, ...)
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
                            hyperNames = c("log_tau","logit_phi","log_tauEps"),
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
    hyperIdx   <- which(parNames %in% hyperNames)
    if(length(hyperIdx) != length(hyperNames))
        stop("Could not find all hyperNames in par.fixed: missing ",
             paste(setdiff(hyperNames, parNames), collapse = ", "))
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

    # ── Step 1: centre evaluation ─────────────────────────────────────────
    cat(sprintf("[inla] axial walk: deltaZ=%.1f deltaPi=%.1f maxSteps=%d\n",
                deltaZ, deltaPi, maxAxialSteps))
    t0 <- proc.time()[3]
    z_centre <- rep(0, p)
    centre   <- .evalCCD(dataList, paramsTemplate, z_to_theta(z_centre), DLL,
                         innerWarm = NULL, verbose = verbose)
    cat(sprintf("  centre   NLL=%.4f  (%.1f s)\n", centre$nll, proc.time()[3]-t0))

    # Use centre's inner mode as the warm start
    extractInnerWarm <- function(mu) {
        list(alpha       = as.numeric(mu[names(mu) == "alpha"]),
             beta        = as.numeric(mu[names(mu) == "beta"]),
             w_bym2Free  = as.numeric(mu[names(mu) == "w_bym2Free"]),
             u_bym2Free  = as.numeric(mu[names(mu) == "u_bym2Free"]))
    }
    innerWarm <- extractInnerWarm(centre$mu)

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
            ev <- tryCatch(.evalCCD(dataList, paramsTemplate,
                                    z_to_theta(z), DLL,
                                    innerWarm = innerWarm, verbose = verbose),
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

    # Pre-Cholesky each Q_k
    L_list <- vector("list", length(ccd))
    for(k in seq_along(ccd)) {
        Q <- ccd[[k]]$Q
        L_list[[k]] <- tryCatch(Matrix::Cholesky(Q, super = TRUE),
                                error = function(e) NULL)
        if(is.null(L_list[[k]])) {
            warning(sprintf("Q at CCD point %s not PD; perturbing", ccd_names[k]))
            n <- nrow(Q)
            Qp <- Q + Matrix::Diagonal(n, 1e-6)
            L_list[[k]] <- Matrix::Cholesky(Qp, super = TRUE)
        }
    }

    # Closest CCD for each draw via squared-distance argmin
    dists <- matrix(0, nrow = length(ccd), ncol = NDRAWS)
    for(k in seq_along(ccd)) {
        diff <- z_draws - ccd_z[, k]
        dists[k, ] <- colSums(diff^2)
    }
    closest <- max.col(-t(dists))    # argmin per draw

    # Sample inner ~ N(mu_k, Q_k^{-1}) via L = Cholesky(Q_k);  x = mu_k + L^{-T} z
    inner_draws <- matrix(0, nrow = nInner, ncol = NDRAWS)
    rownames(inner_draws) <- names(centre$mu)
    for(k in unique(closest)) {
        idx <- which(closest == k)
        nDr <- length(idx)
        z   <- matrix(rnorm(nInner * nDr), nrow = nInner)
        # Reverse permutation+triangular solve. Pattern from predGrid::rmvnorm_prec
        z   <- Matrix::solve(L_list[[k]], z, system = "Lt")
        z   <- Matrix::solve(L_list[[k]], z, system = "Pt")
        z   <- as.matrix(z)
        inner_draws[, idx] <- ccd[[k]]$mu + z
    }

    draws <- rbind(theta_draws, inner_draws)
    attr(draws, "sigma_pos") <- sigma_pos
    attr(draws, "sigma_neg") <- sigma_neg
    attr(draws, "ccd_nll")   <- sapply(ccd, `[[`, "nll")
    attr(draws, "n_ccd")     <- length(ccd)
    attr(draws, "closest_counts") <- as.numeric(table(factor(closest, levels = seq_along(ccd))))
    draws
}
