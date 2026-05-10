#!/usr/bin/env Rscript
# 2D grid/CCD integration over (beta_normPop, log_tauEps)
# FE + nugget model, BYM2 simulated data
# INLA-style numerical integration over outer hyperparameters
#
# Outer (grid/CCD): beta_normPop, log_tauEps
# Random (inner):   alpha, beta_urban, nuggetUrbMICS, nuggetRurMICS
# Mapped out:       u_spatial, log_tau

library(TMB)
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
source("code/modM_MIIDonly.R")

cat("=== 2D Grid/CCD Integration: FE+nugget on BYM2 data ===\n\n")

# ── Data setup (same as debug_split_beta.R) ─────────────────────────
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

thisEnv = environment()
list2env(inputsMDM, envir=thisEnv)

lambdaTau = getLambdaPCprec(u=1, alpha=0.1)
lambdaTauEps = getLambdaPCprec(u=1, alpha=0.1)
nAreas = max(c(areaidxlocUrbanMICS, areaidxlocRuralMICS)) + 1
nUrb = length(ysUrbMICS)
nRur = length(ysRurMICS)

data_full = list(
  y_iUrbanMICS=ysUrbMICS, y_iRuralMICS=ysRurMICS,
  n_iUrbanMICS=nsUrbMICS, n_iRuralMICS=nsRurMICS,
  areaidxlocUrbanMICS=areaidxlocUrbanMICS, areaidxlocRuralMICS=areaidxlocRuralMICS,
  X_betaUrbanMICS=intPtsMICS$XUrb, X_betaRuralMICS=intPtsMICS$XRur,
  wUrbanMICS=intPtsMICS$wUrban, wRuralMICS=intPtsMICS$wRural,
  nAreas=nAreas,
  alpha_pri=c(0, 100^2), beta_pri=c(0, sqrt(1000)),
  lambdaTau=lambdaTau, lambdaTauEps=lambdaTauEps
)

trueAlpha = -1.25; trueUrban = 1.00; trueNormPop = 0.50
trueSigmaEps = sqrt(1.5)

initUrbP = sum(ysUrbMICS)/sum(nsUrbMICS)
initRurP = sum(ysRurMICS)/sum(nsRurMICS)
initAlpha = log(initRurP/(1-initRurP))
initBeta1 = log(initUrbP/(1-initUrbP)) - initAlpha

# ── Compile and load split template ─────────────────────────────────
unloadDynlibs()
compile("code/modM_MIIDonly_split.cpp", framework="TMBad", safebounds=FALSE)
dyn.load(dynlib("code/modM_MIIDonly_split"))
TMB::config(tmbad.sparse_hessian_compress = 1)

# ── Build TMB object ────────────────────────────────────────────────
# Initialize at previously found NM optimum (normPop≈0.89, sigmaEps≈1.225)
initLogTauEps = -2*log(1.225)  # ≈ -0.405
initNormPop = 0.89

tmb_params = list(
  log_tau=0, log_tauEps=initLogTauEps,
  alpha=initAlpha, beta_urban=initBeta1, beta_normPop=initNormPop,
  u_spatial=rep(0, nAreas),
  nuggetUrbMICS=rep(0, nUrb), nuggetRurMICS=rep(0, nRur)
)
mapFE = list(u_spatial=factor(rep(NA, nAreas)), log_tau=factor(NA))

obj = MakeADFun(data=data_full, parameters=tmb_params,
                random=c('alpha','beta_urban','nuggetUrbMICS','nuggetRurMICS'),
                map=mapFE, DLL='modM_MIIDonly_split', silent=TRUE)

cat("Outer params:", paste(names(obj$par), "=", round(obj$par, 4)), "\n")
cat(sprintf("Init: beta_normPop=%.4f, log_tauEps=%.4f\n", initNormPop, initLogTauEps))


# ════════════════════════════════════════════════════════════════════
# STEP 1: Find mode via Nelder-Mead
# ════════════════════════════════════════════════════════════════════
cat("\n── Step 1: Finding mode via Nelder-Mead ──\n")
t0_mode = proc.time()[3]
opt = optim(par=obj$par, fn=obj$fn, method="Nelder-Mead",
            control=list(trace=1, maxit=5000, reltol=1e-8))
time_mode = proc.time()[3] - t0_mode

theta_star = opt$par
nll_star = opt$value
lp = obj$env$last.par.best; rn = names(lp)

cat(sprintf("\nMode:  beta_normPop = %.4f,  log_tauEps = %.4f\n",
            theta_star["beta_normPop"], theta_star["log_tauEps"]))
cat(sprintf("       alpha = %.4f,  beta_urban = %.4f\n",
            lp[rn=="alpha"], lp[rn=="beta_urban"]))
cat(sprintf("       NLL = %.4f,  time = %.1f s\n", nll_star, time_mode))


# ════════════════════════════════════════════════════════════════════
# STEP 2: Hessian at mode → eigen-coordinates
# ════════════════════════════════════════════════════════════════════
cat("\n── Step 2: Hessian and eigendecomposition ──\n")
H = optimHess(theta_star, fn=obj$fn, gr=obj$gr)
cat("Hessian:\n"); print(round(H, 4))

eig = eigen(H)
cat("Eigenvalues:", round(eig$values, 4), "\n")
if(any(eig$values <= 0)) {
  cat("WARNING: Hessian not positive definite! Using absolute eigenvalues.\n")
  eig$values = abs(eig$values)
}

V = eig$vectors                 # rotation matrix
eig_sd = 1/sqrt(eig$values)    # SD along each eigendirection
cat("Eigen-SDs:", round(eig_sd, 4), "\n")

# Hessian-based SEs (Gaussian approx at mode)
H_inv = solve(H)
pn = names(theta_star)
cat(sprintf("SE at mode: beta_normPop = %.4f,  log_tauEps = %.4f\n",
            sqrt(H_inv[pn=="beta_normPop", pn=="beta_normPop"]),
            sqrt(H_inv[pn=="log_tauEps", pn=="log_tauEps"])))

# Transform: z → theta
z_to_theta = function(z) {
  as.vector(theta_star + V %*% (eig_sd * z))
}


# ════════════════════════════════════════════════════════════════════
# STEP 3: Find integration boundaries along eigendirections
# ════════════════════════════════════════════════════════════════════
cat("\n── Step 3: Finding boundaries (delta = 10) ──\n")
delta = 10  # NLL threshold

find_boundary = function(eig_idx, sign, step=0.5, max_z=15) {
  z = rep(0, 2)
  for(s in seq(step, max_z, by=step)) {
    z[eig_idx] = sign * s
    nll_try = obj$fn(z_to_theta(z))
    if(nll_try - nll_star > delta) return(s)
  }
  return(max_z)
}

bnd = matrix(0, 2, 2, dimnames=list(c("eig1","eig2"), c("pos","neg")))
for(k in 1:2) {
  bnd[k, "pos"] = find_boundary(k, +1)
  bnd[k, "neg"] = find_boundary(k, -1)
}
cat("Boundaries (eigen-SD units, pos / neg):\n"); print(round(bnd, 2))


# ── Helpers ──────────────────────────────────────────────────────────
wquantile = function(x, w, probs=c(0.025, 0.975)) {
  ord = order(x)
  x = x[ord]; w = w[ord]
  cw = cumsum(w) / sum(w)
  sapply(probs, function(p) x[which(cw >= p)[1]])
}

# Evaluate a set of z-points, return data.frame with NLL and inner params
eval_points = function(z_mat) {
  n = nrow(z_mat)
  pn = names(theta_star)
  res = data.frame(z1=z_mat[,1], z2=z_mat[,2],
                   beta_normPop=NA_real_, log_tauEps=NA_real_,
                   nll=NA_real_, alpha=NA_real_, beta_urban=NA_real_)

  cat(sprintf("Evaluating %d configurations...\n", n))
  t0 = proc.time()[3]
  best_nll = nll_star

  for(i in 1:n) {
    theta_i = z_to_theta(z_mat[i, ])
    res$beta_normPop[i] = theta_i[pn=="beta_normPop"]
    res$log_tauEps[i]   = theta_i[pn=="log_tauEps"]

    res$nll[i] = obj$fn(theta_i)

    lp_i = obj$env$last.par.best; rn_i = names(lp_i)
    res$alpha[i]      = as.numeric(lp_i[rn_i=="alpha"])
    res$beta_urban[i] = as.numeric(lp_i[rn_i=="beta_urban"])

    if(res$nll[i] < best_nll - 0.1) {
      cat(sprintf("  *** Point %d: NLL=%.2f < mode NLL=%.2f ***\n",
                  i, res$nll[i], best_nll))
      best_nll = res$nll[i]
    }
    if(i %% 100 == 0)
      cat(sprintf("  %d/%d  (%.0f s, %.2f s/eval)\n",
                  i, n, proc.time()[3]-t0, (proc.time()[3]-t0)/i))
  }

  elapsed = proc.time()[3] - t0
  cat(sprintf("  Done: %.0f s total, %.2f s/eval\n", elapsed, elapsed/n))
  attr(res, "time") = elapsed
  attr(res, "best_nll") = best_nll
  res
}

# Summarize weighted results
summarize = function(res, label) {
  log_w = -res$nll
  log_w = log_w - max(log_w)
  w = exp(log_w); w = w / sum(w)

  E = c(normPop    = sum(w * res$beta_normPop),
        alpha      = sum(w * res$alpha),
        urban      = sum(w * res$beta_urban),
        log_tauEps = sum(w * res$log_tauEps))

  # SD: exact for outer params (normPop, log_tauEps);
  #     across-config only for inner params (alpha, urban) — underestimates total
  SD = c(normPop    = sqrt(sum(w * (res$beta_normPop - E["normPop"])^2)),
         alpha      = sqrt(sum(w * (res$alpha - E["alpha"])^2)),
         urban      = sqrt(sum(w * (res$beta_urban - E["urban"])^2)),
         log_tauEps = sqrt(sum(w * (res$log_tauEps - E["log_tauEps"])^2)))

  CI_normPop = wquantile(res$beta_normPop, w)
  CI_alpha   = wquantile(res$alpha, w)
  CI_urban   = wquantile(res$beta_urban, w)
  n_eff = 1 / sum(w^2)

  cat(sprintf("\n── %s ──\n", label))
  cat(sprintf("  N configs: %d,  N_eff: %.1f,  Time: %.0f s\n",
              nrow(res), n_eff, attr(res, "time")))
  cat(sprintf("  normPop:     %.4f (SD=%.4f) CI[%.3f,%.3f]  bias=%+.4f\n",
              E["normPop"], SD["normPop"], CI_normPop[1], CI_normPop[2],
              E["normPop"]-0.5))
  cat(sprintf("  alpha:       %.4f (SD*=%.4f) CI*[%.3f,%.3f] bias=%+.4f\n",
              E["alpha"], SD["alpha"], CI_alpha[1], CI_alpha[2],
              E["alpha"]-(-1.25)))
  cat(sprintf("  beta_urban:  %.4f (SD*=%.4f) CI*[%.3f,%.3f] bias=%+.4f\n",
              E["urban"], SD["urban"], CI_urban[1], CI_urban[2],
              E["urban"]-1.0))
  cat(sprintf("  log_tauEps:  %.4f (SD=%.4f)\n", E["log_tauEps"], SD["log_tauEps"]))
  cat(sprintf("  sigmaEps:    %.4f (truth=%.4f)\n", exp(-0.5*E["log_tauEps"]), sqrt(1.5)))
  cat("  (*SD/CI for alpha,urban = across-config only, underestimates total)\n")

  invisible(list(E=E, SD=SD, CI_normPop=CI_normPop, CI_alpha=CI_alpha,
                 CI_urban=CI_urban, n_eff=n_eff, w=w, time=attr(res, "time")))
}


# ════════════════════════════════════════════════════════════════════
# METHOD A: Regular grid in eigen-coordinates (25 × 25)
# ════════════════════════════════════════════════════════════════════
cat("\n════════════════════════════════════════════════════════════\n")
cat("METHOD A: Regular grid (25 x 25) in eigen-coordinates\n")
cat("════════════════════════════════════════════════════════════\n")

n_per_dim = 25
z1_seq = seq(-bnd["eig1","neg"], bnd["eig1","pos"], length.out=n_per_dim)
z2_seq = seq(-bnd["eig2","neg"], bnd["eig2","pos"], length.out=n_per_dim)
z_grid = as.matrix(expand.grid(z1=z1_seq, z2=z2_seq))

res_grid = eval_points(z_grid)
summ_grid = summarize(res_grid, "Regular grid 25x25")


# ════════════════════════════════════════════════════════════════════
# METHOD B: CCD (8 ray directions, step until threshold)
# ════════════════════════════════════════════════════════════════════
cat("\n════════════════════════════════════════════════════════════\n")
cat("METHOD B: CCD (8 directions, stepped radii)\n")
cat("════════════════════════════════════════════════════════════\n")

# 8 CCD directions (unit vectors in z-space)
ccd_dirs = rbind(
  c( 1,  0), c(-1,  0),                       # axial z1
  c( 0,  1), c( 0, -1),                       # axial z2
  c( 1,  1)/sqrt(2), c( 1, -1)/sqrt(2),       # diagonals
  c(-1,  1)/sqrt(2), c(-1, -1)/sqrt(2)
)

step_ccd = 0.5
max_r_ccd = max(bnd) + 2

# Evaluate center
pn = names(theta_star)
theta_c = z_to_theta(c(0, 0))
nll_c = obj$fn(theta_c)
lp_c = obj$env$last.par.best; rn_c = names(lp_c)

res_ccd = data.frame(z1=0, z2=0,
                     beta_normPop=theta_c[pn=="beta_normPop"],
                     log_tauEps=theta_c[pn=="log_tauEps"],
                     nll=nll_c,
                     alpha=as.numeric(lp_c[rn_c=="alpha"]),
                     beta_urban=as.numeric(lp_c[rn_c=="beta_urban"]))

cat(sprintf("Stepping along 8 CCD directions (step=%.1f, delta=%d)...\n",
            step_ccd, delta))
t0_ccd = proc.time()[3]
best_nll_ccd = nll_star

for(d in 1:nrow(ccd_dirs)) {
  for(r in seq(step_ccd, max_r_ccd, by=step_ccd)) {
    z = r * ccd_dirs[d, ]
    theta_i = z_to_theta(z)
    nll_i = obj$fn(theta_i)
    lp_i = obj$env$last.par.best; rn_i = names(lp_i)

    res_ccd = rbind(res_ccd, data.frame(
      z1=z[1], z2=z[2],
      beta_normPop=theta_i[pn=="beta_normPop"],
      log_tauEps=theta_i[pn=="log_tauEps"],
      nll=nll_i,
      alpha=as.numeric(lp_i[rn_i=="alpha"]),
      beta_urban=as.numeric(lp_i[rn_i=="beta_urban"])))

    if(nll_i < best_nll_ccd - 0.1) {
      cat(sprintf("  *** Dir %d, r=%.1f: NLL=%.2f < mode NLL=%.2f ***\n",
                  d, r, nll_i, best_nll_ccd))
      best_nll_ccd = nll_i
    }
    if(nll_i - nll_star > delta) break
  }
}
time_ccd = proc.time()[3] - t0_ccd
attr(res_ccd, "time") = time_ccd
attr(res_ccd, "best_nll") = best_nll_ccd

cat(sprintf("CCD: %d configurations in %.0f s (%.2f s/eval)\n",
            nrow(res_ccd), time_ccd, time_ccd/nrow(res_ccd)))

summ_ccd = summarize(res_ccd, "CCD (8 directions)")


# ════════════════════════════════════════════════════════════════════
# COMPARISON TABLE
# ════════════════════════════════════════════════════════════════════
cat("\n\n════════════════════════════════════════════════════════════\n")
cat("COMPARISON TABLE (FE + nugget, BYM2 simulated data)\n")
cat("════════════════════════════════════════════════════════════\n\n")

cat(sprintf("%-30s %10s %10s %10s %10s\n",
            "Method", "normPop", "alpha", "urban", "Time(s)"))
cat(strrep("-", 72), "\n")
cat(sprintf("%-30s %10s %10s %10s %10s\n", "Truth", "0.50", "-1.25", "1.00", "--"))
cat(sprintf("%-30s %10s %10s %10s %10s\n", "All random, BFGS",
            "0.13", "-0.94", "1.30", "103"))
cat(sprintf("%-30s %10s %10s %10s %10s\n", "normPop outer, NM",
            "0.89", "-1.68", "0.63", "239"))
cat(sprintf("%-30s %10s %10s %10s %10s\n", "MCMC laplace=FALSE",
            "0.68", "-1.52", "0.86", "24127"))
cat(sprintf("%-30s %10.4f %10.4f %10.4f %10.0f\n",
            sprintf("Grid 25x25 (%d pts)", nrow(res_grid)),
            summ_grid$E["normPop"], summ_grid$E["alpha"],
            summ_grid$E["urban"], summ_grid$time + time_mode))
cat(sprintf("%-30s %10.4f %10.4f %10.4f %10.0f\n",
            sprintf("CCD (%d pts)", nrow(res_ccd)),
            summ_ccd$E["normPop"], summ_ccd$E["alpha"],
            summ_ccd$E["urban"], summ_ccd$time + time_mode))

cat("\nUncertainty comparison:\n")
cat(sprintf("  %-28s normPop SD=0.09 [0.50,0.86]  alpha SD=0.09  urban SD=0.13\n",
            "MCMC laplace=FALSE"))
cat(sprintf("  %-28s normPop SD=%.2f [%.2f,%.2f]   alpha SD*=%.2f  urban SD*=%.2f\n",
            "Grid 25x25",
            summ_grid$SD["normPop"], summ_grid$CI_normPop[1], summ_grid$CI_normPop[2],
            summ_grid$SD["alpha"], summ_grid$SD["urban"]))
cat(sprintf("  %-28s normPop SD=%.2f [%.2f,%.2f]   alpha SD*=%.2f  urban SD*=%.2f\n",
            "CCD",
            summ_ccd$SD["normPop"], summ_ccd$CI_normPop[1], summ_ccd$CI_normPop[2],
            summ_ccd$SD["alpha"], summ_ccd$SD["urban"]))
cat("  (*alpha/urban SDs from grid/CCD = across-config only, underestimate total)\n")


# ════════════════════════════════════════════════════════════════════
# PLOTS
# ════════════════════════════════════════════════════════════════════
dir.create("Figures/gridIntegration", showWarnings=FALSE, recursive=TRUE)

# Plot 1: NLL surface + density + CCD points
pdf("Figures/gridIntegration/nll_surface_grid.pdf", width=14, height=5)
par(mfrow=c(1,3), mar=c(4,4,3,1))

# (a) NLL surface from grid
nll_mat = matrix(res_grid$nll, nrow=n_per_dim)
nll_shifted = nll_mat - min(nll_mat, na.rm=TRUE)
image(z1_seq, z2_seq, pmin(nll_shifted, 20),
      main="NLL - NLL_min (grid, eigen-coords)",
      xlab="z1 (eig-SD)", ylab="z2 (eig-SD)",
      col=hcl.colors(50, "YlOrRd", rev=TRUE))
contour(z1_seq, z2_seq, nll_shifted,
        levels=c(1, 2, 5, 10, 15, 20), add=TRUE)
points(0, 0, pch=4, cex=2, lwd=2)

# (b) Relative density from grid
log_w = -res_grid$nll; log_w = log_w - max(log_w)
w_mat = matrix(exp(log_w), nrow=n_per_dim)
image(z1_seq, z2_seq, w_mat/max(w_mat),
      main="Relative density (grid)",
      xlab="z1 (eig-SD)", ylab="z2 (eig-SD)",
      col=hcl.colors(50, "Blues", rev=TRUE))
contour(z1_seq, z2_seq, w_mat/max(w_mat),
        levels=c(0.001, 0.01, 0.05, 0.1, 0.5), add=TRUE)

# (c) CCD points overlaid on NLL contours
plot(res_ccd$z1, res_ccd$z2, pch=16, cex=0.8,
     main="CCD configurations",
     xlab="z1 (eig-SD)", ylab="z2 (eig-SD)",
     xlim=range(z1_seq), ylim=range(z2_seq))
contour(z1_seq, z2_seq, nll_shifted,
        levels=c(1, 2, 5, 10, 15, 20), add=TRUE, col="grey60")
points(0, 0, pch=4, cex=2, lwd=2, col="red")
dev.off()

# Plot 2: Marginal posterior of normPop
pdf("Figures/gridIntegration/marginal_normPop.pdf", width=8, height=6)

# Kernel density from grid weights
w_grid = exp(-res_grid$nll - max(-res_grid$nll))
w_grid = w_grid / sum(w_grid)
np_seq = seq(min(res_grid$beta_normPop) - 0.1,
             max(res_grid$beta_normPop) + 0.1, length.out=300)
bw = 0.03
dens_grid = sapply(np_seq, function(x)
  sum(w_grid * dnorm(res_grid$beta_normPop, x, bw)))
dens_grid = dens_grid / (sum(dens_grid) * diff(np_seq)[1])

# Kernel density from CCD weights
w_ccd_v = exp(-res_ccd$nll - max(-res_ccd$nll))
w_ccd_v = w_ccd_v / sum(w_ccd_v)
dens_ccd = sapply(np_seq, function(x)
  sum(w_ccd_v * dnorm(res_ccd$beta_normPop, x, bw)))
dens_ccd = dens_ccd / (sum(dens_ccd) * diff(np_seq)[1])

plot(np_seq, dens_grid, type='l', lwd=2, col='blue',
     main="Marginal posterior of normPop",
     xlab=expression(beta[normPop]), ylab="Density",
     ylim=c(0, max(c(dens_grid, dens_ccd))*1.1))
lines(np_seq, dens_ccd, lwd=2, col='darkgreen', lty=2)
abline(v=0.5, col='red', lwd=2, lty=2)
abline(v=summ_grid$E["normPop"], col='blue', lty=3)
abline(v=summ_ccd$E["normPop"], col='darkgreen', lty=3)
legend("topright", c("Grid marginal", "CCD marginal", "Truth (0.50)"),
       col=c("blue", "darkgreen", "red"), lty=c(1, 2, 2), lwd=2)
dev.off()

# Plot 3: NLL surface in original (normPop, log_tauEps) coordinates
pdf("Figures/gridIntegration/nll_surface_original_coords.pdf", width=10, height=5)
par(mfrow=c(1,2), mar=c(4,4,3,1))

# Grid
np_u = sort(unique(round(res_grid$beta_normPop, 6)))
lt_u = sort(unique(round(res_grid$log_tauEps, 6)))
plot(res_grid$beta_normPop, res_grid$log_tauEps, pch=16, cex=0.4,
     col=hcl.colors(50, "YlOrRd", rev=TRUE)[pmin(50, ceiling(
       (res_grid$nll - min(res_grid$nll)) / max(1, max(res_grid$nll)-min(res_grid$nll)) * 49) + 1)],
     main="Grid points (original coords)",
     xlab=expression(beta[normPop]), ylab=expression(log(tau[epsilon])))
abline(v=0.5, col='red', lty=2)
points(theta_star["beta_normPop"], theta_star["log_tauEps"],
       pch=4, cex=2, lwd=2)

# CCD
plot(res_ccd$beta_normPop, res_ccd$log_tauEps, pch=16, cex=0.8,
     col=hcl.colors(50, "YlOrRd", rev=TRUE)[pmin(50, ceiling(
       (res_ccd$nll - min(res_ccd$nll)) / max(1, max(res_ccd$nll)-min(res_ccd$nll)) * 49) + 1)],
     main="CCD points (original coords)",
     xlab=expression(beta[normPop]), ylab=expression(log(tau[epsilon])))
abline(v=0.5, col='red', lty=2)
points(theta_star["beta_normPop"], theta_star["log_tauEps"],
       pch=4, cex=2, lwd=2)
dev.off()

cat("\nPlots saved to Figures/gridIntegration/\n")

# ── Save ─────────────────────────────────────────────────────────────
save(res_grid, summ_grid, res_ccd, summ_ccd,
     theta_star, nll_star, H, eig, bnd, time_mode,
     file="savedOutput/gridIntegration_FEnug_2D.RData")
cat("Results saved to savedOutput/gridIntegration_FEnug_2D.RData\n")
cat("\nDone.\n")
