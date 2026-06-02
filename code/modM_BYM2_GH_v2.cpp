// ============================================================================
// modM_BYM2_GH_v2.cpp  —  MICS-only BYM2 model with Gauss-Hermite nugget
// integration, optimised AD-tape (v2 of modM_BYM2_GH.cpp).
//
// MODEL
// -----
//   y_i | n_i, p_i        ~ Binomial(n_i, p_i)
//   logit p_i              = alpha + x_i^T beta + b_{a(i)} + eps_i
//   b                      ~ BYM2(tau, phi) on `nAreas` admin areas
//   eps_i                  ~ N(0, 1/tau_eps)         (cluster-level nugget)
//
// The BYM2 prior decomposes b = (1/sqrt(tau)) * [sqrt(1-phi) v + sqrt(phi) u]
// with v iid N(0,I) and u ~ ICAR(Q). We parameterise the latent field as
//   w  = b   (the BYM2 effect itself, length n)
//   u  = ICAR component (length n)
// with sum-to-zero constraints on BOTH w and u (n-1 free dims each, so 2n-2
// free params total). The iid component v is implicit (solved from w and u).
//
// We integrate the cluster-level nugget eps_i analytically out of the
// likelihood at each obs by Q-point Gauss-Hermite quadrature, so the inner
// random effects passed to TMB Laplace are only (alpha, beta, w, u).
//
// OPTIMISATIONS vs modM_BYM2_GH.cpp ("v1")
// ----------------------------------------
//   - base_eta computed once (not Q times inside the GH loop)
//   - INLINE binomial likelihood instead of dbinom_robust (see note below)
//   - lchoose factored out as DATA (constant — not on the AD tape)
//   - no if-branching on AD tape (use CppAD::CondExpGt for the argmax)
//   - logSumExp done row-wise on an nObs x Q matrix
//
// NOTE 1 — why not dbinom_robust?
//   TMB's dbinom_robust(y, size, logit_p, give_log) returns the log-binomial
//   PMF computed in a numerically stable way: log_choose(size, y) +
//   y*logit_p - size*log1pexp(logit_p). It's the safe drop-in for a single
//   binomial observation. We DON'T use it here because:
//     (a) Inside the GH nugget integration we need raw log-binomial values
//         for each (obs, integration-point, GH-node) triple, summed into a
//         mixture before being logged again. Calling dbinom_robust per triple
//         adds function-call overhead on the AD tape with no extra stability:
//         we already inline the same numerically-stable formula
//             log_binom = y*eta - n*logspace_add(0, eta)
//         (the lchoose term is constant w.r.t. params and we add it once
//         per obs, outside the GH/integration loops, as DATA).
//     (b) dbinom_robust internally branches and validates inputs, both of
//         which bloat the TMB AD tape.
//     (c) Factoring out lchoose into a DATA_VECTOR keeps it off the AD tape
//         entirely (it's a fixed function of (y, n)).
//   Net result: the inline form generates a smaller, faster tape with
//   identical numerics to dbinom_robust at the points where stability matters.
//
// NOTE 2 — what is logspace_add?
//   logspace_add(a, b) computes log(exp(a) + exp(b)) without intermediate
//   overflow. The standard implementation is:
//       max(a, b) + log1p(exp(-|a - b|))
//   so when a >> b (or vice versa) you don't overflow exp(max). The single
//   call we use is `logspace_add(0, eta)`, which equals log(1 + exp(eta)) —
//   i.e. the binomial log-normaliser, sometimes written log1pexp(eta).
//   With this, `n * logspace_add(0, eta)` is the "n * log(1 + exp(eta))"
//   term in the binomial logit likelihood, stable for any sign of eta.
//
// NOTE 3 — Gauss-Hermite for the nugget
//   We need int_R exp(log_binom(eta + eps)) * N(eps | 0, sigma_eps^2) deps.
//   Substituting eps = sqrt(2)*sigma_eps*z gives the standard GH form
//       (1/sqrt(pi)) * sum_q w_q * exp(log_binom(eta + sqrt(2)*sigma_eps*z_q))
//   where {w_q, z_q} are the Q-point Gauss-Hermite weights/nodes. We then
//   take log of the sum stably via logSumExp (see the matrix construction
//   and argmax loops below).
// ============================================================================

#include <TMB.hpp>
#include <Eigen/Sparse>
using namespace density;
using Eigen::SparseMatrix;

// PC prior on precision (log scale)
template<class Type>
Type dPCPriTau(Type logTau, Type lambda)
{
  Type tau = exp(logTau);
  Type ldensity = log(lambda) - log(2.0) - Type(1.5)*logTau - lambda/sqrt(tau);
  Type ljacobian = logTau;
  return ldensity + ljacobian;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // ── Data ──
  DATA_VECTOR( y_iUrbanMICS );
  DATA_VECTOR( y_iRuralMICS );
  DATA_VECTOR( n_iUrbanMICS );
  DATA_VECTOR( n_iRuralMICS );
  DATA_IVECTOR(areaidxlocUrbanMICS);
  DATA_IVECTOR(areaidxlocRuralMICS);
  DATA_MATRIX( X_betaUrbanMICS );
  DATA_MATRIX( X_betaRuralMICS );
  DATA_ARRAY( wUrbanMICS );
  DATA_ARRAY( wRuralMICS );

  DATA_SPARSE_MATRIX( Q_bym2 );

  // Precomputed log-binomial-coefficient (DATA — not on AD tape)
  DATA_VECTOR( lchoose_urban );
  DATA_VECTOR( lchoose_rural );

  // GH quadrature nodes/weights
  DATA_VECTOR( gh_nodes );
  DATA_VECTOR( gh_weights );

  // Prior parameters
  DATA_VECTOR( alpha_pri );
  DATA_VECTOR( beta_pri );

  // BYM2 precomputed
  DATA_SCALAR( tr );
  DATA_VECTOR( gammaTildesm1 );

  DATA_SCALAR( lambdaPhi );
  DATA_INTEGER( uniformPhiPrior );
  DATA_SCALAR( lambdaTau );
  DATA_SCALAR( lambdaTauEps );

  DATA_SCALAR( options );

  // ── Parameters ──
  // Outer (hyperparameters)
  PARAMETER( log_tau );
  PARAMETER( logit_phi );
  PARAMETER( log_tauEps );

  // Inner (random effects + fixed effects)
  PARAMETER( alpha );
  PARAMETER_VECTOR( beta );

  // BYM2 spatial: n-1 free each, nth = -sum(free)
  PARAMETER_VECTOR( w_bym2Free );
  PARAMETER_VECTOR( u_bym2Free );

  // ── Dimensions ──
  int nUrb = y_iUrbanMICS.size();
  int nRur = y_iRuralMICS.size();
  int nFree = w_bym2Free.size();
  int nAreas = nFree + 1;
  int KUrb = wUrbanMICS.cols();
  int KRur = wRuralMICS.cols();
  int Q = gh_nodes.size();

  // ── Transforms ──
  Type tau = exp(log_tau);
  Type phi = Type(1.0)/(Type(1.0) + exp(-logit_phi));
  Type sigmaEps = exp(Type(-0.5) * log_tauEps);

  // Reconstruct full n-vectors with sum-to-zero
  vector<Type> w_bym2(nAreas);
  vector<Type> u_bym2(nAreas);
  Type wSum = Type(0.0);
  Type uSum = Type(0.0);
  for(int i = 0; i < nFree; i++) {
    w_bym2(i) = w_bym2Free(i);
    u_bym2(i) = u_bym2Free(i);
    wSum += w_bym2Free(i);
    uSum += u_bym2Free(i);
  }
  w_bym2(nFree) = -wSum;
  u_bym2(nFree) = -uSum;

  Type jnll = Type(0.0);

  // ═══════════════════════════════════
  // (1) BYM2 GMRF density
  // ═══════════════════════════════════
  // Joint log density of (w, u) under the BYM2 parameterisation. With
  //   b = (1/sqrt(tau)) [sqrt(1-phi) v + sqrt(phi) u]
  // solving for v gives v = [sqrt(tau)*b - sqrt(phi)*u] / sqrt(1-phi), so
  //   v^T v  =  [tau b^T b  - 2 sqrt(tau*phi) u^T b  +  phi u^T u] / (1-phi)
  // and the log density splits into the three quadratic pieces below plus
  // an ICAR term u^T Q u from p(u). `gammaTildesm1` are the n-1 non-zero
  // eigenvalues of the scaled ICAR precision minus 1; we use them for the
  // PC prior on phi (the log-det piece is absorbed by the data).

  // quadW = tau/(1-phi) * w^T w   (contribution to v^T v from the iid term)
  Type quadW = (w_bym2 * w_bym2).sum() * tau / (Type(1.0) - phi);

  // quadU = u^T Q u + phi/(1-phi) u^T u  (ICAR + iid cross-term on u)
  vector<Type> Qu = (Q_bym2 * u_bym2.matrix()).col(0);
  Type fac = phi / (Type(1.0) - phi);
  Type quadU = (Qu * u_bym2).sum() + fac * (u_bym2 * u_bym2).sum();

  // quadWU = -2 sqrt(phi*tau)/(1-phi) * u^T w  (cross-term from v^T v)
  Type quadWU = (u_bym2 * w_bym2).sum() * (-Type(2.0) * sqrt(phi * tau) / (Type(1.0) - phi));

  Type quadSum = quadW + quadU + quadWU;        // = v^T v  + u^T Q u

  // Log determinant of the scaled ICAR precision at the current phi.
  // sum_i log(1 + phi*(gamma_i - 1))  over the n-1 non-zero eigenvalues.
  // (Used by the PC prior on phi below — NOT added into bym2LogLik because
  // it's independent of the latent field and gets absorbed by the prior.)
  Type logDet = Type(0.0);
  for(int i = 0; i < gammaTildesm1.size(); i++) {
    logDet += log(Type(1.0) + phi * gammaTildesm1(i));
  }
  // Tau/phi-dependent normaliser from the v^T v change of variables:
  // log|det dv/dw| = (n-1)/2 * log(tau/(1-phi))    (sum-to-zero subspace)
  Type logDetTau = Type(nFree) * log((Type(1.0) - phi) / tau);

  // jnll <- jnll - log p(w, u | tau, phi); the (n-1)/2 * log|Q*| term is a
  // constant in tau, phi and drops.
  Type bym2LogLik = Type(-0.5) * logDetTau + Type(-0.5) * quadSum;
  jnll -= bym2LogLik;

  // ═══════════════════════════════════
  // (2) Priors
  // ═══════════════════════════════════

  // PC prior for tau
  Type logPriTau = log(lambdaTau) - log(Type(2.0)) - Type(1.5) * log_tau
                   - lambdaTau / sqrt(tau) + log_tau;
  jnll -= logPriTau;

  // PC prior for tauEps
  jnll -= dPCPriTau(log_tauEps, lambdaTauEps);

  // Prior for phi
  if(uniformPhiPrior) {
    // Uniform prior on phi in [0,1]: Jacobian for logit transform only
    jnll -= -logit_phi - Type(2.0) * log(Type(1.0) + exp(-logit_phi));
  } else {
    // PC prior for phi
    Type n_gam = Type(gammaTildesm1.size());
    Type KLD = Type(0.5) * (phi * tr - phi * n_gam - logDet);
    Type d = sqrt(Type(2.0) * KLD);
    Type lexpDensity = log(lambdaPhi) - lambdaPhi * d;
    Type sumVal = Type(0.0);
    for(int i = 0; i < gammaTildesm1.size(); i++) {
      sumVal += gammaTildesm1(i) / (Type(1.0) + gammaTildesm1(i) * phi);
    }
    Type ljacobian = -log(d) - log(Type(2.0)) + log(tr - n_gam - sumVal);
    ljacobian = ljacobian - logit_phi - Type(2.0) * log(Type(1.0) + exp(-logit_phi));
    jnll -= lexpDensity + ljacobian;
  }

  // Prior for intercept
  jnll -= dnorm(alpha, alpha_pri(0), alpha_pri(1), true);

  // Prior for covariates
  for(int i = 0; i < beta.size(); i++) {
    jnll -= dnorm(beta(i), beta_pri(0), beta_pri(1), true);
  }

  // ═══════════════════════════════════
  // (3) GH-integrated likelihood (vectorized)
  // ═══════════════════════════════════
  // For each obs i with cluster nugget eps_i ~ N(0, sigma_eps^2) we evaluate
  //     L_i = E_{eps} [ Binom(y_i; n_i, expit(eta_i + eps_i)) ]
  // with eta_i = alpha + x_i^T beta + b_{a(i)}. The GH change of variables
  // eps = sqrt(2)*sigma_eps*z gives
  //     L_i = (1/sqrt(pi)) * sum_q w_q * Binom(y_i; n_i, expit(eta_i + sqrt(2)*sigma_eps*z_q))
  // so log L_i = log(1/sqrt(pi)) + logSumExp_q { log(w_q) + log_binom_iq }.
  // For obs whose linear predictor itself is a mixture over K integration
  // points (positional / population integration), the inner sum becomes
  //     log( sum_k w_ik * exp(log_binom_iqk) )
  // which is what the inner K-loop computes.

  // Precompute linear-predictor pieces that don't depend on q (GH node) or
  // eps. fe_* is laid out as length nObs*K, stacked by integration point k.
  vector<Type> fe_urb = X_betaUrbanMICS * beta;
  vector<Type> fe_rur = X_betaRuralMICS * beta;

  // GH abscissae:   eps_q = sqrt(2) * sigma_eps * node_q     (Q values)
  vector<Type> eps_gh(Q);
  for(int q = 0; q < Q; q++) {
    eps_gh(q) = sqrt(Type(2.0)) * sigmaEps * gh_nodes(q);
  }
  Type log_inv_sqrt_pi = -Type(0.5) * log(M_PI);   // = log(1/sqrt(pi))

  // Cache log(w_q) once
  vector<Type> log_gh_w(Q);
  for(int q = 0; q < Q; q++) {
    log_gh_w(q) = log(gh_weights(q));
  }

  // ── base_eta_{i,k} = alpha + x_{ik}^T beta + b_{a(i,k)}  (computed ONCE,
  //    not Q times inside the GH loop). This is the dominant speedup vs v1.

  matrix<Type> base_eta_urb(nUrb, KUrb);
  for(int k = 0; k < KUrb; k++) {
    for(int i = 0; i < nUrb; i++) {
      base_eta_urb(i, k) = alpha + fe_urb(nUrb * k + i)
                            + w_bym2(areaidxlocUrbanMICS(i));
    }
  }

  matrix<Type> base_eta_rur(nRur, KRur);
  for(int k = 0; k < KRur; k++) {
    for(int i = 0; i < nRur; i++) {
      base_eta_rur(i, k) = alpha + fe_rur(nRur * k + i)
                            + w_bym2(areaidxlocRuralMICS(i));
    }
  }

  // ── Urban: for each (obs i, GH node q), compute log( sum_k w_ik * lik )
  //          and store as the (i, q) cell of logMixUrb. Then logSumExp over q.
  matrix<Type> logMixUrb(nUrb, Q);

  for(int q = 0; q < Q; q++) {
    for(int i = 0; i < nUrb; i++) {
      Type mix_lik = Type(0.0);
      Type yi = y_iUrbanMICS(i);
      Type ni = n_iUrbanMICS(i);

      for(int k = 0; k < KUrb; k++) {
        Type w_ik = wUrbanMICS(i, k);
        Type eta = base_eta_urb(i, k) + eps_gh(q);          // logit(p_ikq)

        // INLINE binomial log-likelihood (no lchoose; added once below):
        //   log_binom = y*eta - n*log(1 + exp(eta))
        // We use logspace_add(0, eta) = log(1 + exp(eta)) for numerical
        // stability when eta is large positive (naive exp(eta) would
        // overflow). See header NOTE 1 (why not dbinom_robust) and NOTE 2
        // (what logspace_add is).
        Type log_binom = yi * eta - ni * logspace_add(Type(0), eta);
        mix_lik += w_ik * exp(log_binom);                   // sum over int pts
      }

      // Add log(GH weight) and the constant lchoose(n_i, y_i) once per obs.
      // lchoose is DATA so it's not on the AD tape.
      logMixUrb(i, q) = log_gh_w(q) + lchoose_urban(i) + log(mix_lik);
    }
  }

  // logSumExp over Q nodes (AD-tape-friendly: CondExpGt instead of `if`).
  for(int i = 0; i < nUrb; i++) {
    Type max_val = logMixUrb(i, 0);
    for(int q = 1; q < Q; q++) {
      max_val = CppAD::CondExpGt(logMixUrb(i,q), max_val, logMixUrb(i,q), max_val);
    }
    Type sum_exp = Type(0.0);
    for(int q = 0; q < Q; q++) {
      sum_exp += exp(logMixUrb(i, q) - max_val);
    }
    // log L_i = log(1/sqrt(pi)) + max_q logMix + log( sum_q exp(...) )
    jnll -= log_inv_sqrt_pi + max_val + log(sum_exp);
  }

  // ── Rural: same structure as Urban above — see comments there for
  //          the inline-binomial / logspace_add / logSumExp logic. ──
  matrix<Type> logMixRur(nRur, Q);

  for(int q = 0; q < Q; q++) {
    for(int i = 0; i < nRur; i++) {
      Type mix_lik = Type(0.0);
      Type yi = y_iRuralMICS(i);
      Type ni = n_iRuralMICS(i);

      for(int k = 0; k < KRur; k++) {
        Type w_ik = wRuralMICS(i, k);
        Type eta = base_eta_rur(i, k) + eps_gh(q);
        Type log_binom = yi * eta - ni * logspace_add(Type(0), eta);
        mix_lik += w_ik * exp(log_binom);
      }

      logMixRur(i, q) = log_gh_w(q) + lchoose_rural(i) + log(mix_lik);
    }
  }

  for(int i = 0; i < nRur; i++) {
    Type max_val = logMixRur(i, 0);
    for(int q = 1; q < Q; q++) {
      max_val = CppAD::CondExpGt(logMixRur(i,q), max_val, logMixRur(i,q), max_val);
    }
    Type sum_exp = Type(0.0);
    for(int q = 0; q < Q; q++) {
      sum_exp += exp(logMixRur(i, q) - max_val);
    }
    jnll -= log_inv_sqrt_pi + max_val + log(sum_exp);
  }

  // ── Reports ──
  if(options == 1) {
    ADREPORT(log_tau);
    ADREPORT(logit_phi);
  }

  REPORT(w_bym2);
  REPORT(u_bym2);
  REPORT(log_tau);
  REPORT(logit_phi);
  REPORT(log_tauEps);
  REPORT(alpha);
  REPORT(beta);
  REPORT(quadSum);
  REPORT(logDetTau);
  REPORT(bym2LogLik);
  REPORT(logPriTau);

  return jnll;
}
