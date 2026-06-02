// ============================================================================
// modM_FEnug_GH.cpp  —  Fixed-effects + nugget model with Gauss-Hermite
// integration of the nugget. NO BYM2 spatial structure (an IID spatial effect
// `u_spatial` can be included by leaving `log_tau` un-mapped, or removed by
// mapping it out via `map = list(log_tau = factor(NA), u_spatial = ...)`).
//
// MODEL
// -----
//   y_i | n_i, p_i        ~ Binomial(n_i, p_i)
//   logit p_i              = alpha + x_i^T beta [+ u_{a(i)}] + eps_i
//   u_a                    ~ iid N(0, 1/tau)        (IID area effect; optional)
//   eps_i                  ~ N(0, 1/tau_eps)        (cluster-level nugget)
//
// Like the BYM2 GH templates, the cluster nugget eps_i is integrated out by
// Q-point Gauss-Hermite quadrature so the inner random effects passed to
// Laplace are (alpha, beta [, u_spatial]) only.
//
// SHARED OPTIMISATIONS — see modM_BYM2_GH_v2.cpp's header for full notes:
//   * INLINE binomial likelihood (not dbinom_robust):
//       log_binom_iqk = lchoose(n_i,y_i) + y_i*eta - n_i*logspace_add(0, eta)
//     We add lchoose inside this loop (rather than factoring out as DATA) for
//     clarity here — it's also constant w.r.t. params and gets the same AD
//     treatment either way.
//   * logspace_add(0, eta) = log(1 + exp(eta)) — stable log-of-1-plus-exp,
//     binomial log-normaliser. Computed as max(0, eta) + log1p(exp(-|eta|))
//     to avoid overflow at large positive eta.
//   * Gauss-Hermite change of variables eps = sqrt(2)*sigma_eps*z gives the
//     standard nodes/weights form; (1/sqrt(pi)) is the Jacobian.
//   * logSumExp over Q nodes for the final per-obs log-likelihood.
//
// Differences from the BYM2 templates:
//   * No BYM2 (w, u) decomposition, no sum-to-zero quadratic form.
//   * Spatial effect (if present) is iid, prior dnorm(u, 0, sigma_u).
//   * lchoose is added inside the K-loop rather than factored to DATA — same
//     AD-tape footprint either way since lchoose is constant.
// ============================================================================

#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
using namespace density;
using Eigen::SparseMatrix;

// PC prior on precision parameter (on log scale)
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
  DATA_MATRIX( X_betaUrbanMICS );   // (nUrb * K) x nBeta
  DATA_MATRIX( X_betaRuralMICS );   // (nRur * K) x nBeta
  DATA_ARRAY( wUrbanMICS );         // nUrb x K
  DATA_ARRAY( wRuralMICS );         // nRur x K

  DATA_INTEGER( nAreas );

  // Precomputed log-binomial-coefficient (DATA — not on AD tape)
  DATA_VECTOR( lchoose_urban );      // length nUrb
  DATA_VECTOR( lchoose_rural );      // length nRur

  // Gauss-Hermite quadrature nodes and weights
  DATA_VECTOR( gh_nodes );           // Q nodes
  DATA_VECTOR( gh_weights );         // Q weights

  // Prior parameters
  DATA_VECTOR( alpha_pri );
  DATA_VECTOR( beta_pri );
  DATA_SCALAR( lambdaTau );
  DATA_SCALAR( lambdaTauEps );

  // ── Parameters ──
  PARAMETER( log_tau );       // IID spatial precision (map out for FE-only)
  PARAMETER( log_tauEps );    // nugget precision

  PARAMETER( alpha );         // intercept
  PARAMETER_VECTOR( beta );   // covariate coefficients (generic length)

  PARAMETER_VECTOR( u_spatial ); // IID spatial effects (map out for FE-only)

  // ── Dimensions ──
  int nUrb = y_iUrbanMICS.size();
  int nRur = y_iRuralMICS.size();
  int KUrb = wUrbanMICS.cols();
  int KRur = wRuralMICS.cols();
  int Q = gh_nodes.size();

  // ── Transforms ──
  Type sigma_u = exp(Type(-0.5) * log_tau);
  Type sigmaEps = exp(Type(-0.5) * log_tauEps);

  Type jnll = Type(0.0);

  // ═══════════════════════════════════
  // (1) Priors
  // ═══════════════════════════════════

  // IID spatial random effects prior
  for(int i = 0; i < u_spatial.size(); i++) {
    jnll -= dnorm(u_spatial(i), Type(0), sigma_u, true);
  }

  // PC priors for precision parameters
  jnll -= dPCPriTau(log_tau, lambdaTau);
  jnll -= dPCPriTau(log_tauEps, lambdaTauEps);

  // Fixed effect priors
  jnll -= dnorm(alpha, alpha_pri(0), alpha_pri(1), true);
  for(int i = 0; i < beta.size(); i++) {
    jnll -= dnorm(beta(i), beta_pri(0), beta_pri(1), true);
  }

  // ═══════════════════════════════════
  // (2) Precompute base linear predictor
  // ═══════════════════════════════════
  // fe vectors: length (nObs * K), stacked by integration point
  vector<Type> feUrb = X_betaUrbanMICS * beta;
  vector<Type> feRur = X_betaRuralMICS * beta;

  // base_eta for each obs×intpt: alpha + Xbeta + u_spatial[area]
  // Stored as nObs x K matrix for efficient access
  matrix<Type> base_eta_urb(nUrb, KUrb);
  for(int k = 0; k < KUrb; k++) {
    for(int i = 0; i < nUrb; i++) {
      int idx = nUrb * k + i;  // column-major: obs i, integration point k
      base_eta_urb(i, k) = alpha + feUrb(idx) + u_spatial(areaidxlocUrbanMICS(i));
    }
  }

  matrix<Type> base_eta_rur(nRur, KRur);
  for(int k = 0; k < KRur; k++) {
    for(int i = 0; i < nRur; i++) {
      int idx = nRur * k + i;
      base_eta_rur(i, k) = alpha + feRur(idx) + u_spatial(areaidxlocRuralMICS(i));
    }
  }

  // GH abscissae: eps_q = sqrt(2) * sigmaEps * node_q
  vector<Type> eps_gh(Q);
  for(int q = 0; q < Q; q++) {
    eps_gh(q) = sqrt(Type(2.0)) * sigmaEps * gh_nodes(q);
  }

  // Constant: log(1/sqrt(pi))
  Type log_inv_sqrt_pi = -Type(0.5) * log(M_PI);

  // ═══════════════════════════════════
  // (3) Likelihood — Urban MICS
  // ═══════════════════════════════════
  Type tiny = Type(1e-200);
  for(int i = 0; i < nUrb; i++) {
    Type y_i = y_iUrbanMICS(i);
    Type n_i = n_iUrbanMICS(i);

    // For each GH node, compute the integration-point mixture
    // Then logSumExp over GH nodes
    Type max_log_term = Type(-1e10);
    vector<Type> log_terms(Q);

    for(int q = 0; q < Q; q++) {
      // Sum over integration points k for this GH node, weighted by w_ik.
      Type mix_lik = Type(0.0);
      for(int k = 0; k < KUrb; k++) {
        Type w_ik = wUrbanMICS(i, k);
        if(w_ik > Type(0.0)) {
          Type eta = base_eta_urb(i, k) + eps_gh(q);   // logit(p_iqk)
          // INLINE binomial log-likelihood:
          //   log_binom = lchoose(n,y) + y*eta - n*log(1+exp(eta))
          // logspace_add(0, eta) = log(1+exp(eta)) — stable for any eta
          // (avoids overflow at large positive eta). See header notes for
          // why we don't call dbinom_robust here.
          Type log_lik_k = lchoose_urban(i) + y_i * eta - n_i * logspace_add(Type(0.0), eta);
          mix_lik += w_ik * exp(log_lik_k);
        }
      }
      // Underflow guard: if all log_lik_k are very negative, mix_lik can
      // round to 0 and log(0) = -Inf would poison the AD tape.
      mix_lik = CppAD::CondExpGt(mix_lik, tiny, mix_lik, tiny);
      log_terms(q) = log(gh_weights(q)) + log(mix_lik);
      if(q == 0 || asDouble(log_terms(q)) > asDouble(max_log_term))
        max_log_term = log_terms(q);
    }

    // logSumExp
    Type sum_exp = Type(0.0);
    for(int q = 0; q < Q; q++) {
      sum_exp += exp(log_terms(q) - max_log_term);
    }
    Type log_lik_i = log_inv_sqrt_pi + max_log_term + log(sum_exp);
    jnll -= log_lik_i;
  }

  // ═══════════════════════════════════
  // (4) Likelihood — Rural MICS
  // ═══════════════════════════════════
  for(int i = 0; i < nRur; i++) {
    Type y_i = y_iRuralMICS(i);
    Type n_i = n_iRuralMICS(i);

    Type max_log_term = Type(-1e10);
    vector<Type> log_terms(Q);

    for(int q = 0; q < Q; q++) {
      Type mix_lik = Type(0.0);
      for(int k = 0; k < KRur; k++) {
        Type w_ik = wRuralMICS(i, k);
        if(w_ik > Type(0.0)) {
          Type eta = base_eta_rur(i, k) + eps_gh(q);
          Type log_lik_k = lchoose_rural(i) + y_i * eta - n_i * logspace_add(Type(0.0), eta);
          mix_lik += w_ik * exp(log_lik_k);
        }
      }
      mix_lik = CppAD::CondExpGt(mix_lik, tiny, mix_lik, tiny);
      log_terms(q) = log(gh_weights(q)) + log(mix_lik);
      if(q == 0 || asDouble(log_terms(q)) > asDouble(max_log_term))
        max_log_term = log_terms(q);
    }

    Type sum_exp = Type(0.0);
    for(int q = 0; q < Q; q++) {
      sum_exp += exp(log_terms(q) - max_log_term);
    }
    Type log_lik_i = log_inv_sqrt_pi + max_log_term + log(sum_exp);
    jnll -= log_lik_i;
  }

  // ── Report ──
  REPORT(log_tau);
  REPORT(log_tauEps);
  REPORT(alpha);
  REPORT(beta);
  REPORT(u_spatial);
  REPORT(sigma_u);
  REPORT(sigmaEps);

  return jnll;
}
