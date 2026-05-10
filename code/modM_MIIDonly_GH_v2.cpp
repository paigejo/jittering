// IID spatial + nugget model — nuggets integrated out via Gauss-Hermite quadrature
// OPTIMIZED VERSION: matrix operations replace triple loop
// - base_eta computed once (not per GH node)
// - vectorized binomial log-lik (no dbinom_robust overhead)
// - no if-branching on AD tape
// - nObs x Q matrix multiplied by GH weights

#include <TMB.hpp>
#include <Eigen/Sparse>
using namespace density;

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
  // Data
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

  // Precomputed log-binomial-coefficient (DATA, not on AD tape)
  DATA_VECTOR( lchoose_urban );  // lchoose(n, y) for each urban obs
  DATA_VECTOR( lchoose_rural );  // lchoose(n, y) for each rural obs

  // Gauss-Hermite quadrature
  DATA_VECTOR( gh_nodes );
  DATA_VECTOR( gh_weights );

  // Prior parameters
  DATA_VECTOR( alpha_pri );
  DATA_VECTOR( beta_pri );
  DATA_SCALAR( lambdaTau );
  DATA_SCALAR( lambdaTauEps );

  // Parameters
  PARAMETER( log_tau );
  PARAMETER( log_tauEps );
  PARAMETER( alpha );
  PARAMETER( beta_urban );
  PARAMETER( beta_normPop );
  PARAMETER_VECTOR( u_spatial );

  // Dimensions
  int nUrb = y_iUrbanMICS.size();
  int nRur = y_iRuralMICS.size();
  int KUrb = wUrbanMICS.cols();
  int KRur = wRuralMICS.cols();
  int Q = gh_nodes.size();

  // Transforms
  Type sigma_u = exp(Type(-0.5) * log_tau);
  Type sigmaEps = exp(Type(-0.5) * log_tauEps);

  Type jnll = 0.0;

  // ═══ (1) Priors ═══
  for(int i = 0; i < u_spatial.size(); i++) {
    jnll -= dnorm(u_spatial(i), Type(0), sigma_u, true);
  }
  jnll -= dPCPriTau(log_tau, lambdaTau);
  jnll -= dPCPriTau(log_tauEps, lambdaTauEps);
  jnll -= dnorm(alpha, alpha_pri(0), alpha_pri(1), true);
  jnll -= dnorm(beta_urban, beta_pri(0), beta_pri(1), true);
  jnll -= dnorm(beta_normPop, beta_pri(0), beta_pri(1), true);

  // ═══ (2) Precompute fixed effects ═══
  vector<Type> fe_urb = X_betaUrbanMICS.col(0) * beta_urban +
                        X_betaUrbanMICS.col(1) * beta_normPop;
  vector<Type> fe_rur = X_betaRuralMICS.col(0) * beta_urban +
                        X_betaRuralMICS.col(1) * beta_normPop;

  // ═══ (3) Build base_eta matrices: nObs x K ═══
  // base_eta(i,k) = alpha + fe(nObs*k + i) + u_spatial(area(i))
  // Computed ONCE, not Q times

  matrix<Type> base_eta_urb(nUrb, KUrb);
  for(int k = 0; k < KUrb; k++) {
    for(int i = 0; i < nUrb; i++) {
      base_eta_urb(i, k) = alpha + fe_urb(nUrb * k + i)
                            + u_spatial(areaidxlocUrbanMICS(i));
    }
  }

  matrix<Type> base_eta_rur(nRur, KRur);
  for(int k = 0; k < KRur; k++) {
    for(int i = 0; i < nRur; i++) {
      base_eta_rur(i, k) = alpha + fe_rur(nRur * k + i)
                            + u_spatial(areaidxlocRuralMICS(i));
    }
  }

  // ═══ (4) GH quadrature — vectorized ═══
  // For each GH node q, compute weighted binomial lik for all obs
  // Then combine across Q nodes via logSumExp

  // GH abscissae
  vector<Type> eps_gh(Q);
  for(int q = 0; q < Q; q++) {
    eps_gh(q) = sqrt(Type(2.0)) * sigmaEps * gh_nodes(q);
  }

  Type log_inv_sqrt_pi = -Type(0.5) * log(M_PI);

  // Log GH weights (precompute once)
  vector<Type> log_gh_w(Q);
  for(int q = 0; q < Q; q++) {
    log_gh_w(q) = log(gh_weights(q));
  }

  // ── Urban: build nUrb x Q matrix of log(weighted mixture lik per GH node) ──
  matrix<Type> logMixUrb(nUrb, Q);

  for(int q = 0; q < Q; q++) {
    // For this GH node, compute integration-point-weighted binomial lik for each obs
    // mix_lik(i) = sum_k w(i,k) * Binom(y_i | n_i, eta(i,k) + eps_q)
    // Using: log Binom = lchoose + y*eta - n*log(1+exp(eta))  [robust form]
    // But we need: sum_k w_k * exp(log_binom_k) — must stay in probability space for weighted sum

    for(int i = 0; i < nUrb; i++) {
      Type mix_lik = Type(0.0);
      Type yi = y_iUrbanMICS(i);
      Type ni = n_iUrbanMICS(i);

      for(int k = 0; k < KUrb; k++) {
        Type w_ik = wUrbanMICS(i, k);
        // No branching: w_ik * 0 = 0 when w=0 and eta is garbage
        // But we need valid eta even when w=0 to avoid NaN on tape
        Type eta = base_eta_urb(i, k) + eps_gh(q);
        // Inline binomial (robust): p = 1/(1+exp(-eta))
        // dbinom = exp(lchoose) * p^y * (1-p)^(n-y)
        //        = exp(lchoose + y*eta - n*log(1+exp(eta)))
        Type log_binom = yi * eta - ni * logspace_add(Type(0), eta);
        mix_lik += w_ik * exp(log_binom);  // stays in prob space for weighted sum
      }

      // lchoose is constant data — add it after the weighted sum
      // mix_lik already has the exp(lchoose) factored out:
      // true_mix_lik = exp(lchoose) * mix_lik
      // log(true_mix_lik) = lchoose + log(mix_lik)
      logMixUrb(i, q) = log_gh_w(q) + lchoose_urban(i) + log(mix_lik);
    }
  }

  // Sum over GH nodes via logSumExp per obs
  for(int i = 0; i < nUrb; i++) {
    // logSumExp over Q columns for row i
    Type max_val = logMixUrb(i, 0);
    for(int q = 1; q < Q; q++) {
      max_val = CppAD::CondExpGt(logMixUrb(i,q), max_val, logMixUrb(i,q), max_val);
    }
    Type sum_exp = Type(0.0);
    for(int q = 0; q < Q; q++) {
      sum_exp += exp(logMixUrb(i, q) - max_val);
    }
    jnll -= log_inv_sqrt_pi + max_val + log(sum_exp);
  }

  // ── Rural: same structure ──
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

  // Report
  REPORT(log_tau);
  REPORT(log_tauEps);
  REPORT(alpha);
  REPORT(beta_urban);
  REPORT(beta_normPop);
  REPORT(u_spatial);
  REPORT(sigma_u);
  REPORT(sigmaEps);

  return jnll;
}
