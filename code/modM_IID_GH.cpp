// IID spatial + GH nuggets for MICS-only
// Same likelihood structure as modM_BYM2_GH_v2.cpp but with
// IID N(0, 1/tau) spatial effects replacing BYM2

#include <TMB.hpp>

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

  // Precomputed log-binomial-coefficient
  DATA_VECTOR( lchoose_urban );
  DATA_VECTOR( lchoose_rural );

  // GH quadrature nodes/weights
  DATA_VECTOR( gh_nodes );
  DATA_VECTOR( gh_weights );

  // Prior parameters
  DATA_VECTOR( alpha_pri );
  DATA_VECTOR( beta_pri );

  DATA_SCALAR( lambdaTau );
  DATA_SCALAR( lambdaTauEps );

  DATA_SCALAR( options );

  // ── Parameters ──
  PARAMETER( log_tau );
  PARAMETER( log_tauEps );

  PARAMETER( alpha );
  PARAMETER_VECTOR( beta );

  PARAMETER_VECTOR( v_iidFree );  // n-1, sum-to-zero constrained

  // ── Dimensions ──
  int nUrb = y_iUrbanMICS.size();
  int nRur = y_iRuralMICS.size();
  int nFree = v_iidFree.size();
  int nAreas = nFree + 1;
  int KUrb = wUrbanMICS.cols();
  int KRur = wRuralMICS.cols();
  int Q = gh_nodes.size();

  // ── Transforms ──
  Type tau = exp(log_tau);
  Type sigmaEps = exp(Type(-0.5) * log_tauEps);

  // Reconstruct full n-vector with sum-to-zero
  vector<Type> v_iid(nAreas);
  Type vSum = Type(0.0);
  for(int i = 0; i < nFree; i++) {
    v_iid(i) = v_iidFree(i);
    vSum += v_iidFree(i);
  }
  v_iid(nFree) = -vSum;

  Type jnll = Type(0.0);

  // ═══════════════════════════════════
  // (1) IID spatial density with sum-to-zero constraint
  //     Constrained (n-1)-dim Gaussian: 0.5*(n-1)*log(tau) - 0.5*tau*v^Tv
  // ═══════════════════════════════════
  Type quadV = (v_iid * v_iid).sum() * tau;
  Type iidLogLik = Type(0.5) * Type(nFree) * log_tau - Type(0.5) * quadV;
  jnll -= iidLogLik;

  // ═══════════════════════════════════
  // (2) Priors
  // ═══════════════════════════════════
  jnll -= dPCPriTau(log_tau, lambdaTau);
  jnll -= dPCPriTau(log_tauEps, lambdaTauEps);

  jnll -= dnorm(alpha, alpha_pri(0), alpha_pri(1), true);
  for(int i = 0; i < beta.size(); i++) {
    jnll -= dnorm(beta(i), beta_pri(0), beta_pri(1), true);
  }

  // ═══════════════════════════════════
  // (3) GH-integrated likelihood
  // ═══════════════════════════════════
  vector<Type> fe_urb = X_betaUrbanMICS * beta;
  vector<Type> fe_rur = X_betaRuralMICS * beta;

  vector<Type> eps_gh(Q);
  vector<Type> log_gh_w(Q);
  for(int q = 0; q < Q; q++) {
    eps_gh(q) = sqrt(Type(2.0)) * sigmaEps * gh_nodes(q);
    log_gh_w(q) = log(gh_weights(q));
  }
  Type log_inv_sqrt_pi = -Type(0.5) * log(M_PI);

  // ── Urban ──
  matrix<Type> base_eta_urb(nUrb, KUrb);
  for(int k = 0; k < KUrb; k++) {
    for(int i = 0; i < nUrb; i++) {
      base_eta_urb(i, k) = alpha + fe_urb(nUrb * k + i)
                            + v_iid(areaidxlocUrbanMICS(i));
    }
  }

  matrix<Type> logMixUrb(nUrb, Q);
  for(int q = 0; q < Q; q++) {
    for(int i = 0; i < nUrb; i++) {
      Type mix_lik = Type(0.0);
      Type yi = y_iUrbanMICS(i);
      Type ni = n_iUrbanMICS(i);
      for(int k = 0; k < KUrb; k++) {
        Type w_ik = wUrbanMICS(i, k);
        Type eta = base_eta_urb(i, k) + eps_gh(q);
        Type log_binom = yi * eta - ni * logspace_add(Type(0), eta);
        mix_lik += w_ik * exp(log_binom);
      }
      logMixUrb(i, q) = log_gh_w(q) + lchoose_urban(i) + log(mix_lik);
    }
  }

  for(int i = 0; i < nUrb; i++) {
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

  // ── Rural ──
  matrix<Type> base_eta_rur(nRur, KRur);
  for(int k = 0; k < KRur; k++) {
    for(int i = 0; i < nRur; i++) {
      base_eta_rur(i, k) = alpha + fe_rur(nRur * k + i)
                            + v_iid(areaidxlocRuralMICS(i));
    }
  }

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
  REPORT(v_iid);
  REPORT(log_tau);
  REPORT(log_tauEps);
  REPORT(alpha);
  REPORT(beta);

  return jnll;
}
