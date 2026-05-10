// DHS-only BYM2 constrained (2n-2) + GH nuggets — OPTIMIZED v2
// Same structure as modM_BYM2_GH_v2.cpp but for DHS data
// DHS uses AprojUrbanDHS/AprojRuralDHS (nObs x nArea) to project spatial effects
// DHS integration points: 11 urban, 16 rural (fewer than MICS 100)

#include <TMB.hpp>
#include <Eigen/Sparse>
using namespace density;
using Eigen::SparseMatrix;

// PC prior on precision (log scale)
template<class Type>
Type dPCPriTau(Type logTau, Type lambda)
{
  Type tau = exp(logTau);
  Type ldensity = log(lambda) - log(Type(2.0)) - Type(1.5)*logTau - lambda/sqrt(tau);
  Type ljacobian = logTau;
  return ldensity + ljacobian;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // ── Data ──
  DATA_VECTOR( y_iUrbanDHS );
  DATA_VECTOR( y_iRuralDHS );
  DATA_VECTOR( n_iUrbanDHS );
  DATA_VECTOR( n_iRuralDHS );
  DATA_MATRIX( AprojUrbanDHS );  // nObsUrban x nArea
  DATA_MATRIX( AprojRuralDHS );  // nObsRural x nArea
  DATA_MATRIX( X_betaUrbanDHS ); // (nObsUrban * K) x nBeta
  DATA_MATRIX( X_betaRuralDHS );
  DATA_ARRAY( wUrbanDHS );       // nObsUrban x K
  DATA_ARRAY( wRuralDHS );

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
  PARAMETER( log_tau );
  PARAMETER( logit_phi );
  PARAMETER( log_tauEps );

  PARAMETER( alpha );
  PARAMETER_VECTOR( beta );

  PARAMETER_VECTOR( w_bym2Free );  // n-1
  PARAMETER_VECTOR( u_bym2Free );  // n-1

  // ── Dimensions ──
  int nUrb = y_iUrbanDHS.size();
  int nRur = y_iRuralDHS.size();
  int nFree = w_bym2Free.size();
  int nAreas = nFree + 1;
  int KUrb = wUrbanDHS.cols();
  int KRur = wRuralDHS.cols();
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
  Type quadW = (w_bym2 * w_bym2).sum() * tau / (Type(1.0) - phi);

  vector<Type> Qu = (Q_bym2 * u_bym2.matrix()).col(0);
  Type fac = phi / (Type(1.0) - phi);
  Type quadU = (Qu * u_bym2).sum() + fac * (u_bym2 * u_bym2).sum();

  Type quadWU = (u_bym2 * w_bym2).sum() * (-Type(2.0) * sqrt(phi * tau) / (Type(1.0) - phi));

  Type quadSum = quadW + quadU + quadWU;

  Type logDet = Type(0.0);
  for(int i = 0; i < gammaTildesm1.size(); i++) {
    logDet += log(Type(1.0) + phi * gammaTildesm1(i));
  }
  Type logDetTau = Type(nFree) * log((Type(1.0) - phi) / tau);

  jnll -= Type(-0.5) * logDetTau + Type(-0.5) * quadSum;

  // ═══════════════════════════════════
  // (2) Priors
  // ═══════════════════════════════════
  jnll -= log(lambdaTau) - log(Type(2.0)) - Type(1.5) * log_tau
          - lambdaTau / sqrt(tau) + log_tau;

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

  jnll -= dnorm(alpha, alpha_pri(0), alpha_pri(1), true);
  for(int i = 0; i < beta.size(); i++) {
    jnll -= dnorm(beta(i), beta_pri(0), beta_pri(1), true);
  }

  // ═══════════════════════════════════
  // (3) GH-integrated likelihood (vectorized)
  // ═══════════════════════════════════

  // Precompute fixed effects
  vector<Type> fe_urb = (X_betaUrbanDHS * beta.matrix()).col(0);
  vector<Type> fe_rur = (X_betaRuralDHS * beta.matrix()).col(0);

  // Project spatial effects to obs via Aproj (nObs x nArea) * w_bym2 (nArea)
  vector<Type> proj_urb = (AprojUrbanDHS * w_bym2.matrix()).col(0);
  vector<Type> proj_rur = (AprojRuralDHS * w_bym2.matrix()).col(0);

  // GH abscissae
  vector<Type> eps_gh(Q);
  for(int q = 0; q < Q; q++) {
    eps_gh(q) = sqrt(Type(2.0)) * sigmaEps * gh_nodes(q);
  }
  Type log_inv_sqrt_pi = -Type(0.5) * log(M_PI);

  vector<Type> log_gh_w(Q);
  for(int q = 0; q < Q; q++) {
    log_gh_w(q) = log(gh_weights(q));
  }

  // ── Build base_eta: computed ONCE, not Q times ──
  // DHS: fe is (nObs*K) long, indexed as fe(nObs*k + i) for obs i, intpt k
  // proj is length nObs

  matrix<Type> base_eta_urb(nUrb, KUrb);
  for(int k = 0; k < KUrb; k++) {
    for(int i = 0; i < nUrb; i++) {
      base_eta_urb(i, k) = alpha + fe_urb(nUrb * k + i) + proj_urb(i);
    }
  }

  matrix<Type> base_eta_rur(nRur, KRur);
  for(int k = 0; k < KRur; k++) {
    for(int i = 0; i < nRur; i++) {
      base_eta_rur(i, k) = alpha + fe_rur(nRur * k + i) + proj_rur(i);
    }
  }

  // ── Urban: nUrb x Q matrix of log(GH-weighted mixture lik) ──
  matrix<Type> logMixUrb(nUrb, Q);

  for(int q = 0; q < Q; q++) {
    for(int i = 0; i < nUrb; i++) {
      Type mix_lik = Type(0.0);
      Type yi = y_iUrbanDHS(i);
      Type ni = n_iUrbanDHS(i);

      for(int k = 0; k < KUrb; k++) {
        Type w_ik = wUrbanDHS(i, k);
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
  matrix<Type> logMixRur(nRur, Q);

  for(int q = 0; q < Q; q++) {
    for(int i = 0; i < nRur; i++) {
      Type mix_lik = Type(0.0);
      Type yi = y_iRuralDHS(i);
      Type ni = n_iRuralDHS(i);

      for(int k = 0; k < KRur; k++) {
        Type w_ik = wRuralDHS(i, k);
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

  return jnll;
}
