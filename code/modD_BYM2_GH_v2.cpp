// ============================================================================
// modD_BYM2_GH_v2.cpp  —  DHS-only BYM2 with Gauss-Hermite nugget integration.
//
// This is EXACTLY modMDM_BYM2_GH_v2.cpp with the MICS data blocks removed: one
// shared BYM2 (w, u) plus one shared nugget variance, fit to DHS data only.
// The DHS likelihood block is byte-for-byte the DHS block of modMDM (same
// add_dataset lambda, same areaidx field lookup, same inline binomial + GH
// integration). See modMDM_BYM2_GH_v2.cpp / modM_BYM2_GH_v2.cpp for the
// canonical comments on:
//   * BYM2 (w, u) parameterisation, sum-to-zero, density derivation
//   * the inline binomial (y*eta - n*logspace_add(0, eta)) and why we don't use
//     dbinom_robust
//   * logspace_add(0, eta) = log(1 + exp(eta)), the stable log(1+exp) form
//   * the Gauss-Hermite change of variables and logSumExp over Q nodes
//
// IMPORTANT (the M_D bug fix, 2026-06-29): the DHS spatial effect is applied
// PER CLUSTER at its nominal stratum via areaidxlocUrban/RuralDHS (a DATA_IVECTOR
// of 0-based stratum indices, length nObs), NOT via a dense Aproj projection
// matrix. DHS jitter integration points never cross strata with nonzero weight
// (areasUrban/Rural are per-cluster), so the area lookup is exact; and unlike
// the old Aproj path it is correctly subsetted by subsetMDMInputs in held-out
// fits, so the field stays aligned with the data in cross-validation.
// ============================================================================

#include <TMB.hpp>
#include <Eigen/Sparse>
using namespace density;
using Eigen::SparseMatrix;

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
  // DHS data
  DATA_VECTOR( y_iUrbanDHS );
  DATA_VECTOR( y_iRuralDHS );
  DATA_VECTOR( n_iUrbanDHS );
  DATA_VECTOR( n_iRuralDHS );
  DATA_IVECTOR(areaidxlocUrbanDHS);
  DATA_IVECTOR(areaidxlocRuralDHS);
  DATA_MATRIX( X_betaUrbanDHS );
  DATA_MATRIX( X_betaRuralDHS );
  DATA_ARRAY( wUrbanDHS );
  DATA_ARRAY( wRuralDHS );

  DATA_SPARSE_MATRIX( Q_bym2 );

  // Precomputed log choose
  DATA_VECTOR( lchoose_urban_dhs );
  DATA_VECTOR( lchoose_rural_dhs );

  // GH data
  DATA_VECTOR( gh_nodes );
  DATA_VECTOR( gh_weights );

  // Priors
  DATA_VECTOR( alpha_pri );
  DATA_VECTOR( beta_pri );
  DATA_SCALAR( tr );
  DATA_VECTOR( gammaTildesm1 );
  DATA_SCALAR( lambdaPhi );
  DATA_INTEGER( uniformPhiPrior );
  DATA_SCALAR( lambdaTau );
  DATA_SCALAR( lambdaTauEps );
  DATA_SCALAR( options );

  // Outer
  PARAMETER( log_tau );
  PARAMETER( logit_phi );
  PARAMETER( log_tauEps );

  // Inner
  PARAMETER( alpha );
  PARAMETER_VECTOR( beta );
  PARAMETER_VECTOR( w_bym2Free );
  PARAMETER_VECTOR( u_bym2Free );

  int nUrbD = y_iUrbanDHS.size();
  int nRurD = y_iRuralDHS.size();
  int nFree = w_bym2Free.size();
  int nAreas = nFree + 1;
  int KUrbD = wUrbanDHS.cols();
  int KRurD = wRuralDHS.cols();
  int Q = gh_nodes.size();

  Type tau = exp(log_tau);
  Type phi = Type(1.0)/(Type(1.0) + exp(-logit_phi));
  Type sigmaEps = exp(Type(-0.5) * log_tauEps);

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

  // BYM2 density
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
  Type bym2LogLik = Type(-0.5) * logDetTau + Type(-0.5) * quadSum;
  jnll -= bym2LogLik;

  // Priors
  Type logPriTau = log(lambdaTau) - log(Type(2.0)) - Type(1.5) * log_tau
                   - lambdaTau / sqrt(tau) + log_tau;
  jnll -= logPriTau;
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

  vector<Type> fe_urb_d = X_betaUrbanDHS * beta;
  vector<Type> fe_rur_d = X_betaRuralDHS * beta;

  vector<Type> eps_gh(Q);
  vector<Type> log_gh_w(Q);
  for(int q = 0; q < Q; q++) {
    eps_gh(q) = sqrt(Type(2.0)) * sigmaEps * gh_nodes(q);
    log_gh_w(q) = log(gh_weights(q));
  }
  Type log_inv_sqrt_pi = -Type(0.5) * log(M_PI);

  // add_dataset: identical to modMDM_BYM2_GH_v2.cpp. Runs the GH-integrated
  // log-likelihood for one DHS data block (Urban or Rural):
  //   * base_eta_{i,k} = alpha + Xbeta + w_bym2[area_i]  (area_i constant over k)
  //   * for each (i, GH node q) sum the inline binomial over the K integration
  //     points with weights w(i, k), guarding against underflow
  //   * logSumExp over Q nodes to integrate out the cluster nugget eps_i
  auto add_dataset = [&](int nObs, int Kint,
                         const vector<Type>& y,
                         const vector<Type>& n,
                         const vector<int>& areaidx,
                         const vector<Type>& fe,
                         const matrix<Type>& w,
                         const vector<Type>& lchoose_vec)
  {
    matrix<Type> base_eta(nObs, Kint);
    for(int k = 0; k < Kint; k++) {
      for(int i = 0; i < nObs; i++) {
        base_eta(i, k) = alpha + fe(nObs * k + i) + w_bym2(areaidx(i));
      }
    }

    matrix<Type> logMix(nObs, Q);
    for(int q = 0; q < Q; q++) {
      for(int i = 0; i < nObs; i++) {
        Type mix_lik = Type(0.0);
        Type yi = y(i);
        Type ni = n(i);
        for(int k = 0; k < Kint; k++) {
          Type w_ik = w(i, k);
          Type eta = base_eta(i, k) + eps_gh(q);
          Type log_binom = yi * eta - ni * logspace_add(Type(0), eta);
          mix_lik += w_ik * exp(log_binom);
        }
        // Guard against mix_lik being zero or negative due to underflow/rounding
        Type tiny = Type(1e-200);
        if(mix_lik <= Type(0.0)) mix_lik = tiny;
        logMix(i, q) = log_gh_w(q) + lchoose_vec(i) + log(mix_lik);
      }
    }

    for(int i = 0; i < nObs; i++) {
      Type max_val = logMix(i, 0);
      for(int q = 1; q < Q; q++) {
        max_val = CppAD::CondExpGt(logMix(i,q), max_val, logMix(i,q), max_val);
      }
      Type sum_exp = Type(0.0);
      for(int q = 0; q < Q; q++) {
        sum_exp += exp(logMix(i, q) - max_val);
      }
      jnll -= log_inv_sqrt_pi + max_val + log(sum_exp);
    }
  };

  // Convert DATA_ARRAY weight objects into proper matrices (nObs x Kint) before passing
  matrix<Type> wUrbD(nUrbD, KUrbD);
  for(int i=0;i<nUrbD;i++) for(int k=0;k<KUrbD;k++) wUrbD(i,k) = wUrbanDHS(i,k);
  matrix<Type> wRurD(nRurD, KRurD);
  for(int i=0;i<nRurD;i++) for(int k=0;k<KRurD;k++) wRurD(i,k) = wRuralDHS(i,k);

  add_dataset(nUrbD, KUrbD, y_iUrbanDHS, n_iUrbanDHS, areaidxlocUrbanDHS,
              fe_urb_d, wUrbD, lchoose_urban_dhs);
  add_dataset(nRurD, KRurD, y_iRuralDHS, n_iRuralDHS, areaidxlocRuralDHS,
              fe_rur_d, wRurD, lchoose_rural_dhs);

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
