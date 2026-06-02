// ============================================================================
// modMDM_BYM2_GH_v2.cpp  —  MICS + DHS data fusion with a shared BYM2 spatial
// effect and Gauss-Hermite nugget integration.
//
// This is the joint version of the two single-source templates: it takes both
// MICS and DHS data (with their respective integration-point schemes) and
// fits one shared BYM2 (w, u) plus one shared nugget variance to both. See
// modM_BYM2_GH_v2.cpp for the canonical comments on:
//   * BYM2 (w, u) parameterisation, sum-to-zero, density derivation
//   * inline binomial (y*eta - n*logspace_add(0, eta)) and why we don't use
//     dbinom_robust
//   * logspace_add(0, eta) = log(1 + exp(eta)), the stable log(1+exp) form
//   * Gauss-Hermite change of variables, logSumExp over Q nodes
//
// What's different here:
//   * Both data sources contribute additive likelihood terms — MICS uses
//     areaidxloc lookup (effectively K_MICS = 100 integration points drawn
//     from the stratum × urban polygon), DHS uses its smaller K integration
//     points (~16-21) from the jittering disc.
//   * One BYM2 effect, one nugget variance, one alpha+beta shared across
//     both sources. The shared alpha is the cross-source identifiability
//     anchor; if you need source-specific intercepts, see modMDM_BYM2_GH_comm.
//   * Computational footprint is ~2× the single-source templates because
//     both inner data loops run, but the BYM2 GMRF density is computed once.
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
  // MICS data
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
  DATA_VECTOR( lchoose_urban_mics );
  DATA_VECTOR( lchoose_rural_mics );
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

  int nUrbM = y_iUrbanMICS.size();
  int nRurM = y_iRuralMICS.size();
  int nUrbD = y_iUrbanDHS.size();
  int nRurD = y_iRuralDHS.size();
  int nFree = w_bym2Free.size();
  int nAreas = nFree + 1;
  int KUrbM = wUrbanMICS.cols();
  int KRurM = wRuralMICS.cols();
  int KUrbD = wUrbanDHS.cols();
  int KRurD = wRuralDHS.cols();
  int Q = gh_nodes.size();

    // Quick runtime debug prints to help localize Eigen bounds errors during MakeADFun
    Rprintf("[modMDM_BYM2_GH_v2] sizes: nUrbM=%d nRurM=%d nUrbD=%d nRurD=%d nAreas=%d KUrbM=%d KRurM=%d KUrbD=%d KRurD=%d Q=%d\n",
      nUrbM, nRurM, nUrbD, nRurD, nAreas, KUrbM, KRurM, KUrbD, KRurD, Q);

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
  Rprintf("[modMDM_BYM2_GH_v2] about to multiply Q_bym2 by u_bym2 (nAreas=%d)\n", nAreas);
  vector<Type> Qu = (Q_bym2 * u_bym2.matrix()).col(0);
  Rprintf("[modMDM_BYM2_GH_v2] completed Q_bym2 * u_bym2 multiply\n");
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

  vector<Type> fe_urb_m = X_betaUrbanMICS * beta;
  vector<Type> fe_rur_m = X_betaRuralMICS * beta;
  vector<Type> fe_urb_d = X_betaUrbanDHS * beta;
  vector<Type> fe_rur_d = X_betaRuralDHS * beta;

  vector<Type> eps_gh(Q);
  vector<Type> log_gh_w(Q);
  for(int q = 0; q < Q; q++) {
    eps_gh(q) = sqrt(Type(2.0)) * sigmaEps * gh_nodes(q);
    log_gh_w(q) = log(gh_weights(q));
  }
  Type log_inv_sqrt_pi = -Type(0.5) * log(M_PI);

  // add_dataset: a generic lambda that runs the GH-integrated log-likelihood
  // for one data block (MICS Urban, MICS Rural, DHS Urban, or DHS Rural).
  // Each block:
  //   * builds base_eta_{i,k} = alpha + Xbeta + w_bym2[area_i] once
  //   * for each (i, GH node q) sums the inline binomial likelihood over the
  //     K integration points with weights w(i, k), guarding against underflow
  //   * logSumExp over Q nodes to integrate out the cluster nugget eps_i
  // See the canonical comments in modM_BYM2_GH_v2.cpp for full notes on the
  // inline binomial (y*eta - n*logspace_add(0, eta)) and why we don't use
  // dbinom_robust.
  auto add_dataset = [&](int nObs, int Kint,
                         const vector<Type>& y,
                         const vector<Type>& n,
                         const vector<int>& areaidx,
                         const vector<Type>& fe,
                         const matrix<Type>& w,
                         const vector<Type>& lchoose_vec)
  {
    int area_min = 0;
    int area_max = 0;
    if(areaidx.size() > 0) {
      area_min = areaidx(0);
      area_max = areaidx(0);
      for(int ii = 0; ii < areaidx.size(); ii++){
        if(areaidx(ii) < area_min) area_min = areaidx(ii);
        if(areaidx(ii) > area_max) area_max = areaidx(ii);
      }
    }
    Rprintf("[modMDM_BYM2_GH_v2] add_dataset start: nObs=%d Kint=%d fe.size=%d w.rows=%d w.cols=%d lchoose=%d areaidx.len=%d areaidx.range=%d..%d\n",
            nObs, Kint, (int)fe.size(), (int)w.rows(), (int)w.cols(), (int)lchoose_vec.size(), (int)areaidx.size(), area_min, area_max);

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
  matrix<Type> wUrbM(nUrbM, KUrbM);
  for(int i=0;i<nUrbM;i++) for(int k=0;k<KUrbM;k++) wUrbM(i,k) = wUrbanMICS(i,k);
  matrix<Type> wRurM(nRurM, KRurM);
  for(int i=0;i<nRurM;i++) for(int k=0;k<KRurM;k++) wRurM(i,k) = wRuralMICS(i,k);
  matrix<Type> wUrbD(nUrbD, KUrbD);
  for(int i=0;i<nUrbD;i++) for(int k=0;k<KUrbD;k++) wUrbD(i,k) = wUrbanDHS(i,k);
  matrix<Type> wRurD(nRurD, KRurD);
  for(int i=0;i<nRurD;i++) for(int k=0;k<KRurD;k++) wRurD(i,k) = wRuralDHS(i,k);

  add_dataset(nUrbM, KUrbM, y_iUrbanMICS, n_iUrbanMICS, areaidxlocUrbanMICS,
              fe_urb_m, wUrbM, lchoose_urban_mics);
  add_dataset(nRurM, KRurM, y_iRuralMICS, n_iRuralMICS, areaidxlocRuralMICS,
              fe_rur_m, wRurM, lchoose_rural_mics);
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
