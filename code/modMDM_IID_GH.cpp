// IID spatial + GH nuggets for MDM (MICS + DHS fusion)
// Same likelihood structure as modMDM_BYM2_GH_v2.cpp but with
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
  DATA_SCALAR( lambdaTau );
  DATA_SCALAR( lambdaTauEps );
  DATA_SCALAR( options );

  // Parameters
  PARAMETER( log_tau );
  PARAMETER( log_tauEps );

  PARAMETER( alpha );
  PARAMETER_VECTOR( beta );
  PARAMETER_VECTOR( v_iidFree );  // n-1, sum-to-zero constrained

  int nUrbM = y_iUrbanMICS.size();
  int nRurM = y_iRuralMICS.size();
  int nUrbD = y_iUrbanDHS.size();
  int nRurD = y_iRuralDHS.size();
  int nFree = v_iidFree.size();
  int nAreas = nFree + 1;
  int KUrbM = wUrbanMICS.cols();
  int KRurM = wRuralMICS.cols();
  int KUrbD = wUrbanDHS.cols();
  int KRurD = wRuralDHS.cols();
  int Q = gh_nodes.size();

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
        base_eta(i, k) = alpha + fe(nObs * k + i) + v_iid(areaidx(i));
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

  // Convert DATA_ARRAY weight objects into proper matrices
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

  // ── Reports ──
  REPORT(v_iid);
  REPORT(log_tau);
  REPORT(log_tauEps);
  REPORT(alpha);
  REPORT(beta);

  return jnll;
}
