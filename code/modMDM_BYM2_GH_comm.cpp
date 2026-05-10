// BYM2 constrained (2n-2) + GH nuggets for MDM (MICS + DHS) — COMMENSURATE PRIORS
//
// Separate alpha/beta for MICS and DHS:
//   - alpha and alpha_M have independent N(alpha_pri[0], alpha_pri[1]^2) priors
//   - beta has N(beta_pri[0], beta_pri[1]^2) base prior
//   - (beta_M - beta) ~ N(0, sigma_comm^2 I)  [commensurate prior]
// PC-prior on sigma_comm: P(sigma_comm > u) = alpha (passed as lambdaSigmaComm = -log(alpha)/u).
//
// Designed for empirical-Bayes selection of sigma_comm via marginal likelihood:
// alpha, alpha_M, beta, beta_M are typically declared as random effects in R, so TMB
// Laplace-marginalizes them and the outer optimizer maximizes the marginal likelihood
// over (log_tauEps, log_sigma_comm) [+ (log_tau, logit_phi) when BYM2 is enabled].

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
  DATA_SCALAR( lambdaSigmaComm );  // PC prior rate on sigma_comm
  DATA_SCALAR( options );

  // Outer hyperparameters
  PARAMETER( log_tau );
  PARAMETER( logit_phi );
  PARAMETER( log_tauEps );
  PARAMETER( log_sigma_comm );

  // Latent / inner (typically declared random in R for Laplace marginalization)
  PARAMETER( alpha );
  PARAMETER( alpha_M );
  PARAMETER_VECTOR( beta );
  PARAMETER_VECTOR( beta_M );
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
  int nBeta = beta.size();

  Type tau = exp(log_tau);
  Type phi = Type(1.0)/(Type(1.0) + exp(-logit_phi));
  Type sigmaEps = exp(Type(-0.5) * log_tauEps);
  Type sigma_comm = exp(log_sigma_comm);

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

  // BYM2 density (vanishes when BYM2 params are mapped to NA at zero)
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

  // Hyperpriors
  Type logPriTau = log(lambdaTau) - log(Type(2.0)) - Type(1.5) * log_tau
                   - lambdaTau / sqrt(tau) + log_tau;
  jnll -= logPriTau;
  jnll -= dPCPriTau(log_tauEps, lambdaTauEps);

  // PC prior on sigma_comm: density lambda * exp(-lambda * sigma) on sigma > 0;
  // Jacobian for log_sigma adds + log(sigma_comm) = + log_sigma_comm.
  Type logPriSigmaComm = log(lambdaSigmaComm) - lambdaSigmaComm * sigma_comm
                         + log_sigma_comm;
  jnll -= logPriSigmaComm;

  // Prior for phi
  if(uniformPhiPrior) {
    jnll -= -logit_phi - Type(2.0) * log(Type(1.0) + exp(-logit_phi));
  } else {
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

  // Independent priors on alpha and alpha_M (no commensurate shrinkage on intercepts)
  jnll -= dnorm(alpha,   alpha_pri(0), alpha_pri(1), true);
  jnll -= dnorm(alpha_M, alpha_pri(0), alpha_pri(1), true);

  // Base prior on beta (DHS-side); beta_M's marginal prior is determined by beta + commensurate
  for(int i = 0; i < nBeta; i++) {
    jnll -= dnorm(beta(i), beta_pri(0), beta_pri(1), true);
  }
  // Commensurate prior: (beta_M - beta) ~ N(0, sigma_comm^2)
  for(int i = 0; i < nBeta; i++) {
    jnll -= dnorm(beta_M(i) - beta(i), Type(0.0), sigma_comm, true);
  }

  // Linear predictors using survey-specific alpha and beta
  vector<Type> fe_urb_m = X_betaUrbanMICS * beta_M;
  vector<Type> fe_rur_m = X_betaRuralMICS * beta_M;
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
                         Type alpha_use,
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
        base_eta(i, k) = alpha_use + fe(nObs * k + i) + w_bym2(areaidx(i));
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

  matrix<Type> wUrbM(nUrbM, KUrbM);
  for(int i=0;i<nUrbM;i++) for(int k=0;k<KUrbM;k++) wUrbM(i,k) = wUrbanMICS(i,k);
  matrix<Type> wRurM(nRurM, KRurM);
  for(int i=0;i<nRurM;i++) for(int k=0;k<KRurM;k++) wRurM(i,k) = wRuralMICS(i,k);
  matrix<Type> wUrbD(nUrbD, KUrbD);
  for(int i=0;i<nUrbD;i++) for(int k=0;k<KUrbD;k++) wUrbD(i,k) = wUrbanDHS(i,k);
  matrix<Type> wRurD(nRurD, KRurD);
  for(int i=0;i<nRurD;i++) for(int k=0;k<KRurD;k++) wRurD(i,k) = wRuralDHS(i,k);

  add_dataset(nUrbM, KUrbM, alpha_M, y_iUrbanMICS, n_iUrbanMICS, areaidxlocUrbanMICS,
              fe_urb_m, wUrbM, lchoose_urban_mics);
  add_dataset(nRurM, KRurM, alpha_M, y_iRuralMICS, n_iRuralMICS, areaidxlocRuralMICS,
              fe_rur_m, wRurM, lchoose_rural_mics);
  add_dataset(nUrbD, KUrbD, alpha,   y_iUrbanDHS, n_iUrbanDHS, areaidxlocUrbanDHS,
              fe_urb_d, wUrbD, lchoose_urban_dhs);
  add_dataset(nRurD, KRurD, alpha,   y_iRuralDHS, n_iRuralDHS, areaidxlocRuralDHS,
              fe_rur_d, wRurD, lchoose_rural_dhs);

  if(options == 1) {
    ADREPORT(log_tau);
    ADREPORT(logit_phi);
    ADREPORT(log_sigma_comm);
  }

  REPORT(w_bym2);
  REPORT(u_bym2);
  REPORT(log_tau);
  REPORT(logit_phi);
  REPORT(log_tauEps);
  REPORT(log_sigma_comm);
  REPORT(alpha);
  REPORT(alpha_M);
  REPORT(beta);
  REPORT(beta_M);

  return jnll;
}
