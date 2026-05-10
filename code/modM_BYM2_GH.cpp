// BYM2 constrained (2n-2) model with nuggets integrated out via Gauss-Hermite quadrature
// w_bym2Free: n-1 free params (nth = -sum), u_bym2Free: n-1 free params (nth = -sum)
// No nugget parameters — marginalized by GH
// All covariates via PARAMETER_VECTOR(beta), plus explicit PARAMETER(alpha)
// Inner (random): alpha, beta, w_bym2Free, u_bym2Free
// Outer: log_tau, logit_phi, log_tauEps

#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
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

// Log-sum-exp for tmbutils::vector
template<class Type>
Type logSumExp(vector<Type> x) {
  Type max_x = x.maxCoeff();
  Type sum_exp = Type(0.0);
  for(int i = 0; i < x.size(); i++) {
    sum_exp += exp(x(i) - max_x);
  }
  return max_x + log(sum_exp);
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

  // NO nugget parameters — integrated out by GH

  // ── Dimensions ──
  int num_iUrbanMICS = y_iUrbanMICS.size();
  int num_iRuralMICS = y_iRuralMICS.size();
  int nFree = w_bym2Free.size();
  int nAreas = nFree + 1;
  int n_integrationPointsUrbanMICS = wUrbanMICS.cols();
  int n_integrationPointsRuralMICS = wRuralMICS.cols();
  int Q = gh_nodes.size();

  // ── Transforms ──
  Type tau = exp(log_tau);
  Type phi = Type(1.0)/(Type(1.0) + exp(-logit_phi));
  Type sigmaEps = exp(Type(-0.5) * log_tauEps);

  // Reconstruct full n-vectors with sum-to-zero
  vector<Type> w_bym2(nAreas);
  vector<Type> u_bym2(nAreas);
  Type wSum = 0.0;
  Type uSum = 0.0;
  for(int i = 0; i < nFree; i++) {
    w_bym2(i) = w_bym2Free(i);
    u_bym2(i) = u_bym2Free(i);
    wSum += w_bym2Free(i);
    uSum += u_bym2Free(i);
  }
  w_bym2(nFree) = -wSum;
  u_bym2(nFree) = -uSum;

  Type jnll = 0.0;

  // ═══════════════════════════════════
  // (1) BYM2 GMRF density
  // ═══════════════════════════════════

  // tau/(1-phi) w^T w
  Type quadW = 0.0;
  for(int i = 0; i < nAreas; i++) {
    quadW += w_bym2(i) * w_bym2(i);
  }
  quadW = quadW * tau/(Type(1.0)-phi);

  // u^T Q u + phi/(1-phi) u^T u
  matrix<Type> transformedU = Q_bym2 * u_bym2.matrix();
  Type quadU = 0.0;
  Type fac = phi/(Type(1.0)-phi);
  for(int i = 0; i < nAreas; i++) {
    quadU += transformedU(i) * u_bym2(i) + fac * u_bym2(i) * u_bym2(i);
  }

  // -2 sqrt(phi*tau)/(1-phi) u^T w
  Type quadWU = 0.0;
  for(int i = 0; i < nAreas; i++) {
    quadWU += u_bym2(i) * w_bym2(i);
  }
  quadWU = quadWU * (-Type(2.0) * sqrt(phi * tau)/(Type(1.0)-phi));

  Type quadSum = quadW + quadU + quadWU;

  // Log determinant (constrained: n-1)
  Type logDet = Type(0.0);
  for(int i = 0; i < gammaTildesm1.size(); i++) {
    logDet += log(Type(1.0) + phi*gammaTildesm1(i));
  }
  Type logDetTau = Type(nFree) * log((Type(1.0) - phi)/tau);

  Type bym2LogLik = Type(-0.5) * logDetTau;
  bym2LogLik += Type(-0.5) * quadSum;
  jnll -= bym2LogLik;

  // ═══════════════════════════════════
  // (2) Priors
  // ═══════════════════════════════════

  // PC prior for tau
  Type firstPt = log(lambdaTau) - log(2.0);
  Type secondPt = -Type(1.5)*log_tau;
  Type thirdPt = - lambdaTau/sqrt(tau);
  Type logPriTau = firstPt + secondPt + thirdPt + log_tau;
  jnll -= logPriTau;

  // PC prior for tauEps
  jnll -= dPCPriTau(log_tauEps, lambdaTauEps);

  // PC prior for phi
  Type n = gammaTildesm1.size();
  Type KLD = Type(0.5) * (phi * tr - phi * n - logDet);
  Type d = sqrt(Type(2.0) * KLD);
  Type lexpDensity = log(lambdaPhi) - lambdaPhi * d;
  Type sumVal = 0.0;
  for(int i = 0; i < n; i++) {
    sumVal += gammaTildesm1(i)/(1 + gammaTildesm1(i) * phi);
  }
  Type ljacobian = - log(d) - log(2.0) + log(tr - n - sumVal);
  ljacobian = ljacobian - logit_phi - Type(2.0) * log(Type(1.0) + exp(-logit_phi));
  Type logPriPhi = lexpDensity + ljacobian;
  jnll -= logPriPhi;

  // Prior for intercept
  jnll -= dnorm(alpha, alpha_pri(0), alpha_pri(1), true);

  // Prior for covariates
  for(int i = 0; i < beta.size(); i++) {
    jnll -= dnorm(beta(i), beta_pri(0), beta_pri(1), true);
  }

  // ═══════════════════════════════════
  // (3) GH-integrated likelihood
  // ═══════════════════════════════════

  // Precompute fixed effects at all integration points
  vector<Type> fe_iUrbanMICS = X_betaUrbanMICS * beta;
  vector<Type> fe_iRuralMICS = X_betaRuralMICS * beta;

  // GH abscissae: eps_q = sqrt(2) * sigmaEps * node_q
  vector<Type> eps_gh(Q);
  for(int q = 0; q < Q; q++) {
    eps_gh(q) = sqrt(Type(2.0)) * sigmaEps * gh_nodes(q);
  }
  Type log_inv_sqrt_pi = -Type(0.5) * log(M_PI);

  int thisIndex;
  Type thisWeight;

  // ── Urban MICS ──
  for (int obsI = 0; obsI < num_iUrbanMICS; obsI++) {
    vector<Type> log_terms(Q);

    for(int q = 0; q < Q; q++) {
      Type mix_lik = Type(0.0);

      for (int intI = 0; intI < n_integrationPointsUrbanMICS; intI++) {
        thisIndex = num_iUrbanMICS * intI + obsI;
        thisWeight = wUrbanMICS(obsI, intI);

        if(thisWeight > Type(0.0)) {
          Type eta = alpha + fe_iUrbanMICS(thisIndex)
                     + w_bym2(areaidxlocUrbanMICS(obsI)) + eps_gh(q);
          mix_lik += thisWeight * dbinom_robust(y_iUrbanMICS(obsI),
                                                n_iUrbanMICS(obsI), eta, false);
        }
      }

      if(mix_lik > Type(0.0)) {
        log_terms(q) = log(gh_weights(q)) + log(mix_lik);
      } else {
        log_terms(q) = Type(-1000.0);
      }
    }

    jnll -= log_inv_sqrt_pi + logSumExp(log_terms);
  }

  // ── Rural MICS ──
  for (int obsI = 0; obsI < num_iRuralMICS; obsI++) {
    vector<Type> log_terms(Q);

    for(int q = 0; q < Q; q++) {
      Type mix_lik = Type(0.0);

      for (int intI = 0; intI < n_integrationPointsRuralMICS; intI++) {
        thisIndex = num_iRuralMICS * intI + obsI;
        thisWeight = wRuralMICS(obsI, intI);

        if(thisWeight > Type(0.0)) {
          Type eta = alpha + fe_iRuralMICS(thisIndex)
                     + w_bym2(areaidxlocRuralMICS(obsI)) + eps_gh(q);
          mix_lik += thisWeight * dbinom_robust(y_iRuralMICS(obsI),
                                                n_iRuralMICS(obsI), eta, false);
        }
      }

      if(mix_lik > Type(0.0)) {
        log_terms(q) = log(gh_weights(q)) + log(mix_lik);
      } else {
        log_terms(q) = Type(-1000.0);
      }
    }

    jnll -= log_inv_sqrt_pi + logSumExp(log_terms);
  }

  // ── Reports ──
  if(options==1){
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
  REPORT(logPriPhi);

  return jnll;
}
