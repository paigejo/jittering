// IID spatial effects only model - simplified for testing
// M_M (MICS-only) data
// Replaces BYM2 with IID spatial random effects + cluster-level nugget
// Uses same integration point structure as the full models

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
  DATA_INTEGER( nAreas );
  
  // Prior parameters
  DATA_VECTOR( alpha_pri );
  DATA_VECTOR( beta_pri );
  DATA_SCALAR( lambdaTau );
  DATA_SCALAR( lambdaTauEps );
  
  // Fixed effects
  PARAMETER( log_tau );
  PARAMETER( log_tauEps );
  
  // Intercept
  PARAMETER( alpha );
  PARAMETER_VECTOR( beta );
  
  // IID spatial random effects
  PARAMETER_VECTOR( u_spatial );
  
  // Cluster-level nugget effects
  PARAMETER_VECTOR( nuggetUrbMICS );
  PARAMETER_VECTOR( nuggetRurMICS );
  
  // Dimensions
  int num_iUrbanMICS = y_iUrbanMICS.size();
  int num_iRuralMICS = y_iRuralMICS.size();
  int n_integrationPointsUrbanMICS = wUrbanMICS.cols();
  int n_integrationPointsRuralMICS = wRuralMICS.cols();
  
  // Transform
  Type tau = exp(log_tau);
  Type sigma_u = exp(Type(-0.5) * log_tau);
  Type sigmaEps = exp(Type(-0.5) * log_tauEps);
  
  Type jnll = 0.0;
  
  // (1) IID spatial random effects prior
  for(int i = 0; i < u_spatial.size(); i++) {
    jnll -= dnorm(u_spatial(i), Type(0), sigma_u, true);
  }
  
  // (2) Priors
  // PC prior for tau (IID spatial)
  jnll -= dPCPriTau(log_tau, lambdaTau);
  
  // PC prior for tauEps (nugget)
  jnll -= dPCPriTau(log_tauEps, lambdaTauEps);
  
  // Nugget priors
  for(int i = 0; i < nuggetUrbMICS.size(); i++) {
    jnll -= dnorm(nuggetUrbMICS(i), Type(0), sigmaEps, true);
  }
  for(int i = 0; i < nuggetRurMICS.size(); i++) {
    jnll -= dnorm(nuggetRurMICS(i), Type(0), sigmaEps, true);
  }
  
  // Prior for intercept
  jnll -= dnorm(alpha, alpha_pri(0), alpha_pri(1), true);
  
  // Prior for covariates
  for(int i = 0; i < beta.size(); i++) {
    jnll -= dnorm(beta(i), beta_pri(0), beta_pri(1), true);
  }
  
  // (3) Likelihood with integration points (same structure as Marg model)
  vector<Type> fe_iUrbanMICS = X_betaUrbanMICS * beta;
  vector<Type> fe_iRuralMICS = X_betaRuralMICS * beta;
  
  int thisIndex;
  Type thisWeight;
  Type thislik;
  
  // Urban MICS
  for (int obsI = 0; obsI < num_iUrbanMICS; obsI++) {
    thislik = 0.0;
    
    for (int intI = 0; intI < n_integrationPointsUrbanMICS; intI++) {
      thisIndex = num_iUrbanMICS * intI + obsI;
      
      Type eta = alpha + fe_iUrbanMICS(thisIndex) + u_spatial(areaidxlocUrbanMICS(obsI)) + nuggetUrbMICS(obsI);
      
      thisWeight = wUrbanMICS(obsI, intI);
      
      if(thisWeight > 0.0) {
        thislik += thisWeight * dbinom_robust(y_iUrbanMICS(obsI), n_iUrbanMICS(obsI), eta, false);
      }
    }
    
    if(thislik > 0.0) {
      jnll -= log(thislik);
    }
  }
  
  // Rural MICS
  for (int obsI = 0; obsI < num_iRuralMICS; obsI++) {
    thislik = 0.0;
    
    for (int intI = 0; intI < n_integrationPointsRuralMICS; intI++) {
      thisIndex = num_iRuralMICS * intI + obsI;
      
      Type eta = alpha + fe_iRuralMICS(thisIndex) + u_spatial(areaidxlocRuralMICS(obsI)) + nuggetRurMICS(obsI);
      
      thisWeight = wRuralMICS(obsI, intI);
      
      if(thisWeight > 0.0) {
        thislik += thisWeight * dbinom_robust(y_iRuralMICS(obsI), n_iRuralMICS(obsI), eta, false);
      }
    }
    
    if(thislik > 0.0) {
      jnll -= log(thislik);
    }
  }
  
  // Report
  REPORT(log_tau);
  REPORT(log_tauEps);
  REPORT(alpha);
  REPORT(beta);
  REPORT(u_spatial);
  REPORT(sigma_u);
  REPORT(sigmaEps);
  REPORT(nuggetUrbMICS);
  REPORT(nuggetRurMICS);
  
  return jnll;
}
