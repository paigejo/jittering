// IID spatial + nugget model — nuggets integrated out via Gauss-Hermite quadrature
// No PARAMETER_VECTOR for nuggets: they are analytically marginalized
// 
// For FE-only (u_spatial mapped out): zero random effects, pure optimization
// For IID: only u_spatial + alpha + beta in inner Newton (~40 params)

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

// Log-sum-exp for a tmbutils::vector
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
  
  // Gauss-Hermite quadrature nodes and weights
  DATA_VECTOR( gh_nodes );    // Q nodes for standard Gauss-Hermite
  DATA_VECTOR( gh_weights );  // Q weights
  
  // Prior parameters
  DATA_VECTOR( alpha_pri );
  DATA_VECTOR( beta_pri );
  DATA_SCALAR( lambdaTau );
  DATA_SCALAR( lambdaTauEps );
  
  // Hyperparameters
  PARAMETER( log_tau );
  PARAMETER( log_tauEps );
  
  // Intercept
  PARAMETER( alpha );
  
  // Split beta
  PARAMETER( beta_urban );
  PARAMETER( beta_normPop );
  
  // IID spatial random effects
  PARAMETER_VECTOR( u_spatial );
  
  // NO nugget parameters — integrated out by GH quadrature
  
  // Dimensions
  int num_iUrbanMICS = y_iUrbanMICS.size();
  int num_iRuralMICS = y_iRuralMICS.size();
  int n_integrationPointsUrbanMICS = wUrbanMICS.cols();
  int n_integrationPointsRuralMICS = wRuralMICS.cols();
  int Q = gh_nodes.size();
  
  // Transform
  Type tau = exp(log_tau);
  Type sigma_u = exp(Type(-0.5) * log_tau);
  Type sigmaEps = exp(Type(-0.5) * log_tauEps);
  
  Type jnll = 0.0;
  
  // (1) IID spatial random effects prior
  for(int i = 0; i < u_spatial.size(); i++) {
    jnll -= dnorm(u_spatial(i), Type(0), sigma_u, true);
  }
  
  // (2) Hyperparameter priors
  jnll -= dPCPriTau(log_tau, lambdaTau);
  jnll -= dPCPriTau(log_tauEps, lambdaTauEps);
  
  // (3) Fixed effect priors
  jnll -= dnorm(alpha, alpha_pri(0), alpha_pri(1), true);
  jnll -= dnorm(beta_urban, beta_pri(0), beta_pri(1), true);
  jnll -= dnorm(beta_normPop, beta_pri(0), beta_pri(1), true);
  
  // (4) Precompute fixed effects for all integration points
  vector<Type> fe_iUrbanMICS = X_betaUrbanMICS.col(0) * beta_urban +
                                X_betaUrbanMICS.col(1) * beta_normPop;
  vector<Type> fe_iRuralMICS = X_betaRuralMICS.col(0) * beta_urban +
                                X_betaRuralMICS.col(1) * beta_normPop;
  
  // Precompute GH abscissae: eps_q = sqrt(2) * sigmaEps * node_q
  vector<Type> eps_gh(Q);
  for(int q = 0; q < Q; q++) {
    eps_gh(q) = sqrt(Type(2.0)) * sigmaEps * gh_nodes(q);
  }
  
  // Constant: log(1/sqrt(pi))
  Type log_inv_sqrt_pi = -Type(0.5) * log(M_PI);
  
  int thisIndex;
  Type thisWeight;
  
  // ── Urban MICS ──
  for (int obsI = 0; obsI < num_iUrbanMICS; obsI++) {
    // For each obs, compute log L_i = log integral over eps of [mixture likelihood * N(eps)] deps
    // Using GH: L_i = (1/sqrt(pi)) * sum_q w_q^GH * [sum_k w_k * Binom(y|n, logit^-1(eta_k + eps_q))]
    
    // log-sum-exp over GH nodes
    vector<Type> log_terms(Q);
    
    for(int q = 0; q < Q; q++) {
      // Integration-point mixture for this GH node
      Type mix_lik = Type(0.0);
      
      for (int intI = 0; intI < n_integrationPointsUrbanMICS; intI++) {
        thisIndex = num_iUrbanMICS * intI + obsI;
        thisWeight = wUrbanMICS(obsI, intI);
        
        if(thisWeight > Type(0.0)) {
          Type eta = alpha + fe_iUrbanMICS(thisIndex) 
                     + u_spatial(areaidxlocUrbanMICS(obsI)) + eps_gh(q);
          mix_lik += thisWeight * dbinom_robust(y_iUrbanMICS(obsI), 
                                                n_iUrbanMICS(obsI), eta, false);
        }
      }
      
      // log(w_q^GH * mix_lik)
      if(mix_lik > Type(0.0)) {
        log_terms(q) = log(gh_weights(q)) + log(mix_lik);
      } else {
        log_terms(q) = Type(-1000.0);  // effectively zero
      }
    }
    
    // log L_i = log(1/sqrt(pi)) + logSumExp(log_terms)
    Type log_lik_i = log_inv_sqrt_pi + logSumExp(log_terms);
    jnll -= log_lik_i;
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
                     + u_spatial(areaidxlocRuralMICS(obsI)) + eps_gh(q);
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
    
    Type log_lik_i = log_inv_sqrt_pi + logSumExp(log_terms);
    jnll -= log_lik_i;
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
