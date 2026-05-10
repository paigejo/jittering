// Marginalized BYM2 model for M_M (MICS-only) data
// Uses single n-dimensional Epsilon_bym2 with eigendecomposition V_bym2
// Eliminates the (w,u) funnel from the separable parameterization
//
// include libraries
#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
using namespace density;
using Eigen::SparseMatrix;

// helper function for detecting NAs in the data supplied from R
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// PC prior on BYM2 tau (precision) parameter (on log scale)
template<class Type>
Type dPCPriTau(Type logTau, Type lambda)
{
  Type tau = exp(logTau);
  
  Type ldensity = log(lambda) - log(2.0) -Type(1.5)*logTau - lambda/sqrt(tau);
  
  // add in log Jacobian
  Type ljacobian = logTau;
  
  // get final log density
  Type lPriorTau = ldensity + ljacobian;
  
  return lPriorTau;
}

///////////////////////////
// the main function     //
// to calculate the jnll //
///////////////////////////
template<class Type>
Type objective_function<Type>::operator() ()
{
  
  // ~~~~~~~~~------------------------------------------------------~~
  // FIRST, we define params/values/data that will be passed in from R
  // ~~~~~~~~~~~------------------------------------------------------
  
  // Data (all except for X_ij is a vector of length num_i)
  DATA_VECTOR( y_iUrbanMICS );   // obs per binomial experiment at point i (clust)
  DATA_VECTOR( y_iRuralMICS );
  DATA_VECTOR( n_iUrbanMICS );   // Trials per cluster
  DATA_VECTOR( n_iRuralMICS );
  DATA_IVECTOR(areaidxlocUrbanMICS); // vector of areal indices
  DATA_IVECTOR(areaidxlocRuralMICS);
  DATA_MATRIX( X_betaUrbanMICS );  // design matrix
  DATA_MATRIX( X_betaRuralMICS );
  DATA_ARRAY( wUrbanMICS ); // nObsUrban x nIntegrationPointsUrban weight matrix
  DATA_ARRAY( wRuralMICS );
  
  DATA_MATRIX( V_bym2 ); // Q_bym2 = V Lambda V^T is the eigendecomposition of Q_bym2
  
  // prior parameters
  DATA_VECTOR( alpha_pri );
  DATA_VECTOR( beta_pri );
  
  // BYM2 and prior precomputed values
  DATA_SCALAR( tr );
  DATA_VECTOR( gammaTildesm1 );
  
  DATA_SCALAR( lambdaPhi );
  DATA_SCALAR( lambdaTau );
  DATA_SCALAR( lambdaTauEps );
  
  // Options
  DATA_SCALAR( options );
  // options == 1 : adreport spatial params
  
  // Fixed effects
  PARAMETER( log_tau );
  PARAMETER( logit_phi );
  PARAMETER( log_tauEps );
  
  // Intercept (explicit, not absorbed into spatial field)
  PARAMETER( alpha );
  PARAMETER_VECTOR( beta );
  
  // BYM2 effects: marginalized n-dimensional representation
  // Epsilon_bym2 ~ N(0, (1/tau) * [(1-phi)*I + phi*Q_besag^+])
  PARAMETER_VECTOR( Epsilon_bym2 );
  PARAMETER_VECTOR( nuggetUrbMICS );
  PARAMETER_VECTOR( nuggetRurMICS );
  
  // ~~~~~~~~~------------------------------------------------~~
  // SECOND, we define all other objects that we need internally
  // ~~~~~~~~~------------------------------------------------~~
  
  int num_iUrbanMICS = y_iUrbanMICS.size();
  int num_iRuralMICS = y_iRuralMICS.size();
  int nAreas = Epsilon_bym2.size();
  int n_integrationPointsUrbanMICS = wUrbanMICS.cols();
  int n_integrationPointsRuralMICS = wRuralMICS.cols();
  
  // Transform parameters
  Type sqrtTau = exp(Type(0.5) * log_tau);
  Type tau = exp(log_tau);
  Type phi = Type(1.0)/(Type(1.0) + exp(-logit_phi));
  Type sigmaEps = exp(Type(-0.5) * log_tauEps);
  
  // objective function -- joint negative log-likelihood/posterior
  Type jnll = 0.0;
  
  // ~~~~~~~~~------------------------------------------------~~-
  // THIRD, we calculate the contribution to the likelihood from:
  // 1) BYM2 GMRF density
  // 2) priors
  // 3) likelihood
  // ~~~~~~~~~------------------------------------------------~~-
  
  /////////
  // (1) //
  /////////
  
  // Marginalized BYM2 density using eigendecomposition
  // Epsilon ~ N(0, (1/tau) * [(1-phi)*I + phi*Q^+])
  // Precision: tau * [(1-phi)*I + phi*Q^+]^{-1}
  //          = tau * V * diag(1/(1 + phi*gammaTildesm1)) * V^T
  // Quadratic form: sum_i (sqrt(tau) * (Eps^T V)_i / sqrt(1 + phi*gammaTildesm1_i))^2
  
  matrix<Type> transformedEpsilon = Epsilon_bym2.matrix().transpose() * V_bym2;
  Type quad = 0.0;
  for(int i = 0; i < transformedEpsilon.size(); i++) {
    quad += pow(transformedEpsilon(0,i) * sqrtTau / sqrt(Type(1.0) + phi*gammaTildesm1(i)), 2.0);
  }
  
  // Log determinant of covariance: log|(1/tau)[(1-phi)I + phi*Q^+]|
  //   = sum(log(1 + phi*gammaTildesm1)) + n*(-log_tau)
  Type logDet = Type(0.0);
  for(int i = 0; i < gammaTildesm1.size(); i++) {
    logDet += log(Type(1.0) + phi*gammaTildesm1(i));
  }
  Type logDetTau = logDet + Type(nAreas) * (-log_tau);
  
  // BYM2 log density (dropping 2*pi constant)
  Type bym2LogLik = Type(-0.5) * logDetTau;
  bym2LogLik += Type(-0.5) * quad;
  jnll -= bym2LogLik;
  
  // add in iid nugget model
  Type nuggetLogLik = 0.0;
  for(int i = 0; i< nuggetUrbMICS.size(); i++) {
    nuggetLogLik += dnorm(nuggetUrbMICS(i), Type(0), sigmaEps, true);
  }
  for(int i = 0; i< nuggetRurMICS.size(); i++) {
    nuggetLogLik += dnorm(nuggetRurMICS(i), Type(0), sigmaEps, true);
  }
  jnll -= nuggetLogLik;
  
  /////////
  // (2) //
  /////////
  // Prior contributions to joint likelihood
  
  // PC prior for tau
  Type logPriTauEps = dPCPriTau(log_tauEps, lambdaTauEps);
  
  Type firstPt = log(lambdaTau) - log(2.0);
  Type secondPt = -Type(1.5)*log_tau;
  Type thirdPt = - lambdaTau/sqrt(tau);
  Type ldensityTau = firstPt + secondPt + thirdPt;
  Type ljacobianTau = log_tau;
  Type logPriTau = ldensityTau + ljacobianTau;
  
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
  
  jnll -= logPriTau;
  jnll -= logPriTauEps;
  jnll -= logPriPhi;
  
  // prior for intercept
  jnll -= dnorm(alpha, alpha_pri(0), alpha_pri(1), true);
  
  // prior for covariates
  for(int i = 0; i < beta.size(); i++) {
    jnll -= dnorm(beta(i), beta_pri(0), beta_pri(1), true);
  }
  
  /////////
  // (3) //
  /////////
  // jnll contribution from each datapoint i
  
  // Fixed effects at data points (includes alpha)
  vector<Type> fe_iUrbanMICS = X_betaUrbanMICS * beta;
  vector<Type> fe_iRuralMICS = X_betaRuralMICS * beta;
  
  vector<Type> latentFieldUrbMICS(num_iUrbanMICS * n_integrationPointsUrbanMICS);
  vector<Type> latentFieldRurMICS(num_iRuralMICS * n_integrationPointsRuralMICS);
  int thisIndex;
  Type thisWeight;
  Type thislik;
  array<Type> liksUrbMICS(num_iUrbanMICS, n_integrationPointsUrbanMICS);
  array<Type> liksRurMICS(num_iRuralMICS, n_integrationPointsRuralMICS);
  
  for (int obsI = 0; obsI < num_iUrbanMICS; obsI++) {
    thislik = 0.0;
    
    for (int intI = 0; intI < n_integrationPointsUrbanMICS; intI++) {
      thisIndex = num_iUrbanMICS * intI + obsI;
      
      // latent field: alpha + X*beta + Epsilon_bym2[area] + nugget
      latentFieldUrbMICS(thisIndex) = alpha + fe_iUrbanMICS(thisIndex) + Epsilon_bym2(areaidxlocUrbanMICS(obsI)) + nuggetUrbMICS(obsI);
      
      thisWeight = wUrbanMICS(obsI,intI);
      
      if(thisWeight > 0.0) {
        liksUrbMICS(obsI, intI) = dbinom_robust( y_iUrbanMICS(obsI), n_iUrbanMICS(obsI), latentFieldUrbMICS(thisIndex), false);
        thislik += thisWeight*liksUrbMICS(obsI, intI);
      }
      else {
        liksUrbMICS(obsI, intI) = thisWeight;
      }
      
    }
    
    if(thislik > 0.0) {
      jnll -= log(thislik);
    }
  }
  
  for (int obsI = 0; obsI < num_iRuralMICS; obsI++) {
    thislik = 0.0;
    
    for (int intI = 0; intI < n_integrationPointsRuralMICS; intI++) {
      thisIndex = num_iRuralMICS * intI + obsI;
      
      latentFieldRurMICS(thisIndex) = alpha + fe_iRuralMICS(thisIndex) + Epsilon_bym2(areaidxlocRuralMICS(obsI)) + nuggetRurMICS(obsI);
      
      thisWeight = wRuralMICS(obsI,intI);
      
      if(thisWeight > 0.0) {
        liksRurMICS(obsI, intI) = dbinom_robust( y_iRuralMICS(obsI), n_iRuralMICS(obsI), latentFieldRurMICS(thisIndex), false);
        thislik += thisWeight*liksRurMICS(obsI, intI);
      }
      else {
        liksRurMICS(obsI, intI) = thisWeight;
      }
      
    }
    
    if(thislik > 0.0) {
      jnll -= log(thislik);
    }
  }
  
  // ~~~~~~~~~~~
  // ADREPORT
  // ~~~~~~~~~~~
  if(options==1){
    ADREPORT(log_tau);
    ADREPORT(logit_phi);
  }
  
  REPORT(Epsilon_bym2);
  REPORT(log_tau);
  REPORT(logit_phi);
  REPORT(log_tauEps);
  REPORT(alpha);
  REPORT(beta);
  REPORT(quad);
  REPORT(logDetTau);
  REPORT(bym2LogLik);
  REPORT(nuggetLogLik);
  REPORT(logPriTau);
  REPORT(ldensityTau);
  REPORT(ljacobianTau);
  REPORT(firstPt);
  REPORT(secondPt);
  REPORT(thirdPt);
  REPORT(logPriTauEps);
  REPORT(KLD);
  REPORT(lexpDensity);
  REPORT(ljacobian);
  REPORT(d);
  REPORT(logPriPhi);
  REPORT(liksUrbMICS);
  REPORT(liksRurMICS);
  REPORT(latentFieldUrbMICS);
  REPORT(latentFieldRurMICS);
  REPORT(fe_iUrbanMICS);
  REPORT(fe_iRuralMICS);
  
  
  return jnll;
  
}
