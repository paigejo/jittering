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

// // PC prior on BYM2 phi parameter (on logit scale)
// template<class Type>
// Type dPCPriPhi(Type logitPhi, Type lambda, Type tr, vector<Type> gammaTildesm1, Type logDet, 
//                int give_log=0)
// {
//   Type phi = (Type) (1 / (1 + exp(-logitPhi)));
//   
//   // calculate log determinant, KLD(phi), and d(phi)
//   // Type logDet = sum(log(1 + phi*gammaTildesm1));
//   
//   Type n = gammaTildesm1.size();
//   Type KLD = 0.5 * (phi * tr - phi * n - logDet);
//   Type d = sqrt(2*KLD);
//   
//   // calculate exponential density of d(phi)
//   Type lexpDensity = log(lambda) - lambda * d;
//   
//   // add Jacobian with respect to d(phi) and logit(phi)
//   Type sumVal = 0;
//   for(int i = 0; i < n; i++) {
//     sumVal += gammaTildesm1(i)/(1 + gammaTildesm1(i) * phi);
//   }
//   //Type ljacobian = - log(d) - log(2) + log(abs(tr - n - sumVal));
//   Type ljacobian = - log(d) - log(2) + log(tr - n - sumVal); // this should be fine since argument should be positive
//   
//   ljacobian = ljacobian - logitPhi - 2*log(1 + exp(-logitPhi));
//   
//   // get final log density
//   Type lPriorPhi = lexpDensity + ljacobian;
//   
//   return lPriorPhi;
// }

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
  
  // normalization flag - used for speed-up
  // DATA_INTEGER( flag ); // flag == 0 => no data contribution added to jnll
  
  // Data (all except for X_ij is a vector of length num_i)
  DATA_VECTOR( y_iUrbanDHS );   // obs per binomial experiment at point i (clust)
  DATA_VECTOR( y_iRuralDHS );
  DATA_VECTOR( n_iUrbanDHS );   // Trials per cluster
  DATA_VECTOR( n_iRuralDHS );
  DATA_SPARSE_MATRIX( AprojUrbanDHS ); // nObsUrban x nArea matrix with ij-th entry = 1 if pt i associated with area j and 0 o.w.
  DATA_SPARSE_MATRIX( AprojRuralDHS ); // nObsRural x nArea matrix with ij-th entry = 1 if pt i associated with area j and 0 o.w.
  DATA_MATRIX( X_betaUrbanDHS );  // nObsUrban x nPar design matrix
  DATA_MATRIX( X_betaRuralDHS );  // nObsRural x nPar design matrix
  
  DATA_MATRIX( V_bym2 ); // Q_bym2 = V Lamda V^T is the eigendecomposition of Q_bym2
  DATA_SPARSE_MATRIX( Q_bym2 );
  
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
  PARAMETER( alpha ); // Intercept
  PARAMETER_VECTOR( beta ); // fixed effect/covariate effect sizes
  // Log of INLA tau param (generalized precision of BYM2)
  PARAMETER( log_tau );
  // Logit of phi (proportion of variance that is structured)
  PARAMETER( logit_phi );
  // log of nugget precision
  PARAMETER( log_tauEps );
  
  // BYM2 effects (structured + unstructured)
  // NOTE: may need to modify parameterization later
  PARAMETER_VECTOR( Epsilon_bym2 );
  PARAMETER_VECTOR( nuggetUrbDHS );
  PARAMETER_VECTOR( nuggetRurDHS );
  
  // ~~~~~~~~~------------------------------------------------~~
  // SECOND, we define all other objects that we need internally
  // ~~~~~~~~~------------------------------------------------~~
  
  // calculated values from data
  int num_iUrbanDHS = y_iUrbanDHS.size();   // Number of urban data points in space
  int num_iRuralDHS = y_iRuralDHS.size();   // Number of rural data points in space
  
  // objective function -- joint negative log-likelihood/posterior
  Type jnll = 0.0;
  
  // Transform some of our parameters
  Type sqrtTau = exp(Type(0.5) * log_tau);
  Type phi = 1.0/(1.0 + exp(-logit_phi));
  Type sigmaEps = exp(-0.5 * log_tauEps);
  
  vector<Type> projepsilon_iUrbanDHS(num_iUrbanDHS);
  vector<Type> projepsilon_iRuralDHS(num_iRuralDHS);
  
  // Project GP approx from mesh points to data points
  projepsilon_iUrbanDHS = (AprojUrbanDHS * Epsilon_bym2.matrix()).col(0);
  projepsilon_iRuralDHS = (AprojRuralDHS * Epsilon_bym2.matrix()).col(0);
  
  // ~~~~~~~~~------------------------------------------------~~-
  // THIRD, we calculate the contribution to the likelihood from:
  // 1) GP field first, for 'flag' normalization purposes
  // 2) priors
  // 3) GP field
  // ~~~~~~~~~------------------------------------------------~~-
  
  /////////
  // (1) //
  /////////
  
  // add in BYM2 latent model component to the posterior
  
  // Calculate the quadratic form Eps tau [(1-phi) I + phi Q_besag^+]^(-1) Eps^T
  // = Eps tau [I + phi (Q_besag^+ - I)]^(-1) Eps^T
  // = Eps tau V [I + phi (GammaTilde - I)]^(-1) V^T Eps^T
  // i.e. the sum of squares of tau^0.5 Eps V diag(1/sqrt(1 + phi*gammaTildesm1))
  matrix<Type> transformedEpsilon = Epsilon_bym2.matrix().transpose() * V_bym2;
  Type quad = 0.0;
  for(int i = 0; i < transformedEpsilon.size(); i++) {
    quad += pow(transformedEpsilon(0,i) * sqrtTau / sqrt(1.0 + phi*gammaTildesm1(i)), 2.0);
  }
  
  // calculate log determinant of (1/tau) [(1-phi) I + phi Q_besag^+]
  Type logDet = Type(0);
  for(int i = 0; i < gammaTildesm1.size(); i++) {
    logDet += log(1.0 + phi*gammaTildesm1(i));
  }
  Type logDetTau = logDet + gammaTildesm1.size() * (-log_tau);
  
  // add quadratic form and normalizing constant to the negative posterior log density
  // first add denominator
  Type bym2LogLik = -0.5 * (gammaTildesm1.size() * log(2.0*M_PI) + logDetTau);
  // now add exponent
  bym2LogLik += -0.5 * quad;
  // add all to the negative posterior log density
  jnll -= bym2LogLik;
  
  // add in iid nugget model
  Type nuggetLogLik = 0.0;
  for(int i = 0; i< nuggetUrbDHS.size(); i++) {
    nuggetLogLik += dnorm(nuggetUrbDHS(i), Type(0.0), sigmaEps, true); // N(mean, sd);
  }
  for(int i = 0; i< nuggetRurDHS.size(); i++) {
    nuggetLogLik += dnorm(nuggetRurDHS(i), Type(0.0), sigmaEps, true); // N(mean, sd);
  }
  jnll -= nuggetLogLik;
  
  /////////
  // (2) //
  /////////
  // Prior contributions to joint likelihood
  
  // add in priors for BYM2 and nugget precision parameters
  Type logPriTauEps = dPCPriTau(log_tauEps, lambdaTauEps);
  // Type logPriTau = dPCPriTau(log_tau, lambdaTau);
  // Type logPriPhi = dPCPriPhi(logit_phi, lambdaPhi, tr, gammaTildesm1, logDet, true);
  
  Type tau = exp(log_tau);
  
  Type firstPt = log(lambdaTau) - log(2.0);
  Type secondPt = -Type(1.5)*log_tau;
  Type thirdPt = - lambdaTau/sqrt(tau);
  
  // Type ldensityTau = log(lambdaTau) - log(2) -Type(3/2)*log_tau - lambdaTau/sqrt(tau);
  Type ldensityTau = firstPt + secondPt + thirdPt;
  
  // add in log Jacobian
  Type ljacobianTau = log_tau;
  
  // get final log density
  Type logPriTau = ldensityTau + ljacobianTau;
  
  // calculate log determinant, KLD(phi), and d(phi)
  // Type logDet = sum(log(1 + phi*gammaTildesm1));
  
  Type n = gammaTildesm1.size();
  Type KLD = 0.5 * (phi * tr - phi * n - logDet);
  Type d = sqrt(2.0*KLD);
  
  // calculate exponential density of d(phi)
  Type lexpDensity = log(lambdaPhi) - lambdaPhi * d;
  
  // add Jacobian with respect to d(phi) and logit(phi)
  Type sumVal = 0.0;
  for(int i = 0; i < n; i++) {
    sumVal += gammaTildesm1(i)/(1.0 + gammaTildesm1(i) * phi);
  }
  //Type ljacobian = - log(d) - log(2) + log(abs(tr - n - sumVal));
  Type ljacobian = - log(d) - log(2.0) + log(tr - n - sumVal); // this should be fine since argument should be positive
  
  ljacobian = ljacobian - logit_phi - 2.0*log(1.0 + exp(-logit_phi));
  
  Type logPriPhi = lexpDensity + ljacobian;
  
  // get final log density
  jnll -= logPriTau;
  jnll -= logPriTauEps;
  jnll -= logPriPhi;
  
  // prior for intercept (removed to be same as INLA)
  // jnll -= dnorm(alpha, alpha_pri(0), alpha_pri(1), true); // N(mean, sd)
  
  // calculate fixed effects
  vector<Type> feUrban = X_betaUrbanDHS * beta;
  vector<Type> feRural = X_betaRuralDHS * beta;
  
  /////////
  // (3) //
  /////////
  // jnll contribution from each datapoint i
  vector<Type> latentFieldUrbDHS(num_iUrbanDHS);
  vector<Type> latentFieldRurDHS(num_iRuralDHS);
  vector<Type> liksUrbDHS(num_iUrbanDHS);
  vector<Type> liksRurDHS(num_iRuralDHS);
  
  for (int obsI = 0; obsI < num_iUrbanDHS; obsI++) {
    
    // latent field estimate at each obs
    latentFieldUrbDHS(obsI) = alpha + feUrban(obsI) + projepsilon_iUrbanDHS(obsI) + nuggetUrbDHS(obsI);
    
    // Uses the dbinom_robust function, which takes the logit probability
    liksUrbDHS(obsI) = dbinom_robust( y_iUrbanDHS(obsI), n_iUrbanDHS(obsI), latentFieldUrbDHS(obsI), false);
    jnll -= log(liksUrbDHS(obsI));
    
  } // for( obsI )
  
  for (int obsI = 0; obsI < num_iRuralDHS; obsI++) {
    
    // latent field estimate at each obs
    latentFieldRurDHS(obsI) = alpha + feRural(obsI) + projepsilon_iRuralDHS(obsI) + nuggetRurDHS(obsI);
    
    // Uses the dbinom_robust function, which takes the logit probability
    liksRurDHS(obsI) = dbinom_robust( y_iRuralDHS(obsI), n_iRuralDHS(obsI), latentFieldRurDHS(obsI), false);
    jnll -= log(liksRurDHS(obsI));
    
  } // for( obsI )
  
  // ~~~~~~~~~~~
  // ADREPORT: used to return estimates and cov for transforms?
  // ~~~~~~~~~~~
  if(options==1){
    ADREPORT(log_tau);
    ADREPORT(logit_phi);
  }
  
  REPORT(log_tau);
  REPORT(quad);
  REPORT(logDet);
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
  REPORT(liksUrbDHS);
  REPORT(liksRurDHS);
  REPORT(transformedEpsilon);
  REPORT(feUrban);
  REPORT(feRural);
  REPORT(latentFieldUrbDHS);
  REPORT(latentFieldRurDHS);
  
  return jnll;
  
}
