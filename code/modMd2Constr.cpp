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
  DATA_SPARSE_MATRIX( AprojUrbanDHS ); // (nObsUrban * nIntegrationPointsUrban) x nArea matrix with ij-th entry = 1 if pt i associated with area j and 0 o.w.
  DATA_SPARSE_MATRIX( AprojRuralDHS ); // (nObsRural * nIntegrationPointsUrban) x nArea matrix with ij-th entry = 1 if pt i associated with area j and 0 o.w.
  DATA_MATRIX( X_betaUrbanDHS );  // (nObsUrban * nIntegrationPointsUrban) x nPar design matrix. Indexed mod numObsUrban
  DATA_MATRIX( X_betaRuralDHS );  // first nObsRural rows correspond to first int pt
  DATA_ARRAY( wUrbanDHS ); // nObsUrban x nIntegrationPointsUrban weight matrix
  DATA_ARRAY( wRuralDHS ); // nObsRural x nIntegrationPointsRural weight matrix
  
  // Q_bym2 = V Lamda V^T is the eigendecomposition of Q_bym2
  // NOTE: in this version, V is n x (n-1) matrix
  DATA_MATRIX( V_bym2 );
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
  // Log of INLA tau param (generalized precision of BYM2)
  PARAMETER( log_tau );
  // Logit of phi (proportion of variance that is structured)
  PARAMETER( logit_phi );
  // log of nugget precision
  PARAMETER( log_tauEps );
  
  // Random effects
  PARAMETER( alpha ); // Intercept
  PARAMETER_VECTOR( beta ); // fixed effect/covariate effect sizes
  
  // BYM2 effects (structured + unstructured)
  // NOTE: Epsilon_bym2 is iid std Gaussian vector. 
  //       Later transformed to have correct BYM2 distribution
  PARAMETER_VECTOR( Epsilon_bym2 );
  PARAMETER_VECTOR( nuggetUrbDHS );
  PARAMETER_VECTOR( nuggetRurDHS );
  
  // ~~~~~~~~~------------------------------------------------~~
  // SECOND, we define all other objects that we need internally
  // ~~~~~~~~~------------------------------------------------~~
  
  // calculated values from data
  int num_iUrbanDHS = y_iUrbanDHS.size();   // Number of urban data points in space
  int num_iRuralDHS = y_iRuralDHS.size();   // Number of rural data points in space
  int n_integrationPointsUrbanDHS = wUrbanDHS.cols();    // number of integration points for each observation
  int n_integrationPointsRuralDHS = wRuralDHS.cols();
  
  // objective function -- joint negative log-likelihood/posterior
  Type jnll = 0.0;
  
  // Transform some of our parameters
  // Type sqrtTau = exp(Type(0.5) * log_tau);
  Type phi = 1.0/(1.0 + exp(-logit_phi));
  Type sigmaEps = exp(-0.5 * log_tauEps);
  
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
  
  // Get log density of untransformed Epsilon (easy, since iid std normal)
  Type bym2LogLik = 0;
  for(int i = 0; i < Epsilon_bym2.size(); i++) {
    bym2LogLik += dnorm(Epsilon_bym2(i), Type(0.0), Type(1.0), true);
  }
  // sum(dnorm(Epsilon_bym2, Type(0.0), Type(1.0), true)); // vectorized, but we need vector 1s and 0s
  
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
  
  // calculate log determinant of (1/tau) [(1-phi) I + phi Q_besag^+]
  Type logDet = Type(0);
  for(int i = 0; i < gammaTildesm1.size(); i++) {
    logDet += log(1.0 + phi*gammaTildesm1(i));
  }
  
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
  
  // prior for other covariates
  for(int i = 0; i < beta.size(); i++) {
    jnll -= dnorm(beta(i), beta_pri(0), beta_pri(1), true); // N(mean, sd)
  }
  
  /////////
  // (3) //
  /////////
  
  // transform Epsilon from n-1 dimensional vector to n-dimensional with sum to 0 constr:
  // first get sqrt of lambdas, the eigenvalues corresponding to the non-constant eigenvectors 
  // of the bym2 effect
  vector<Type> sqrtLambdas(gammaTildesm1.size()-1);
  for(int i = 0; i < gammaTildesm1.size()-1; i++) {
    sqrtLambdas(i) = sqrt((1.0 + phi * gammaTildesm1(i))/tau);
  }
  // now transform as transformedEpsilon = V Lambdas^1/2 Epsilon_bym2
  vector<Type> transformedEpsilon = V_bym2 * (Epsilon_bym2 * sqrtLambdas).matrix();
  
  vector<Type> projepsilon_iUrbanDHS(num_iUrbanDHS);
  vector<Type> projepsilon_iRuralDHS(num_iRuralDHS);
  
  // Project GP approx from mesh points to data points
  projepsilon_iUrbanDHS = (AprojUrbanDHS * transformedEpsilon.matrix()).col(0);
  projepsilon_iRuralDHS = (AprojRuralDHS * transformedEpsilon.matrix()).col(0);
  
  // Calculate fixed effects at data points
  vector<Type> fe_iUrbanDHS = alpha + X_betaUrbanDHS * beta;
  vector<Type> fe_iRuralDHS = alpha + X_betaRuralDHS * beta;
  
  // jnll contribution from each datapoint i
  // vector<Type> latentFieldUrbDHS = alpha + fe_iUrbanDHS + projepsilon_iUrbanDHS + nuggetUrbDHS;
  // vector<Type> latentFieldRurDHS = alpha + fe_iRuralDHS + projepsilon_iRuralDHS + nuggetRurDHS;
  // vector<Type> liksUrbDHS = dbinom_robust(y_iUrbanDHS, n_iUrbanDHS, latentFieldUrbDHS, false);
  // vector<Type> liksRurDHS = dbinom_robust(y_iRuralDHS, n_iRuralDHS, latentFieldRurDHS, false);
  // 
  // jnll -= sum(log(liksUrbDHS));
  // jnll -= sum(log(liksRurDHS));
  vector<Type> latentFieldUrbDHS(num_iUrbanDHS * n_integrationPointsUrbanDHS);
  vector<Type> latentFieldRurDHS(num_iRuralDHS * n_integrationPointsRuralDHS);
  int thisIndex;
  // Type thisLatentField;
  Type thisWeight;
  Type thislik;
  array<Type> liksUrbDHS(num_iUrbanDHS, n_integrationPointsUrbanDHS);
  array<Type> liksRurDHS(num_iRuralDHS, n_integrationPointsRuralDHS);
  
  for (int obsI = 0; obsI < num_iUrbanDHS; obsI++) {
    thislik = 0.0;
    
    for (int intI = 0; intI < n_integrationPointsUrbanDHS; intI++) {
      thisIndex = num_iUrbanDHS * intI + obsI;
      
      // latent field estimate at each obs
      // thisLatentField = fe_iUrbanDHS(thisIndex) + projepsilon_iUrbanDHS(obsI);
      latentFieldUrbDHS(thisIndex) = fe_iUrbanDHS(thisIndex) + projepsilon_iUrbanDHS(thisIndex) + nuggetUrbDHS(obsI);
      // latentFieldUrbDHS(thisIndex) = fe_iUrbanDHS(thisIndex) + projepsilon_iUrbanDHS(obsI);
      
      // and add data contribution to jnll
      // get integration weight
      thisWeight = wUrbanDHS(obsI,intI);
      
      // Uses the dbinom_robust function, which takes the logit probability
      if(thisWeight > 0.0) {
        // thislik += thisWeight*dbinom_robust( y_iUrbanDHS(obsI), n_iUrbanDHS(obsI), thisLatentField, false);
        liksUrbDHS(obsI, intI) = dbinom_robust( y_iUrbanDHS(obsI), n_iUrbanDHS(obsI), latentFieldUrbDHS(thisIndex), false);
        thislik += thisWeight*liksUrbDHS(obsI, intI);
      }
      else {
        liksUrbDHS(obsI, intI) = thisWeight;
      }
      
    } // for( intI )
    
    if(thislik > 0.0) {
      jnll -= log(thislik);
    }
  } // for( obsI )
  
  for (int obsI = 0; obsI < num_iRuralDHS; obsI++) {
    thislik = 0.0;
    
    for (int intI = 0; intI < n_integrationPointsRuralDHS; intI++) {
      thisIndex = num_iRuralDHS * intI + obsI;
      
      // latent field estimate at each obs
      // thisLatentField = fe_iRuralDHS(thisIndex) + projepsilon_iRuralDHS(obsI);
      latentFieldRurDHS(thisIndex) = fe_iRuralDHS(thisIndex) + projepsilon_iRuralDHS(thisIndex) + nuggetRurDHS(obsI);
      // latentFieldRurDHS(thisIndex) = fe_iRuralDHS(thisIndex) + projepsilon_iRuralDHS(obsI);
      
      // and add data contribution to jnll
      // get integration weight
      thisWeight = wRuralDHS(obsI,intI);
      
      // Uses the dbinom_robust function, which takes the logit probability
      if(thisWeight > 0.0) {
        // thislik += thisWeight*dbinom_robust( y_iRuralDHS(obsI), n_iRuralDHS(obsI), thisLatentField, false);
        liksRurDHS(obsI, intI) = dbinom_robust( y_iRuralDHS(obsI), n_iRuralDHS(obsI), latentFieldRurDHS(thisIndex), false);
        thislik += thisWeight*liksRurDHS(obsI, intI);
      }
      else {
        liksRurDHS(obsI, intI) = thisWeight;
      }
      
    } // for( intI )
    
    if(thislik > 0.0) {
      jnll -= log(thislik);
    }
  } // for( obsI )
  
  
  // ~~~~~~~~~~~
  // ADREPORT: used to return estimates and cov for transforms?
  // ~~~~~~~~~~~
  if(options==1){
    ADREPORT(log_tau);
    ADREPORT(logit_phi);
  }
  
  REPORT(log_tau);
  REPORT(logDet);
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
  REPORT(latentFieldUrbDHS);
  REPORT(latentFieldRurDHS);
  REPORT(fe_iUrbanDHS);
  REPORT(fe_iRuralDHS);
  
  return jnll;
  
}
