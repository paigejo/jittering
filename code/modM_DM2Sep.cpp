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
  DATA_VECTOR( y_iUrbanMICS );   // obs per binomial experiment at point i (clust)
  DATA_VECTOR( y_iRuralMICS );
  DATA_VECTOR( n_iUrbanMICS );   // Trials per cluster
  DATA_VECTOR( n_iRuralMICS );
  DATA_IVECTOR(areaidxlocUrbanMICS); // vector of length n_iUrbanMICS of areal indices
  DATA_IVECTOR(areaidxlocRuralMICS); // vector of length n_iRuralMICS of areal indices
  DATA_MATRIX( X_betaUrbanMICS );  // (nObsUrban * nIntegrationPointsUrban) x nPar design matrix. Indexed mod numObsUrban
  DATA_MATRIX( X_betaRuralMICS );  // first nObsRural rows correspond to first int pt
  DATA_ARRAY( wUrbanMICS ); // nObsUrban x nIntegrationPointsUrban weight matrix
  DATA_ARRAY( wRuralMICS ); // nObsRural x nIntegrationPointsRural weight matrix
  
  DATA_VECTOR( y_iUrbanDHS );   // obs per binomial experiment at point i (clust)
  DATA_VECTOR( y_iRuralDHS );
  DATA_VECTOR( n_iUrbanDHS );   // Trials per cluster
  DATA_VECTOR( n_iRuralDHS );
  DATA_IVECTOR(areaidxlocUrbanDHS); // vector of length n_iUrbanDHS of areal indices
  DATA_IVECTOR(areaidxlocRuralDHS); // vector of length n_iRuralDHS of areal indices
  DATA_MATRIX( X_betaUrbanDHS );  // (nObsUrban * nIntegrationPointsUrban) x nPar design matrix. Indexed mod numObsUrban
  DATA_MATRIX( X_betaRuralDHS );  // first nObsRural rows correspond to first int pt
  DATA_ARRAY( wUrbanDHS ); // nObsUrban x nIntegrationPointsUrban weight matrix
  DATA_ARRAY( wRuralDHS ); // nObsRural x nIntegrationPointsRural weight matrix
  
  // DATA_MATRIX( V_bym2 ); // Q_bym2 = V Lamda V^T is the eigendecomposition of Q_bym2
  DATA_SPARSE_MATRIX( Q_bym2 ); // really, the precision matrix for the structured, Besag component u
  
  // prior parameters
  DATA_VECTOR( alpha_pri );
  DATA_VECTOR( beta_pri );
  
  // BYM2 and prior precomputed values
  DATA_SCALAR( tr );
  DATA_VECTOR( gammaTildesm1 );
  
  // vector of row sums of Q^+, the generalized inv of Q, normalized by sum of 
  // all elements of Q^+
  DATA_VECTOR( QinvSumsNorm ); 
  
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
  
  PARAMETER( alpha ); // Intercept
  PARAMETER_VECTOR( beta ); // fixed effect/covariate effect sizes
  
  // BYM2 effects: structured (u) + unstructured (v)
  // NOTE: this is parameterized as (w=(1/sqrt(tau))(sqrt(phi) u + sqrt(1-phi) v), u)
  // these are labelled "Star" because the constraint has not yet been put on them
  PARAMETER_VECTOR( w_bym2Star );
  PARAMETER_VECTOR( u_bym2Star );
  PARAMETER_VECTOR( nuggetUrbMICS );
  PARAMETER_VECTOR( nuggetRurMICS );
  PARAMETER_VECTOR( nuggetUrbDHS );
  PARAMETER_VECTOR( nuggetRurDHS );
  
  // ~~~~~~~~~------------------------------------------------~~
  // SECOND, we define all other objects that we need internally
  // ~~~~~~~~~------------------------------------------------~~
  
  // calculated values from data
  int num_iUrbanMICS = y_iUrbanMICS.size();   // Number of urban data points in space
  int num_iRuralMICS = y_iRuralMICS.size();   // Number of rural data points in space
  int num_iUrbanDHS = y_iUrbanDHS.size();   // Number of urban data points in space
  int num_iRuralDHS = y_iRuralDHS.size();   // Number of rural data points in space
  int nAreas = u_bym2Star.size();               // Number of areas
  int n_integrationPointsUrbanMICS = wUrbanMICS.cols();    // number of integration points for each observation
  int n_integrationPointsRuralMICS = wRuralMICS.cols();
  int n_integrationPointsUrbanDHS = wUrbanDHS.cols();    // number of integration points for each observation
  int n_integrationPointsRuralDHS = wRuralDHS.cols();
  
  // Transform some of our parameters
  // Type sqrtTau = exp(Type(0.5) * log_tau);
  Type tau = exp(log_tau);
  Type phi = Type(1.0)/(Type(1.0) + exp(-logit_phi));
  // Type phi = exp(logit_phi)/(Type(1.0) + exp(logit_phi));
  // Type logPhi = logit_phi - log((Type(1.0) + exp(logit_phi)));
  // Type logPhi = -log((Type(1.0) + exp(-logit_phi)));
  // oneMphi = exp(-logit_phi)/(Type(1.0) + exp(-logit_phi));
  // oneMphi = Type(1.0)/(Type(1.0) + exp(logit_phi));
  // Type logOneMphi = -log(Type(1.0) + exp(logit_phi));
  Type sigmaEps = exp(Type(-0.5) * log_tauEps);
  
  // make sure u sums to 0 (not w, the full BYM2 effect, according to Andrea), 
  // and then subtract off scaled u from w to get v
  vector<Type> w_bym2 = w_bym2Star;
  // vector<Type> v_bym2 = w_bym2Star;
  // vector<Type> u_bym2 = u_bym2Star;
  Type uSum = u_bym2Star.sum();
  Type scaleUpU = sqrt(phi/tau);
  Type uFac = scaleUpU * uSum;
  Type reduceScaledUpUi = 0.0;
  // Type scaleDownV = sqrt( tau/(Type(1.0)-phi) );
  for(int i = 0; i < nAreas; i++) {
    // ensure unit scaled u sums to 0:
    reduceScaledUpUi = uFac * QinvSumsNorm(i);
    
    // u_bym2(i) = scaleUpU * (u_bym2(i) - uSum * QinvSumsNorm(i));
    // u_bym2(i) = scaleUpU * u_bym2(i) - reduceScaledUpUi;
    
    // no need to calculate constrained v to get w, since v not affected by constraint
    w_bym2(i) -= reduceScaledUpUi;
    // v_bym2(i) = scaleDownV * (w_bym2(i) - scaleUpU * u_bym2(i));
    
    // v not affected by constraint, but must be scaled appropriately:
    // v_bym2(i) = w_bym2Star(i) - scaleUpU * u_bym2Star(i);
  }
  
  // objective function -- joint negative log-likelihood/posterior
  Type jnll = 0.0;
  
  vector<Type> projepsilon_iUrbanDHS(num_iUrbanDHS);
  vector<Type> projepsilon_iRuralDHS(num_iRuralDHS);
  
  // Project GP approx from mesh points to data points
  // projepsilon_iUrbanDHS = (AprojUrbanDHS * Epsilon_bym2.matrix()).col(0);
  // projepsilon_iRuralDHS = (AprojRuralDHS * Epsilon_bym2.matrix()).col(0);
  
  // ~~~~~~~~~------------------------------------------------~~-
  // THIRD, we calculate the contribution to the likelihood from:
  // 1) GMRF first, for 'flag' normalization purposes
  // 2) priors
  // 3) likelihood
  // ~~~~~~~~~------------------------------------------------~~-
  
  /////////
  // (1) //
  /////////
  
  // add in BYM2 latent model component to the posterior
  
  // Calculate the quadratic form Eps tau [(1-phi) I + phi Q_besag^+]^(-1) Eps^T
  // = Eps tau [I + phi (Q_besag^+ - I)]^(-1) Eps^T
  // = Eps tau V [I + phi (GammaTilde - I)]^(-1) V^T Eps^T
  // i.e. the sum of squares of tau^0.5 Eps V diag(1/sqrt(1 + phi*gammaTildesm1))
  // matrix<Type> transformedEpsilon = Epsilon_bym2.matrix().transpose() * V_bym2;
  // Type quad = 0.0;
  // for(int i = 0; i < transformedEpsilon.size(); i++) {
  //   quad += pow(transformedEpsilon(0,i) * sqrtTau / sqrt(1.0 + phi*gammaTildesm1(i)), 2);
  // }
  
  // Calculate the quadratic form: 
  // x^T Q_x x
  //   = tau/(1-phi) w*^T w* + u*^T Q u* + phi/(1-phi) u*^T u* - 2 sqrt(phi tau)/(1-phi) u*^T w*
  // where the * denotes the unconstrained effect
  
  // Calculate the quadratic form tau/(1-phi) w*^T w*
  // where the * denotes the unconstrained effect
  Type quadW = 0.0;
  for(int i = 0; i < nAreas; i++) {
    quadW += w_bym2Star(i) * w_bym2Star(i);
  }
  quadW = quadW * tau/(Type(1.0)-phi);
  
  // Calculate the quadratic form u*^T Q u* + phi/(1-phi) u*^T u*
  // where the * denotes the unconstrained effect
  matrix<Type> transformedU = Q_bym2 * u_bym2Star.matrix();
  Type quadU = 0.0;
  Type fac = phi/(Type(1.0)-phi);
  for(int i = 0; i < nAreas; i++) {
    quadU += transformedU(i) * u_bym2Star(i) + fac * u_bym2Star(i) * u_bym2Star(i);
  }
  
  // Calculate the quadratic form -2 sqrt(phi tau)/(1-phi) u*^T w*
  // where the * denotes the unconstrained effect
  Type quadWU = 0.0;
  for(int i = 0; i < nAreas; i++) {
    quadWU += u_bym2Star(i) * w_bym2Star(i);
  }
  quadWU = quadWU * (-Type(2.0) * sqrt(phi * tau)/(Type(1.0)-phi));
  
  // // Calculate the quadratic form EpsV^T tau [(1-phi) I]^(-1) EpsV
  // // = (tau / (1-phi)) EpsV^T EpsV
  // Type quadV = 0.0;
  // for(int i = 0; i < nAreas; i++) {
  //   quadV += pow(v_bym2(i), 2.0);
  // }
  // quadV = quadV * (tau / (Type(1.0) - phi));
  
  Type quadSum = quadW + quadU + quadWU;
  
  // calculate log determinant of (1/tau) [(1-phi) I + phi Q_besag^+]
  // (used for PC prior for phi and old version of BYM2 density)
  Type logDet = Type(0.0);
  for(int i = 0; i < gammaTildesm1.size(); i++) {
    logDet += log(Type(1.0) + phi*gammaTildesm1(i));
  }
  // Type logDetTau = logDet + gammaTildesm1.size() * (-log_tau);
  
  // leave out the log determinant of the Q term, since it is constant, 
  // but include the term for the parameters. In other words, calculate part of: 
  // log |(1/tau) * phi * Q^+| + log|(1/tau) * (1-phi) * I| 
  //   = nAreas * log(phi / tau) + log|Q^+| + nAreas * log((1-phi) / tau)
  //   = nAreas * log(phi (1 - phi) / tau^2) + log|Q^+|
  //   = nAreas * log(phi (1 - phi) / tau^2) + C
  // Type logDetTau = Type(nAreas) * log(phi * (1.0 - phi) / pow(tau, 2));
  // Type logDetTau = Type(nAreas) * (log(phi) + log(Type(1.0) - phi) - Type(2.0) * log_tau);
  
  // for x = (w^T u^T)^T, 
  // |Q_x^-1| = |AD - BC| (from wikipedia on block matrices)
  //   = ...
  //   = |tau/(1-phi) Q|^-1
  //   = ...
  //   = (1 - phi)^n / tau^n |Q^+|
  //   propto (1 - phi)^n / tau^n
  Type logDetTau = Type(nAreas) * log((Type(1.0) - phi)/tau); // leave out the constant
  
  // add quadratic form and normalizing constant to the negative posterior log density
  // first add denominator
  // Type bym2LogLik = -0.5 * (2.0 * Type(nAreas) * log(2.0*M_PI) + logDetTau);
  // Type bym2LogLik = Type(-0.5) * (Type(2.0) * Type(nAreas) * (log(2.0) + log(M_PI)) + logDetTau);
  Type bym2LogLik = Type(-0.5) * logDetTau; // leave out the constant
  
  // now add exponent
  bym2LogLik += Type(-0.5) * quadSum;
  // add all to the negative posterior log density
  jnll -= bym2LogLik;
  
  // add in iid nugget model
  Type nuggetLogLik = 0.0;
  for(int i = 0; i< nuggetUrbMICS.size(); i++) {
    nuggetLogLik += dnorm(nuggetUrbMICS(i), Type(0), sigmaEps, true); // N(mean, sd);
  }
  for(int i = 0; i< nuggetRurMICS.size(); i++) {
    nuggetLogLik += dnorm(nuggetRurMICS(i), Type(0), sigmaEps, true); // N(mean, sd);
  }
  for(int i = 0; i< nuggetUrbDHS.size(); i++) {
    nuggetLogLik += dnorm(nuggetUrbDHS(i), Type(0), sigmaEps, true); // N(mean, sd);
  }
  for(int i = 0; i< nuggetRurDHS.size(); i++) {
    nuggetLogLik += dnorm(nuggetRurDHS(i), Type(0), sigmaEps, true); // N(mean, sd);
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
  
  // Type tau = exp(log_tau);
  
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
  Type KLD = Type(0.5) * (phi * tr - phi * n - logDet);
  Type d = sqrt(Type(2.0) * KLD);
  
  // calculate exponential density of d(phi)
  Type lexpDensity = log(lambdaPhi) - lambdaPhi * d;
  
  // add Jacobian with respect to d(phi) and logit(phi)
  Type sumVal = 0.0;
  for(int i = 0; i < n; i++) {
    sumVal += gammaTildesm1(i)/(1 + gammaTildesm1(i) * phi);
  }
  //Type ljacobian = - log(d) - log(2) + log(abs(tr - n - sumVal));
  Type ljacobian = - log(d) - log(2.0) + log(tr - n - sumVal); // this should be fine since argument should be positive
  
  ljacobian = ljacobian - logit_phi - Type(2.0) * log(Type(1.0) + exp(-logit_phi));
  
  Type logPriPhi = lexpDensity + ljacobian;
  
  // get final log density
  jnll -= logPriTau;
  jnll -= logPriTauEps;
  jnll -= logPriPhi;
  
  // prior for intercept
  // Type alphaPri = dnorm(alpha, alpha_pri(0), alpha_pri(1), true); // N(mean, sd)
  // jnll -= alphaPri;
  
  // prior for other covariates
  for(int i = 0; i < beta.size(); i++) {
    jnll -= dnorm(beta(i), beta_pri(0), beta_pri(1), true); // N(mean, sd)
  }
  
  /////////
  // (3) //
  /////////
  // jnll contribution from each datapoint i
  
  // Calculate fixed effects at data points
  vector<Type> fe_iUrbanMICS = X_betaUrbanMICS * beta;
  vector<Type> fe_iRuralMICS = X_betaRuralMICS * beta;
  vector<Type> fe_iUrbanDHS = X_betaUrbanDHS * beta;
  vector<Type> fe_iRuralDHS = X_betaRuralDHS * beta;
  
  vector<Type> latentFieldUrbMICS(num_iUrbanMICS * n_integrationPointsUrbanMICS);
  vector<Type> latentFieldRurMICS(num_iRuralMICS * n_integrationPointsRuralMICS);
  vector<Type> latentFieldUrbDHS(num_iUrbanDHS * n_integrationPointsUrbanDHS);
  vector<Type> latentFieldRurDHS(num_iRuralDHS * n_integrationPointsRuralDHS);
  int thisIndex;
  // Type thisLatentField;
  Type thisWeight;
  Type thislik;
  array<Type> liksUrbMICS(num_iUrbanMICS, n_integrationPointsUrbanMICS);
  array<Type> liksRurMICS(num_iRuralMICS, n_integrationPointsRuralMICS);
  array<Type> liksUrbDHS(num_iUrbanDHS, n_integrationPointsUrbanDHS);
  array<Type> liksRurDHS(num_iRuralDHS, n_integrationPointsRuralDHS);
  // Type EpsInd = 0;
  
  for (int obsI = 0; obsI < num_iUrbanMICS; obsI++) {
    thislik = 0.0;
    
    for (int intI = 0; intI < n_integrationPointsUrbanMICS; intI++) {
      thisIndex = num_iUrbanMICS * intI + obsI;
      
      // latent field estimate at each obs
      // thisLatentField = fe_iUrbanMICS(thisIndex) + projepsilon_iUrbanMICS(obsI);
      latentFieldUrbMICS(thisIndex) = alpha + fe_iUrbanMICS(thisIndex) + w_bym2(areaidxlocUrbanMICS(thisIndex)) + nuggetUrbMICS(obsI);
      // latentFieldUrbMICS(thisIndex) = fe_iUrbanMICS(thisIndex) + projepsilon_iUrbanMICS(obsI);
      
      // and add data contribution to jnll
      // get integration weight
      thisWeight = wUrbanMICS(obsI,intI);
      
      // Uses the dbinom_robust function, which takes the logit probability
      if(thisWeight > 0.0) {
        // thislik += thisWeight*dbinom_robust( y_iUrbanMICS(obsI), n_iUrbanMICS(obsI), thisLatentField, false);
        liksUrbMICS(obsI, intI) = dbinom_robust( y_iUrbanMICS(obsI), n_iUrbanMICS(obsI), latentFieldUrbMICS(thisIndex), false);
        thislik += thisWeight*liksUrbMICS(obsI, intI);
      }
      else {
        liksUrbMICS(obsI, intI) = thisWeight;
      }
      
    } // for( intI )
    
    if(thislik > 0.0) {
      jnll -= log(thislik);
    }
  } // for( obsI )
  
  for (int obsI = 0; obsI < num_iRuralMICS; obsI++) {
    thislik = 0.0;
    
    for (int intI = 0; intI < n_integrationPointsRuralMICS; intI++) {
      thisIndex = num_iRuralMICS * intI + obsI;
      
      // latent field estimate at each obs
      // thisLatentField = fe_iRuralMICS(thisIndex) + projepsilon_iRuralMICS(obsI);
      latentFieldRurMICS(thisIndex) = alpha + fe_iRuralMICS(thisIndex) + w_bym2(areaidxlocRuralMICS(thisIndex)) + nuggetRurMICS(obsI);
      // latentFieldRurMICS(thisIndex) = fe_iRuralMICS(thisIndex) + projepsilon_iRuralMICS(obsI);
      
      // and add data contribution to jnll
      // get integration weight
      thisWeight = wRuralMICS(obsI,intI);
      
      // Uses the dbinom_robust function, which takes the logit probability
      if(thisWeight > 0.0) {
        // thislik += thisWeight*dbinom_robust( y_iRuralMICS(obsI), n_iRuralMICS(obsI), thisLatentField, false);
        liksRurMICS(obsI, intI) = dbinom_robust( y_iRuralMICS(obsI), n_iRuralMICS(obsI), latentFieldRurMICS(thisIndex), false);
        thislik += thisWeight*liksRurMICS(obsI, intI);
      }
      else {
        liksRurMICS(obsI, intI) = thisWeight;
      }
      
    } // for( intI )
    
    if(thislik > 0.0) {
      jnll -= log(thislik);
    }
  } // for( obsI )
  
  for (int obsI = 0; obsI < num_iUrbanDHS; obsI++) {
    thislik = 0.0;
    
    for (int intI = 0; intI < n_integrationPointsUrbanDHS; intI++) {
      thisIndex = num_iUrbanDHS * intI + obsI;
      
      // latent field estimate at each obs
      // thisLatentField = fe_iUrbanDHS(thisIndex) + projepsilon_iUrbanDHS(obsI);
      latentFieldUrbDHS(thisIndex) = alpha + fe_iUrbanDHS(thisIndex) + w_bym2(areaidxlocUrbanDHS(thisIndex)) + nuggetUrbDHS(obsI);
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
      latentFieldRurDHS(thisIndex) = alpha + fe_iRuralDHS(thisIndex) + w_bym2(areaidxlocRuralDHS(thisIndex)) + nuggetRurDHS(obsI);
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
    ADREPORT(w_bym2);
  }
  
  REPORT(w_bym2);
  REPORT(log_tau);
  // REPORT(alphaPri);
  REPORT(quadU);
  REPORT(quadW);
  REPORT(quadWU);
  REPORT(quadSum);
  // REPORT(logDet);
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
  REPORT(liksUrbDHS);
  REPORT(liksRurDHS);
  // REPORT(transformedEpsilon);
  REPORT(latentFieldUrbMICS);
  REPORT(latentFieldRurMICS);
  REPORT(latentFieldUrbDHS);
  REPORT(latentFieldRurDHS);
  REPORT(fe_iUrbanMICS);
  REPORT(fe_iRuralMICS);
  REPORT(fe_iUrbanDHS);
  REPORT(fe_iRuralDHS);
  
  
  return jnll;
  
}
