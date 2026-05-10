// Constrained 2n-2 BYM2 model for M_DM (DHS+MICS fusion) data
// Uses n-1 free params for w and n-1 free params for u,
// deriving the nth element as -sum(free) to enforce sum-to-zero constraints.
// Explicit alpha intercept parameter.
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
  
  DATA_SPARSE_MATRIX( Q_bym2 ); // precision matrix for Besag component
  
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
  
  // Fixed effects
  PARAMETER( log_tau );
  PARAMETER( logit_phi );
  PARAMETER( log_tauEps );
  
  // Intercept (explicit)
  PARAMETER( alpha );
  PARAMETER_VECTOR( beta );
  
  // BYM2 effects: n-1 free parameters each
  PARAMETER_VECTOR( w_bym2Free ); // length n-1
  PARAMETER_VECTOR( u_bym2Free ); // length n-1
  PARAMETER_VECTOR( nuggetUrbMICS );
  PARAMETER_VECTOR( nuggetRurMICS );
  PARAMETER_VECTOR( nuggetUrbDHS );
  PARAMETER_VECTOR( nuggetRurDHS );
  
  // ~~~~~~~~~------------------------------------------------~~
  // SECOND, we define all other objects that we need internally
  // ~~~~~~~~~------------------------------------------------~~
  
  int num_iUrbanMICS = y_iUrbanMICS.size();
  int num_iRuralMICS = y_iRuralMICS.size();
  int num_iUrbanDHS = y_iUrbanDHS.size();
  int num_iRuralDHS = y_iRuralDHS.size();
  int nFree = w_bym2Free.size();
  int nAreas = nFree + 1;
  int n_integrationPointsUrbanMICS = wUrbanMICS.cols();
  int n_integrationPointsRuralMICS = wRuralMICS.cols();
  int n_integrationPointsUrbanDHS = wUrbanDHS.cols();
  int n_integrationPointsRuralDHS = wRuralDHS.cols();
  
  // Transform parameters
  Type tau = exp(log_tau);
  Type phi = Type(1.0)/(Type(1.0) + exp(-logit_phi));
  Type sigmaEps = exp(Type(-0.5) * log_tauEps);
  
  // Reconstruct full n-vectors with sum-to-zero constraint
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
  
  // ~~~~~~~~~------------------------------------------------~~-
  // THIRD, we calculate the contribution to the likelihood from:
  // 1) BYM2 GMRF density
  // 2) priors
  // 3) likelihood
  // ~~~~~~~~~------------------------------------------------~~-
  
  /////////
  // (1) //
  /////////
  
  // Quadratic form using full n-vectors
  
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
  
  // Log determinant for the constrained model: (n-1)*log((1-phi)/tau)
  Type logDet = Type(0.0);
  for(int i = 0; i < gammaTildesm1.size(); i++) {
    logDet += log(Type(1.0) + phi*gammaTildesm1(i));
  }
  Type logDetTau = Type(nFree) * log((Type(1.0) - phi)/tau);
  
  Type bym2LogLik = Type(-0.5) * logDetTau;
  bym2LogLik += Type(-0.5) * quadSum;
  jnll -= bym2LogLik;
  
  // iid nugget model
  Type nuggetLogLik = 0.0;
  for(int i = 0; i< nuggetUrbMICS.size(); i++) {
    nuggetLogLik += dnorm(nuggetUrbMICS(i), Type(0), sigmaEps, true);
  }
  for(int i = 0; i< nuggetRurMICS.size(); i++) {
    nuggetLogLik += dnorm(nuggetRurMICS(i), Type(0), sigmaEps, true);
  }
  for(int i = 0; i< nuggetUrbDHS.size(); i++) {
    nuggetLogLik += dnorm(nuggetUrbDHS(i), Type(0), sigmaEps, true);
  }
  for(int i = 0; i< nuggetRurDHS.size(); i++) {
    nuggetLogLik += dnorm(nuggetRurDHS(i), Type(0), sigmaEps, true);
  }
  jnll -= nuggetLogLik;
  
  /////////
  // (2) //
  /////////
  
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
  
  // Fixed effects at data points
  vector<Type> fe_iUrbanMICS = X_betaUrbanMICS * beta;
  vector<Type> fe_iRuralMICS = X_betaRuralMICS * beta;
  vector<Type> fe_iUrbanDHS = X_betaUrbanDHS * beta;
  vector<Type> fe_iRuralDHS = X_betaRuralDHS * beta;
  
  vector<Type> latentFieldUrbMICS(num_iUrbanMICS * n_integrationPointsUrbanMICS);
  vector<Type> latentFieldRurMICS(num_iRuralMICS * n_integrationPointsRuralMICS);
  vector<Type> latentFieldUrbDHS(num_iUrbanDHS * n_integrationPointsUrbanDHS);
  vector<Type> latentFieldRurDHS(num_iRuralDHS * n_integrationPointsRuralDHS);
  int thisIndex;
  Type thisWeight;
  Type thislik;
  array<Type> liksUrbMICS(num_iUrbanMICS, n_integrationPointsUrbanMICS);
  array<Type> liksRurMICS(num_iRuralMICS, n_integrationPointsRuralMICS);
  array<Type> liksUrbDHS(num_iUrbanDHS, n_integrationPointsUrbanDHS);
  array<Type> liksRurDHS(num_iRuralDHS, n_integrationPointsRuralDHS);
  
  // MICS Urban
  for (int obsI = 0; obsI < num_iUrbanMICS; obsI++) {
    thislik = 0.0;
    for (int intI = 0; intI < n_integrationPointsUrbanMICS; intI++) {
      thisIndex = num_iUrbanMICS * intI + obsI;
      latentFieldUrbMICS(thisIndex) = alpha + fe_iUrbanMICS(thisIndex) + w_bym2(areaidxlocUrbanMICS(obsI)) + nuggetUrbMICS(obsI);
      thisWeight = wUrbanMICS(obsI,intI);
      if(thisWeight > 0.0) {
        liksUrbMICS(obsI, intI) = dbinom_robust( y_iUrbanMICS(obsI), n_iUrbanMICS(obsI), latentFieldUrbMICS(thisIndex), false);
        thislik += thisWeight*liksUrbMICS(obsI, intI);
      } else {
        liksUrbMICS(obsI, intI) = thisWeight;
      }
    }
    if(thislik > 0.0) {
      jnll -= log(thislik);
    }
  }
  
  // MICS Rural
  for (int obsI = 0; obsI < num_iRuralMICS; obsI++) {
    thislik = 0.0;
    for (int intI = 0; intI < n_integrationPointsRuralMICS; intI++) {
      thisIndex = num_iRuralMICS * intI + obsI;
      latentFieldRurMICS(thisIndex) = alpha + fe_iRuralMICS(thisIndex) + w_bym2(areaidxlocRuralMICS(obsI)) + nuggetRurMICS(obsI);
      thisWeight = wRuralMICS(obsI,intI);
      if(thisWeight > 0.0) {
        liksRurMICS(obsI, intI) = dbinom_robust( y_iRuralMICS(obsI), n_iRuralMICS(obsI), latentFieldRurMICS(thisIndex), false);
        thislik += thisWeight*liksRurMICS(obsI, intI);
      } else {
        liksRurMICS(obsI, intI) = thisWeight;
      }
    }
    if(thislik > 0.0) {
      jnll -= log(thislik);
    }
  }
  
  // DHS Urban
  for (int obsI = 0; obsI < num_iUrbanDHS; obsI++) {
    thislik = 0.0;
    for (int intI = 0; intI < n_integrationPointsUrbanDHS; intI++) {
      thisIndex = num_iUrbanDHS * intI + obsI;
      latentFieldUrbDHS(thisIndex) = alpha + fe_iUrbanDHS(thisIndex) + w_bym2(areaidxlocUrbanDHS(obsI)) + nuggetUrbDHS(obsI);
      thisWeight = wUrbanDHS(obsI,intI);
      if(thisWeight > 0.0) {
        liksUrbDHS(obsI, intI) = dbinom_robust( y_iUrbanDHS(obsI), n_iUrbanDHS(obsI), latentFieldUrbDHS(thisIndex), false);
        thislik += thisWeight*liksUrbDHS(obsI, intI);
      } else {
        liksUrbDHS(obsI, intI) = thisWeight;
      }
    }
    if(thislik > 0.0) {
      jnll -= log(thislik);
    }
  }
  
  // DHS Rural
  for (int obsI = 0; obsI < num_iRuralDHS; obsI++) {
    thislik = 0.0;
    for (int intI = 0; intI < n_integrationPointsRuralDHS; intI++) {
      thisIndex = num_iRuralDHS * intI + obsI;
      latentFieldRurDHS(thisIndex) = alpha + fe_iRuralDHS(thisIndex) + w_bym2(areaidxlocRuralDHS(obsI)) + nuggetRurDHS(obsI);
      thisWeight = wRuralDHS(obsI,intI);
      if(thisWeight > 0.0) {
        liksRurDHS(obsI, intI) = dbinom_robust( y_iRuralDHS(obsI), n_iRuralDHS(obsI), latentFieldRurDHS(thisIndex), false);
        thislik += thisWeight*liksRurDHS(obsI, intI);
      } else {
        liksRurDHS(obsI, intI) = thisWeight;
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
  
  REPORT(w_bym2);
  REPORT(u_bym2);
  REPORT(log_tau);
  REPORT(quadU);
  REPORT(quadW);
  REPORT(quadWU);
  REPORT(quadSum);
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
