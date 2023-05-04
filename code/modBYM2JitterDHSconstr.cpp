// include libraries
#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
using namespace density;
using Eigen::SparseMatrix;

// PC prior on BYM2 phi parameter (on logit scale)
template<class Type>
Type dPCPriPhi(Type logitPhi, Type lambda, Type tr, vector<Type> gammaTildesm1, Type logDet,
               int give_log=1)
{
  Type phi = (1 / (1 + exp(-logitPhi)));
  
  // calculate log determinant, KLD(phi), and d(phi)
  // Type logDet = sum(log(1 + phi*gammaTildesm1));
  
  Type n = gammaTildesm1.size();
  Type KLD = 0.5 * (phi * tr - phi * n - logDet);
  Type d = sqrt(2*KLD);
  
  // calculate exponential density of d(phi)
  Type lexpDensity = log(lambda) - lambda * d;
  
  // add Jacobian with respect to d(phi) and logit(phi)
  Type sumVal = 0;
  for(int i = 0; i < n; i++) {
    sumVal += gammaTildesm1(i)/(1 + gammaTildesm1(i) * phi);
  }
  //Type ljacobian = - log(d) - log(2) + log(abs(tr - n - sumVal));
  Type ljacobian = - log(d) - log(2) + log(tr - n - sumVal); // this should be fine since argument should be positive
  
  ljacobian = ljacobian - logitPhi - 2*log(1 + exp(-logitPhi));
  
  // get final log density
  Type lPriorPhi = lexpDensity + ljacobian;
  
  if(give_log) {
    return lPriorPhi;
  }
  else {
    return exp(lPriorPhi);
  }
}

// PC prior on BYM2 tau (precision) parameter (on log scale)
template<class Type>
Type dPCPriTau(Type logTau, Type lambda)
{
  Type sqrtTau = exp(0.5 * logTau);
  
  Type ldensity = log(lambda) - log(2) -(3/2)*logTau - lambda/sqrtTau;
  
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
  DATA_MATRIX( AprojUrbanDHS ); // nObsUrban x nArea matrix with ij-th entry = 1 if cluster i associated with area j and 0 o.w.
  DATA_MATRIX( AprojRuralDHS ); // nObsRural x nArea matrix with ij-th entry = 1 if cluster i associated with area j and 0 o.w.
  DATA_MATRIX( X_betaUrbanDHS );  // (nObsUrban * nIntegrationPointsUrban) x nPar design matrix. Indexed mod numObsUrban
  DATA_MATRIX( X_betaRuralDHS );  // first nObsRural rows correspond to first int pt
  DATA_ARRAY( wUrbanDHS ); // nObsUrban x nIntegrationPointsUrban weight matrix
  DATA_ARRAY( wRuralDHS ); // nObsRural x nIntegrationPointsRural weight matrix
  
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
  
  // BYM2 effects (structured + unstructured), but minus 1 effect since there is
  // a constraint NOTE: may need to modify parameterization later
  PARAMETER_VECTOR( Epsilon_bym2 );
  
  // ~~~~~~~~~------------------------------------------------~~
  // SECOND, we define all other objects that we need internally
  // ~~~~~~~~~------------------------------------------------~~
  
  // calculated values from data
  int num_iUrbanDHS = y_iUrbanDHS.size();   // Number of urban data points in space
  int num_iRuralDHS = y_iRuralDHS.size();   // Number of rural data points in space
  int n_integrationPointsUrbanDHS = wUrbanDHS.cols();    // number of integration points for each observation
  int n_integrationPointsRuralDHS = wRuralDHS.cols();
  int n = gammaTildesm1.size();
  
  // objective function -- joint negative log-likelihood/posterior
  Type jnll = 0;
  
  // Transform some of our parameters
  Type sqrtTau = exp(0.5 * log_tau);
  Type phi = 1/(1 + exp(-logit_phi));
  
  // normalize spatial field using precision tau
  // vector<Type> unscaledEpsilon(Epsilon_bym2.size());
  // for(int i = 0; i < unscaledEpsilon.size(); i++) {
  //   unscaledEpsilon(i) = Epsilon_bym2(i) * sqrtTau;
  // }
  
  // Define objects for derived values
  vector<Type> fe_iUrbanDHS(num_iUrbanDHS * n_integrationPointsUrbanDHS); // main effect: alpha
  vector<Type> fe_iRuralDHS(num_iRuralDHS * n_integrationPointsRuralDHS);
  // Logit estimated prob for each cluster i
  //vector<Type> latent_field_iUrban(num_iUrban);
  //vector<Type> latent_field_iRural(num_iRural);
  // value of gmrf at data points
  vector<Type> projepsilon_iUrbanDHS(num_iUrbanDHS);
  vector<Type> projepsilon_iRuralDHS(num_iRuralDHS);
  
  // linear combination of fixed effects
  // fe_iUrbanDHS = X_betaUrbanDHS * beta + alpha;
  // fe_iRuralDHS = X_betaRuralDHS * beta + alpha;
  fe_iUrbanDHS = X_betaUrbanDHS * beta.matrix();
  fe_iRuralDHS = X_betaRuralDHS * beta.matrix();
  for(int i = 0; i < fe_iUrbanDHS.size(); i++) {
    fe_iUrbanDHS(i) = fe_iUrbanDHS(i) + alpha;
  }
  for(int i = 0; i < fe_iRuralDHS.size(); i++) {
    fe_iRuralDHS(i) = fe_iRuralDHS(i) + alpha;
  }
  
  // transform Epsilon via V* diag(sqrt(gammaTildem1* + 1)) Epsilon
  // Where V* is V but with constant eigenvector removed and gammaTildem1* + 1 are 
  // the eigenvalues corresponding to the eigenvectors of V*
  vector<Type> scales(n-1);
  for(int i = 0; i < n-1; i++) {
    scales(i) = sqrt(gammaTildem1(i) + 1);
  }
  matrix<Type> transformedEpsilon = V_bym2 * (scales * Epsilon_bym2).matrix();
  
  // Project GP approx from mesh points to data points
  // projepsilon_iUrban = AprojUrban * Epsilon_s.matrix();
  // projepsilon_iRural = AprojRural * Epsilon_s.matrix();
  projepsilon_iUrbanDHS = AprojUrbanDHS * Epsilon_bym2.matrix();
  projepsilon_iRuralDHS = AprojRuralDHS * Epsilon_bym2.matrix();
  
  // ~~~~~~~~~------------------------------------------------~~-
  // THIRD, we calculate the contribution to the likelihood from:
  // 1) GP field first, for 'flag' normalization purposes
  // 2) priors
  // 3) GP field
  // ~~~~~~~~~------------------------------------------------~~-
  
  /////////
  // (1) //
  /////////
  
  // add in BYM2 latent model component to the posterior (the log density without the normalizing constant)
  // jnll += GMRF(Q_bym2, false)(unscaledEpsilon);
  
  // Calculate the quadratic form Eps tau [(1-phi) I + phi Q_besag]^(-1) Eps^T
  // Assume transformedEpsilon has been premultiplied by V^T
  Type quad = 0;
  for(int i = 0; i < transformedEpsilon.size(); i++) {
    quad += pow(transformedEpsilon(0,i) * sqrtTau / (1 + phi*gammaTildesm1(i)), 2);
  }
  
  // calculate log determinant of (1/tau) [(1-phi) I + phi Q_besag^(-1)]
  Type logDet = Type(0);
  for(int i = 0; i < gammaTildesm1.size(); i++) {
    logDet += log(1 + phi*gammaTildesm1(i));
  }
  Type logDetTau = logDet + gammaTildesm1.size() * (-log_tau);
  
  // add quadratic form and normalizing constant to the negative posterior log density
  // first add denominator
  Type bym2LogLik = -((0.5 * gammaTildesm1.size()) * (log(2) + log(M_PI)) + 0.5 * logDetTau);
  // now add exponent
  bym2LogLik += -0.5 * quad;
  // add all to the negative posterior log density
  jnll -= bym2LogLik;
  
  /////////
  // (2) //
  /////////
  // Prior contributions to joint likelihood
  
  // add in priors for BYM2 parameter
  Type logPriTau = dPCPriTau(log_tau, lambdaTau);
  // Type logPriPhi = dPCPriPhi(logit_phi, lambdaPhi, tr, gammaTildesm1, logDet, true);
  
  // calculate log determinant, KLD(phi), and d(phi)
  // Type logDet = sum(log(1 + phi*gammaTildesm1));
  
  Type n = gammaTildesm1.size();
  Type KLD = 0.5 * (phi * tr - phi * n - logDet);
  Type d = sqrt(2*KLD);
  
  // calculate exponential density of d(phi)
  Type lexpDensity = log(lambdaPhi) - lambdaPhi * d;
  
  // add Jacobian with respect to d(phi) and logit(phi)
  Type sumVal = 0;
  for(int i = 0; i < n; i++) {
    sumVal += gammaTildesm1(i)/(1 + gammaTildesm1(i) * phi);
  }
  //Type ljacobian = - log(d) - log(2) + log(abs(tr - n - sumVal));
  Type ljacobian = - log(d) - log(2) + log(tr - n - sumVal); // this should be fine since argument should be positive
  
  ljacobian = ljacobian - logit_phi - 2*log(1 + exp(-logit_phi));
  
  // get final log density
  Type logPriPhi = lexpDensity + ljacobian;
  jnll -= logPriTau;
  jnll -= logPriPhi;
  
  // prior for intercept
  jnll -= dnorm(alpha, alpha_pri(0), alpha_pri(1), true); // N(mean, sd)
  
  // prior for other covariates
  for(int i = 0; i < beta.size(); i++) {
    jnll -= dnorm(beta(i), beta_pri(0), beta_pri(1), true); // N(mean, sd)
  }
  
  /////////
  // (3) //
  /////////
  // jnll contribution from each datapoint i
  int thisIndex;
  Type thisLatentField;
  Type thisWeight;
  Type thislik;
  array<Type> liksUrbDHS(num_iUrbanDHS, n_integrationPointsUrbanDHS);
  
  for (int intI = 0; intI < n_integrationPointsUrbanDHS; intI++) {
    thislik = 0;
    
    for (int obsI = 0; obsI < num_iUrbanDHS; obsI++) {
      thisIndex = num_iUrbanDHS * intI + obsI;
      
      // latent field estimate at each obs
      thisLatentField = fe_iUrbanDHS(thisIndex) + projepsilon_iUrbanDHS(obsI);
      
      // and add data contribution to jnll
      // get integration weight
      thisWeight = wUrbanDHS(obsI,intI);
      
      // Uses the dbinom_robust function, which takes the logit probability
      if(thisWeight > 0) {
        // thislik += thisWeight*dbinom_robust( y_iUrbanDHS(obsI), n_iUrbanDHS(obsI), thisLatentField, false);
        liksUrbDHS(obsI, intI) = dbinom_robust( y_iUrbanDHS(obsI), n_iUrbanDHS(obsI), thisLatentField, false);
        thislik += thisWeight*liksUrbDHS(obsI, intI);
      }
      else {
        liksUrbDHS(obsI, intI) = thisWeight;
      }
      
    } // for( intI )
    
    jnll -= log(thislik);
  } // for( obsI )
  
  for (int intI = 0; intI < n_integrationPointsRuralDHS; intI++) {
    thislik = 0;
    
    for (int obsI = 0; obsI < num_iRuralDHS; obsI++) {
      thisIndex = num_iRuralDHS * intI + obsI;
      
      // latent field estimate at each obs
      thisLatentField = fe_iRuralDHS(thisIndex) + projepsilon_iRuralDHS(obsI);
      
      // and add data contribution to jnll
      // get integration weight
      thisWeight = wRuralDHS(obsI,intI);
      
      // Uses the dbinom_robust function, which takes the logit probability
      if(thisWeight > 0) {
        // could try Kahan summation, but who knows if this would affect the gradient...
        // https://scicomp.stackexchange.com/questions/21483/how-to-avoid-the-round-off-errors-in-the-larger-calculations
        thislik += thisWeight*dbinom_robust( y_iRuralDHS(obsI), n_iRuralDHS(obsI), thisLatentField, false);
      }
      
    } // for( intI )
    
    jnll -= log(thislik);
  } // for( obsI )
  
  // ~~~~~~~~~~~
  // ADREPORT: used to return estimates and cov for transforms?
  // ~~~~~~~~~~~
  if(options==1){
    ADREPORT(log_tau);
    ADREPORT(logit_phi);
  }
  
  REPORT(quad);
  REPORT(logDet);
  REPORT(logDetTau);
  REPORT(bym2LogLik);
  REPORT(logPriTau);
  REPORT(KLD);
  REPORT(lexpDensity);
  REPORT(ljacobian);
  REPORT(d);
  REPORT(logPriPhi);
  REPORT(liksUrbDHS);
  REPORT(transformedEpsilon);
  
  return jnll;
  
}
