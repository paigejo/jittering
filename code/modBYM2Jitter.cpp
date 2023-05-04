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

// PC prior on BYM2 phi parameter (on logit scale)
template<class Type>
Type dPCPriPhi(Type logitPhi, Type lambda, Type tr, std::vector<Type> gammaTildesm1, Type logDet, 
               int give_log=0)
{
  Type phi = (Type) (1 / (1 + exp(logitPhi)));
  
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
  Type ljacobian = - log(d) - log(2) + log(abs(tr - n - sumVal));
  
  ljacobian = ljacobian - logitPhi - 2*log(1 + exp(-logitPhi));
  
  // get final log density
  Type lPriorPhi = lexpDensity + ljacobian;
  
  return lPriorPhi;
}

// PC prior on BYM2 tau (precision) parameter (on log scale)
template<class Type>
Type dPCPriTau(Type logTau, Type lambda)
{
  Type tau = exp(logTau);
  
  Type ldensity = log(lambda) - log(2) -(3/2)*logTau - lambda/sqrt(tau);
  
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
  DATA_INTEGER( flag ); // flag == 0 => no data contribution added to jnll
  
  // Indices
  DATA_INTEGER( num_iUrban );   // Number of urban data points in space
  DATA_INTEGER( num_iRural );   // Number of rural data points in space
  DATA_INTEGER( n_integrationPointsUrban );   // number of integration points for each observation
  DATA_INTEGER( n_integrationPointsRural );
  DATA_INTEGER( num_a );   // Number of areas
  
  // Data (all except for X_ij is a vector of length num_i)
  DATA_VECTOR( y_iUrban );   // obs per binomial experiment at point i (clust)
  DATA_VECTOR( y_iRural );
  DATA_VECTOR( n_iUrban );   // Trials per cluster
  DATA_VECTOR( n_iRural );
  DATA_VECTOR( AprojUrban ); // nObsUrban x nArea matrix with ij-th entry = 1 if cluster i associated with area j and 0 o.w.
  DATA_VECTOR( AprojRural ); // nObsRural x nArea matrix with ij-th entry = 1 if cluster i associated with area j and 0 o.w.
  DATA_MATRIX( X_betaUrban );  // (nObsUrban * nIntegrationPointsUrban) x nPar design matrix. Indexed mod numObsUrban
  DATA_MATRIX( X_betaRural );  // first nObsRural rows correspond to first int pt
  DATA_MATRIX( wUrban ); // nObsUrban x nIntegrationPointsUrban weight matrix
  DATA_MATRIX( wRural ); // nObsRural x nIntegrationPointsRural weight matrix
  
  DATA_MATRIX( Q_bym2 );
  
  // BYM2 and prior precomputed values
  DATA_SCALAR( tr );
  DATA_VECTOR( gammaTildesm1 );
  DATA_SCALAR( lambdaPhi );
  DATA_SCALAR( lambdaTau );
  
  // Options
  DATA_VECTOR( options );
  // options[0] == 1 : use normalization trick
  // options[1] == 1 : adreport transformed params
  
  // Fixed effects
  PARAMETER( alpha ); // Intercept
  PARAMETER_VECTOR( beta ); // fixed effect/covariate effect sizes
  // Log of INLA tau param (generalized precision of BYM2)
  PARAMETER( log_tau );
  // Logit of phi (proportion of variance that is structured)
  PARAMETER( logit_phi );
  
  // BYM2 effects (structured + unstructured)
  // NOTE: may need to modify parameterization later
  PARAMETER_VECTOR( Epsilon_bym2 );
  
  // ~~~~~~~~~------------------------------------------------~~
  // SECOND, we define all other objects that we need internally
  // ~~~~~~~~~------------------------------------------------~~
  
  // objective function -- joint negative log-likelihood/posterior
  Type jnll = 0;
  
  // Transform some of our parameters
  Type tau = exp(log_tau);
  Type phi = 1/(1 + exp(-logit_phi));
  
  // Define objects for derived values
  vector<Type> fe_iUrban(num_iUrban * n_integrationPointsUrban); // main effect: alpha
  vector<Type> fe_iRural(num_iRural * n_integrationPointsRural);
  // Logit estimated prob for each cluster i
  //vector<Type> latent_field_iUrban(num_iUrban);
  //vector<Type> latent_field_iRural(num_iRural);
  // value of gmrf at data points
  vector<Type> projepsilon_iUrban(num_iUrban);
  vector<Type> projepsilon_iRural(num_iRural);
  
  // linear combination of fixed effects
  fe_iUrban = X_betaUrban * Type(beta) + Type(alpha); // initialize
  fe_iRural = X_betaRural * Type(beta) + Type(alpha);
  
  // Project GP approx from mesh points to data points
  // projepsilon_iUrban = AprojUrban * Epsilon_s.matrix();
  // projepsilon_iRural = AprojRural * Epsilon_s.matrix();
  projepsilon_iUrban = AprojUrban * Epsilon_bym2;
  projepsilon_iRural = AprojRural * Epsilon_bym2;
  
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
  jnll += GMRF(Q_bym2)(Epsilon_s, false);
  
  // calculate log determinant, and add it to the posterior
  Type logDet = sum(log(1 + phi*gammaTildesm1));
  jnll += logDet
  
  /////////
  // (2) //
  /////////
  // Prior contributions to joint likelihood
  
  // add in priors for BYM2 parameter
  jnll -= dPCPriTau(log_tau, lambdaTau);
  jnll -= dPCPriPhi(logit_phi, lambdaPhi, tr, gammaTildesm1, logDet);
  
  // prior for intercept
  jnll -= dnorm(alpha, alpha_pri[0], alpha_pri[1], true); // N(mean, sd)
  
  // prior for other covariates
  for(int i = 0; i < beta.size(); i++) {
    jnll -= dnorm(beta[i], beta_pri[0], beta_pri[1], true); // N(mean, sd)
  }
  
  /////////
  // (3) //
  /////////
  // jnll contribution from each datapoint i
  int thisIndex;
  Type thisLatentField;
  Type thisWeight;
  Type thislik;
  for (int intI = 0; intI < n_integrationPointsUrban; intI++) {
    thislik = 0;
    
    for (int obsI = 0; obsI < num_iUrban; obsI++) {
      thisIndex = num_iUrban * intI + obsI;
      
      // latent field estimate at each obs
      thisLatentField = fe_iUrban(thisIndex) + projepsilon_iUrban(obsI);
      
      // and add data contribution to jnll
      if(!isNA(y_iUrban(obsI))){
        // get integration weight
        thisWeight = wUrban(obsI,intI);
        
        // Uses the dbinom_robust function, which takes the logit probability
        if(thisWeight > 0) {
          thislik += thisWeight*dbinom_robust( y_iUrban(obsI), n_iUrban(obsI), thisLatentField, false);
        }
        
      } // !isNA
      
    } // for( intI )
    
    jnll -= log(thislik);
  } // for( obsI )
  
  for (int intI = 0; intI < n_integrationPointsRural; intI++) {
    thislik = 0;
    
    for (int obsI = 0; obsI < num_iRural; obsI++) {
      thisIndex = num_iRural * intI + obsI;
      
      // latent field estimate at each obs
      thisLatentField = fe_iRural(thisIndex) + projepsilon_iRural(thisIndex);
      
      // and add data contribution to jnll
      if(!isNA(y_iRural(obsI))){
        // get integration weight
        thisWeight = wRural(obsI,intI);
        
        // Uses the dbinom_robust function, which takes the logit probability
        if(thisWeight > 0) {
          thislik += thisWeight*dbinom_robust( y_iRural(obsI), n_iRural(obsI), thisLatentField, false);
        }
        
      } // !isNA
      
    } // for( intI )
    
    jnll -= log(thislik);
  } // for( obsI )
  
  // ~~~~~~~~~~~
  // ADREPORT: used to return estimates and cov for transforms?
  // ~~~~~~~~~~~
  if(options[1]==1){
    ADREPORT(sp_range);
    ADREPORT(sp_sigma);
  }
  
  return jnll;
  
}
