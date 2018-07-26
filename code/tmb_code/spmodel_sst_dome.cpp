
// Source TMB package
#include <TMB.hpp>


// Setup objective function
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_INTEGER(Nstocks);
  DATA_INTEGER(Nobs);
  DATA_FACTOR(StockID);
  DATA_VECTOR(B_t);
  DATA_VECTOR(P_t);
  DATA_VECTOR(Temp_t);

  // Parameters
  PARAMETER_VECTOR(ln_B0);
  PARAMETER_VECTOR(ln_r);
  PARAMETER_VECTOR(ln_BetaT);
  PARAMETER_VECTOR(Z);
  PARAMETER_VECTOR(ln_sigmaP);
  PARAMETER(mu_T);
  PARAMETER(ln_sd_T);

  // Transform back to real space
  Type nll=0;
  Type sd_T = exp(ln_sd_T);
  vector<Type> B0(Nstocks);
  vector<Type> r(Nstocks);
  vector<Type> BetaT(Nstocks);
  vector<Type> sigmaP(Nstocks);
  for(int i=0; i<Nstocks; i++){
    B0(i) = exp(ln_B0(i));
    r(i) = exp(ln_r(i));
    BetaT(i) = exp(ln_BetaT(i));
    sigmaP(i) = exp(ln_sigmaP(i));
  }
  
  // Likelihood contribution from observations (additive process error in production)
  vector<Type> P_t_pred(Nobs);
  for(int i=0; i<Nobs; i++){
    // predicted productivity
    P_t_pred(i) = r(StockID(i))*B_t(i)*(1-B_t(i)/B0(StockID(i))) *
      exp(-pow(Temp_t(i) - Z(StockID(i)), 2) * BetaT(StockID(i)));
    nll -= dnorm( P_t(i), P_t_pred(i), sigmaP(StockID(i)), true);
  }

  // Probability of random effects
  for(int i=0; i<Nstocks; i++){
    nll -= dnorm( BetaT(i), mu_T, sd_T, true);
  }
  
  // Calculate the standard erros for the "derived" parameters
  // These are derived because we estimated them in log space
  ADREPORT( r );
  ADREPORT( B0 );
  ADREPORT( BetaT );
  ADREPORT( sigmaP );
  ADREPORT( sd_T );
  
  return nll;
}
