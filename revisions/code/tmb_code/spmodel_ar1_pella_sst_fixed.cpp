// Space time 
#include <TMB.hpp>

// Posfun - penalize when below 0
template<class Type>
Type posfun(Type x, Type eps, Type &pen){
  pen += CppAD::CondExpLt(x,eps,Type(0.01)*pow(x-eps,2),Type(0));
  return CppAD::CondExpGe(x,eps,x,eps/(Type(2)-x/eps));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_INTEGER(Nstocks);
  DATA_INTEGER(Nobs);
  DATA_SCALAR(p);
  DATA_FACTOR(StockID);
  DATA_VECTOR(B_t);
  DATA_VECTOR(P_t);
  DATA_VECTOR(Temp_t);

  // Parameters
  PARAMETER_VECTOR(ln_B0);
  PARAMETER_VECTOR(ln_r);
  PARAMETER_VECTOR(BetaT);
  PARAMETER_VECTOR(ln_sigmaP);
  PARAMETER_VECTOR(rho); 

  // Transform parameters
  Type nll=0;
  vector<Type> B0(Nstocks);
  vector<Type> r(Nstocks);
  vector<Type> sigmaP(Nstocks);
  for(int i=0; i<Nstocks; i++){
    B0(i) = exp(ln_B0(i));
    r(i) = exp(ln_r(i));
    sigmaP(i) = exp(ln_sigmaP(i));
  }
  
  // Likelihood contribution from observations (additive process error in production)
  vector<Type> eps_t(Nobs); // residuals
  vector<Type> P_t_exp(Nobs); // predicted productivity
  // Calculate first prediction and residual
  P_t_exp(0) = r(StockID(0))/p*B_t(0)*(1-pow(B_t(0)/B0(StockID(0)),p)) * exp(Temp_t(0)*BetaT(StockID(0)));
  eps_t(0) = P_t(0) - P_t_exp(0);
  nll -= dnorm(eps_t(0), Type(0.0), sigmaP(StockID(0)), true); 
  for(int i=1; i<Nobs; i++){
    // exp for expected productivity
    P_t_exp(i) = r(StockID(i))/p*B_t(i)*(1-pow(B_t(i)/B0(StockID(i)),p)) * exp(Temp_t(i)*BetaT(StockID(i)));
    eps_t(i) = P_t(i) - P_t_exp(i);
    if(StockID(i) != StockID(i-1)){
      nll -= dnorm(eps_t(i), Type(0.0), sigmaP(StockID(i)), true); 
    }else{
      nll -= dnorm(eps_t(i), rho(StockID(i))*eps_t(i-1), sigmaP(StockID(i)), true); // rho for each stock
    }
  }

  // Probability of random effects
  Type pen=0;
  Type eps=1e-3;
  vector<Type> B0_check(Nstocks);
  for(int i=0; i<Nstocks; i++){
    // penalize for B0 > 5
    B0_check(i) = 5 - B0(i);
    B0_check(i) = posfun(B0_check(i), eps, pen);
    nll += pen;
  }

  // Report transformed parameters  
  ADREPORT( r );
  ADREPORT( B0 );
  ADREPORT( sigmaP );

  return nll;
}
