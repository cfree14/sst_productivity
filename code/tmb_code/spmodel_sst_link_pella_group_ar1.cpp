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
  DATA_INTEGER(Ngroups);
  DATA_SCALAR(p);
  DATA_FACTOR(StockID);
  DATA_FACTOR(Group);
  DATA_VECTOR(B_t);
  DATA_VECTOR(P_t);
  DATA_VECTOR(Temp_t);

  // Parameters
  PARAMETER_VECTOR(ln_B0);
  PARAMETER_VECTOR(ln_r);
  PARAMETER_VECTOR(BetaT);
  PARAMETER_VECTOR(ln_sigmaP);
  PARAMETER_VECTOR(logit_phi);
  PARAMETER(mu_T);
  PARAMETER(ln_sd_T);
  PARAMETER_VECTOR(mu_group);
  PARAMETER(ln_sd_group);

  int i;
  Type nll=0;
  Type sd_T = exp(ln_sd_T);
  Type sd_group = exp(ln_sd_group);
  vector<Type> B0(Nstocks);
  vector<Type> r(Nstocks);
  vector<Type> sigmaP(Nstocks);
  vector<Type> phi(Nstocks);
  for(int i=0; i<Nstocks; i++){
    B0(i) = exp(ln_B0(i));
    r(i) = exp(ln_r(i));
    sigmaP(i) = exp(ln_sigmaP(i));
    // phi(i) = exp(logit_phi(i)) / (1 + exp(logit_phi(i))); // Kiva
    phi(i) = 1 / (1 + exp(-logit_phi(i))); // Jim
  }
  
  // Likelihood contribution from observations (additive process error in production)
  vector<Type> eps_t(Nobs);
  vector<Type> P_t_exp(Nobs);
  P_t_exp(0) = r(StockID(0))/p*B_t(0)*(1-pow(B_t(0)/B0(StockID(0)),p)) * exp(Temp_t(0)*BetaT(StockID(0)));
  eps_t(0) = P_t(0) - P_t_exp(0);
  // nll -= dnorm(eps_t(0), Type(0.0), sigmaP(StockID(0))/sqrt(1-pow(phi(StockID(0)), 2.0)), true); // Kiva
  nll -= dnorm(eps_t(0), Type(0.0), phi(StockID(0)), true); // Jim
  for(int i=1; i<Nobs; i++){
    // exp for expected productivity
    P_t_exp(i) = r(StockID(i))/p*B_t(i)*(1-pow(B_t(i)/B0(StockID(i)),p)) * exp(Temp_t(i)*BetaT(StockID(i)));
    eps_t(i) = P_t(i) - P_t_exp(i);
    // nll -= dnorm(eps_t(i), phi(StockID(i))*eps_t(i-1), sigmaP(StockID(i)), true); // Kiva
    nll -= dnorm(eps_t(i), phi(StockID(i))*eps_t(i-1), sigmaP(StockID(i)), true); // Jim
  }

  // Probability of random effects
  for(int i=0; i<Ngroups; i++){
    nll -= dnorm( mu_group(i), mu_T, sd_T, true);
  }
  
  // Probability of random effects
  Type pen=0;
  Type eps=1e-3;
  vector<Type> B0_check(Nstocks);
  for(int i=0; i<Nstocks; i++){
    nll -= dnorm( BetaT(i), mu_group(Group(i)), sd_group, true);
    // penalize for B0 > 5
    B0_check(i) = 5 - B0(i);
    B0_check(i) = posfun(B0_check(i), eps, pen);
    nll += pen;
  }
  
  // Report productivity predictions/observations and B observations
  // matrix<Type> resids(Nobs, 3);
  // for(int i=0; i<Nobs; i++){
  //   resids(i,1) = B_t(i);
  //   resids(i,2) = P_t(i);
  //   resids(i,3) = P_t_exp(i);
  // }
  // REPORT( resids );
  
  ADREPORT( r );
  ADREPORT( B0 );
  ADREPORT( sigmaP );
  ADREPORT( sd_T );
  ADREPORT( sd_group );
  
  return nll;
}
