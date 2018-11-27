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
  DATA_FACTOR(StockID);
  DATA_VECTOR(B_t);
  DATA_VECTOR(P_t);

  // Parameters
  PARAMETER_VECTOR(ln_B0);
  PARAMETER_VECTOR(ln_r);
  PARAMETER_VECTOR(ln_sigmaP);
  PARAMETER_VECTOR(logit_phi);

  int i;
  Type nll=0;
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
  P_t_exp(0) = r(StockID(0))*B_t(0)*(1-B_t(0)/B0(StockID(0)));
  eps_t(0) = P_t(0) - P_t_exp(0);
  // nll -= dnorm(eps_t(0), Type(0.0), sigmaP(StockID(0))/sqrt(1-pow(phi(StockID(0)), 2.0)), true); // Kiva
  nll -= dnorm(eps_t(0), Type(0.0), pow(sigmaP(StockID(0)), 0.5), true); // Jim
  for(int i=1; i<Nobs; i++){
    // exp for expected productivity
    P_t_exp(i) = r(StockID(i))*B_t(i)*(1-B_t(i)/B0(StockID(i)));
    eps_t(i) = P_t(i) - P_t_exp(i);
    // nll -= dnorm(eps_t(i), phi(StockID(i))*eps_t(i-1), sigmaP(StockID(i)), true); // Kiva
    // nll -= dnorm(eps_t(i), phi(StockID(i))*eps_t(i-1), pow(sigmaP(StockID(i)), 0.5), true); // Jim
    // nll -= dnorm(eps_t(i), phi(StockID(i))*eps_t(i-1)+pow(1-phi(StockID(i)),0.5), pow(sigmaP(StockID(i)), 0.5), true); // Jim 2.0
    // nll -= dnorm( P_t(i), P_t_exp(i), sigmaP(StockID(i)), true);
    nll += SCALE( AR1(phi(StockID(i))), pow(sigmaP(StockID(i)) / (1-pow(phi(StockID(i)),2)),0.5))( eps_t(i) );
  }
  
  // penalize for B0 > 5
  Type pen=0;
  Type eps=1e-3;
  vector<Type> B0_check(Nstocks);
  for(int i=0; i<Nstocks; i++){
    B0_check(i) = 5 - B0(i);
    B0_check(i) = posfun(B0_check(i), eps, pen);
    nll += pen;
  }
  
  ADREPORT( r );
  ADREPORT( B0 );
  ADREPORT( sigmaP );
  ADREPORT( phi );
  
  return nll;
}
