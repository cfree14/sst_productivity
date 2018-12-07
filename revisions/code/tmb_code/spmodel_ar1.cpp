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
  // PARAMETER_VECTOR(logit_rho); // comment out when fixing rho=0
  PARAMETER_VECTOR(rho); // uncomment when fixing rho=0
  // PARAMETER(rho); // global rho
  
  Type nll=0;
  vector<Type> B0(Nstocks);
  vector<Type> r(Nstocks);
  vector<Type> sigmaP(Nstocks);
  // vector<Type> rho(Nstocks); // comment out when fixing rho=0
  for(int i=0; i<Nstocks; i++){
    B0(i) = exp(ln_B0(i));
    r(i) = exp(ln_r(i));
    sigmaP(i) = exp(ln_sigmaP(i));
    // rho(i) = exp(logit_rho(i)) / (1 + exp(logit_rho(i))); // inverse logit, comment out when fixing rho=0
  }
  
  // Likelihood contribution from observations (additive process error in production)
  vector<Type> eps_t(Nobs); // residuals
  vector<Type> P_t_exp(Nobs); // predicted productivity
  // Calculate first prediction and residual
  P_t_exp(0) = r(StockID(0))*B_t(0)*(1-B_t(0)/B0(StockID(0)));
  eps_t(0) = P_t(0) - P_t_exp(0);
  nll -= dnorm(eps_t(0), Type(0.0), sigmaP(StockID(0)), true); 
  // nll -= dnorm(eps_t(0), Type(0.0), sqrt(pow(sigmaP(StockID(0)), 2.0)/1-pow(rho(StockID(0)), 2.0)), true); // dnorm(value, mean, sd), sd=sqrt(var); var=sd^2, var of AR1 = sigmaP^2/(1-rho^2), sd = sqrt(sigmaP^2/(1-rho^2))
  for(int i=1; i<Nobs; i++){
    // exp for expected productivity
    P_t_exp(i) = r(StockID(i))*B_t(i)*(1-B_t(i)/B0(StockID(i)));
    eps_t(i) = P_t(i) - P_t_exp(i);
    if(StockID(i) != StockID(i-1)){
      nll -= dnorm(eps_t(i), Type(0.0), sigmaP(StockID(i)), true); 
      // nll -= dnorm(eps_t(i), Type(0.0), sqrt(pow(sigmaP(StockID(i)), 2.0)/1-pow(rho(StockID(i)), 2.0)), true); // from above
    }else{
      nll -= dnorm(eps_t(i), rho(StockID(i))*eps_t(i-1), sigmaP(StockID(i)), true); // rho for each stock
      // nll -= dnorm(eps_t(i), rho*eps_t(i-1), sigmaP(StockID(i)), true); // global rho
    }
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
  
  // penalize for rhos deviating from 0
  // Type rho_pen=0; 
  // vector<Type> rho2(Nstocks);
  // vector<Type> rho2_check(Nstocks);
  // for(int i=0; i<Nstocks; i++){
  //   rho2(i) = pow(rho(i), 2.0);
  //   rho2_check(i) = pow(0.8, 2.0) - rho2(i);
  //   rho2_check(i) = posfun(rho2_check(i), eps, rho_pen);
  //   nll += rho_pen;
  // }
  
  ADREPORT( r );
  ADREPORT( B0 );
  ADREPORT( sigmaP );
  // ADREPORT( rho );

  return nll;
}
