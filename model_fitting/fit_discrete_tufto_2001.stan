// This stan program fits a discrete time genrealization of the model 
// described in tufto 2001. some simplificaitons and aproximation were made 
// inorder to produce a model that could be fit to "data" the main idea 
// behind this model is not neccisarly to make inferences about real
// existing data sets, but to understand if and how information about 
// migrations and abundance can be used to make inferences about 
// the underlying populaiton genetic process. 
//
// In addition to tracking the abundnace of the population we track
// the survival of 
data {
  
  // meta data 
  int<lower=0> N;
  
  // data
  vector[N] y; // log(n_{t+1}-m_t) - log(n_t)
  vector[N] N_t; 
  vector[N] m_t;
  int n[N,3];
  
  // priors
  real mu_a;
  real sigma_a;
  real alpha_b;
  real beta_b;
  real alpha_s;
  real beta_s;
  real lambda_z;
  real lambda_sigma;
  real a_rho;
  real b_rho;
  real sigma_z;
  real alpha_sigma_nu;
  real beta_sigma_nu;
}

transformed data{
  vector[N] mu;
  mu[1] = 0.0;
  for(i in 2:N){
    mu[i] = m_t[i-1]/N_t[i-1];
  }
}

parameters {
  real a;
  real<lower=0> b;
  real<lower=0> s;
  real<lower=0> z1;
  real<lower=0> sigma;
  real<lower=0,upper=1> rho;
  real z0; // inital value fo z
  vector[N-1] nu;
  real<lower=0>sigma_nu;

}

transformed parameters{
  vector[N] z;
  vector[N] a_t;
  vector[N] resids;
  vector[3] p[N];
  real K;
  z[1] = z0;//epsilon_z[1];
  a_t[1] = a;
  for(i in 2:N){
    z[i] = (N_t[i] -m_t[i-1])*z[i-1]/(N_t[i]*(1+s)) + m_t[i-1]/N_t[i]*z1 + nu[i-1];
    a_t[i] = a - 0.5*s*z[i]^2/(1+s);
  }
  
  for(i in 1:N){
    
    p[i,1] = (1-mu[i])^2*exp(-0.5*s*z[i]^2/(1+s));
    p[i,2] = 2*mu[i]*(1-mu[i])*exp(-0.5*s*((z[i] + z1)/2)^2/(1+s));
    p[i,3] = mu[i]^2*exp(-0.5*s*z1^2/(1+s));

    p[i] = p[i]/sum(p[i]);
  }

  for(i in 1:N){
    resids[i] = y[i] - a_t[i] + b*N_t[i];
  }
  K = b/a;

}

model {
  // priors
  a ~ normal(mu_a,sigma_a);
  target += gamma_lpdf(K|alpha_b,beta_b);
  s ~ gamma(alpha_s,beta_s);
  z1 ~ exponential(lambda_z);
  sigma ~ exponential(lambda_sigma);
  rho~ beta(a_rho,b_rho);
  z0 ~ normal(0, sigma_z);
  nu ~ normal(0, sigma_nu);
  sigma_nu ~ gamma(alpha_sigma_nu,beta_sigma_nu);
  //epsilon_z ~ normal(0,0.05);
  
  //likelihood
  target += normal_lpdf(resids[1]|0,sigma);
  for(i in 2:N){
    target += normal_lpdf(resids[i]|rho*resids[i-1],sigma);
  }

  for(i in 1:N){
    n[i] ~ multinomial(p[i]);
  }
  
}

generated quantities{
  real r;
  r = a;
}

