//modified competing risk model

functions {
  real sev_logpdf(real y, real mu, real sigma){
    real z;
    z = (y - mu) / sigma;
    return -log(sigma) + z - exp(z);
  }
  
  real sev_logccdf(real y, real mu, real sigma){
    return -exp((y - mu) / sigma);
  }
  
  real sev_logcdf(real y, real mu, real sigma){
    return log1m_exp(-exp((y - mu) / sigma));
  }
}

data {
  int M;
  int N_obs;
  int N_cens;
  real endtime_obs[N_obs];
  real endtime_cens[N_cens];
  real starttime_obs[N_obs];
  real starttime_cens[N_cens];
  int<lower=1, upper=M> dm_obs[N_obs];
  int<lower=1, upper=M> dm_cens[N_cens];
  vector<lower=0, upper=1>[2] p; # quantiles to model
}
transformed data{
  vector[2] z_corr;
  for(i in 1:2)
    z_corr[i] = log(-1.0 * log1m(p[i])); # used to get location(mu) from quantile(t_p)
}
parameters{
  real log_tp1;
  real<lower=0> sigma1;
  real<lower=0, upper=1> sigma2[M];
  real eta_tp2;
  real<lower=0> tau_tp2;
  real eta_s2;
  real<lower=0> tau_s2;
  vector[M] log_tp2_raw;
  real logit_pi;
}

transformed parameters{
  real mu1;
  vector[M] mu2;
  real log_pi;
  mu1 = log_tp1 - (sigma1 * z_corr[1]); //change
  for(m in 1:M){
    mu2[m] = (tau_tp2*log_tp2_raw[m] + eta_tp2) - (sigma2[m] * z_corr[2]);
  }
  log_pi = log_inv_logit(logit_pi);
}

model{
  real tmp[2];
  real logpi;
  real m1;
  real m2;
  real s1;
  real s2;
  int m;
  
  logpi = log_pi;
  m1 = mu1;
  s1 = sigma1;
  for(i in 1:N_obs){
  m = dm_obs[i];
  m2 = mu2[m];
  s2 = sigma2[m];
    // numerator:   log( p * f1 * (1 - F2) + f2 * (1 - p * F1) )
    //            = log( exp(log(p) + log(f1) + log(1 - F2)) + exp(log(f2) + log(1 - exp(log(p) + log(F1)))) )
    tmp[1] = log_sum_exp(logpi + sev_logpdf(endtime_obs[i], m1, s1) +
               sev_logccdf(endtime_obs[i], m2, s2),
               sev_logpdf( endtime_obs[i], m2, s2) + 
               log1m_exp(logpi + sev_logcdf(endtime_obs[i], m1, s1)
               )
             );
    // denominator:  log((1 - p * F1) * (1 - F2))
    //            =  log(1 - p * F1) + log(1 - F2)
    tmp[2] = log1m_exp(logpi + sev_logcdf(starttime_obs[i], m1, s1)) + 
             sev_logccdf(starttime_obs[i], m2, s2);
             
    target += tmp[1] - tmp[2];
  }
  
  for(i in 1:N_cens){
    m = dm_cens[i];
    m2 = mu2[m];
    s2 = sigma2[m];
    // numerator:   log((1 - p * F1) * (1 - F2))
    //            =  log(1 - p * F1) + log(1 - F2)
    tmp[1] = log1m_exp(logpi + sev_logcdf(endtime_cens[i], m1, s1)) + 
             sev_logccdf(endtime_cens[i], m2, s2);
    // denominator:  log((1 - p * F1) * (1 - F2))
    //            =  log(1 - p * F1) + log(1 - F2)
    tmp[2] = log1m_exp(logpi + sev_logcdf(starttime_cens[i], m1, s1)) + 
             sev_logccdf(starttime_cens[i], m2, s2);
             
    target += tmp[1] - tmp[2];
  }
  //priors
  log_tp1     ~ normal(7, 2);
  log_tp2_raw ~ normal(0, 1);
  eta_tp2     ~ normal(9, 2);
  tau_tp2     ~ cauchy(0, 1);
  sigma1      ~ lognormal(0, 1);
  sigma2      ~ lognormal(eta_s2, tau_s2);
  eta_s2      ~ normal(0, 2);
  tau_s2      ~ cauchy(0, 1);
  logit_pi   ~ normal(-3, 2);
}
