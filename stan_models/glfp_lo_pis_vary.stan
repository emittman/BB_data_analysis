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
  real log_tp2;
  real<lower=0> sigma1;
  real<lower=0> sigma2;
  vector[M] logit_pi_raw;
  real eta_pi;
  real<lower=0> tau_pi;
}

transformed parameters{
  real mu1;
  real mu2;
  vector[M] log_pi;
  mu1 = log_tp1 - (sigma1 * z_corr[1]);
  mu2 = log_tp2 - (sigma2 * z_corr[2]);
  for(m in 1:M)
    log_pi[m] = log_inv_logit(eta_pi + tau_pi * logit_pi_raw[m]);
}

model{
  real tmp[2];
  int m;
  real logpi;
  for(i in 1:N_obs){
    m = dm_obs[i];
    logpi = log_pi[m];
    // numerator:   log( p * f1 * (1 - F2) + f2 * (1 - p * F1) )
    //            = log( exp(log(p) + log(f1) + log(1 - F2)) + exp(log(f2) + log(1 - exp(log(p) + log(F1)))) )
    tmp[1] = log_sum_exp(log_pi + sev_logpdf(endtime_obs[i], mu1, sigma1) +
               sev_logccdf(endtime_obs[i], mu2, sigma2),
               sev_logpdf( endtime_obs[i], mu2, sigma2) + 
               log1m_exp(logpi + sev_logcdf(endtime_obs[i], mu1, sigma1)
               )
             );
    // denominator:  log((1 - p * F1) * (1 - F2))
    //            =  log(1 - p * F1) + log(1 - F2)
    tmp[2] = log1m_exp(logpi + sev_logcdf(starttime_obs[i], mu1, sigma1)) + 
             sev_logccdf(starttime_obs[i], mu2, sigma2);
             
    target += tmp[1] - tmp[2];
  }
  
  for(i in 1:N_cens){
    // numerator:   log((1 - p * F1) * (1 - F2))
    //            =  log(1 - p * F1) + log(1 - F2)
    tmp[1] = log1m_exp(logpi + sev_logcdf(endtime_cens[i], mu1, sigma1)) + 
             sev_logccdf(endtime_cens[i], mu2, sigma2);
    // denominator:  log((1 - p * F1) * (1 - F2))
    //            =  log(1 - p * F1) + log(1 - F2)
    tmp[2] = log1m_exp(logpi + sev_logcdf(starttime_cens[i], mu1, sigma1)) + 
             sev_logccdf(starttime_cens[i], mu2, sigma2);
             
    target += tmp[1] - tmp[2];
  }
  
  log_tp1        ~ normal(8, 2);
  log_tp2        ~ normal(10, 2);
  sigma1         ~ lognormal(0, 2.5);
  sigma2         ~ lognormal(0, 2.5);
  logit_pi_raw   ~ normal(0, 1);
  eta_pi         ~ normal(-3, 2);
  tau_pi         ~ cauchy(0, 1);

}
