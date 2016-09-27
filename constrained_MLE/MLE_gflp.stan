//modified competing risk model

functions {
  real sev_logpdf(real y, real mu, real log_sigma){
    real z;
    z = (y - mu) / exp(log_sigma);
    return -log_sigma + z - exp(z);
  }
  
  real sev_logccdf(real y, real mu, real log_sigma){
    return -exp((y - mu) / exp(log_sigma));
  }
  
  real sev_logcdf(real y, real mu, real log_sigma){
    return log1m_exp(-exp((y - mu) / exp(log_sigma)));
  }
}

data {
  int N_obs;
  int N_cens;
  int N_obs2;
  int N_cens2;
  int M;
  real endtime_obs[N_obs];
  real endtime_cens[N_cens];
  real starttime_obs[N_obs];
  real starttime_cens[N_cens];
  real endtime_obs2[N_obs2];
  real endtime_cens2[N_cens2];
  real starttime_obs2[N_obs2];
  real starttime_cens2[N_cens2];
  int<lower=1, upper=M> model_obs[N_obs];
  int<lower=1, upper=M> model_cens[N_cens];
  int<lower=1, upper=M> model_obs2[N_obs2];
  int<lower=1, upper=M> model_cens2[N_cens2];
  vector<lower=0, upper=1>[2] p;
}
transformed data{
  vector[2] z_corr;
  real<lower=0> a;
  real<lower=0> b;
  for(i in 1:2)
    z_corr[i] = log(-1.0 * log1m(p[i]));
  a = 1;
  b = 9;
}
parameters{
  real eta_ltp1;
  real eta_ltp2;
  real eta_ls1;
  real eta_ls2;
  real<lower=0> tau_ltp1;
  real<lower=0> tau_ltp2;
  real<lower=0> tau_ls1;
  real<lower=0> tau_ls2;
  vector[M] log_tp1_raw;
  vector[M] log_tp2_raw;
  vector[M] log_sigma1_raw;
  vector[M] log_sigma2_raw;
  real<lower=0, upper=1> pi[M];
}

transformed parameters{
  vector[M] mu1;
  vector[M] mu2;
  vector[M] log_sigma1;
  vector[M] log_sigma2;
  log_sigma1 = eta_ls1 + tau_ls1 * log_sigma1_raw;
  log_sigma2 = eta_ls2 + tau_ls2 * log_sigma2_raw;
  mu1 = (eta_ltp1 + tau_ltp1 * log_tp1_raw) - (exp(log_sigma1) * z_corr[1]);
  mu2 = (eta_ltp2 + tau_ltp2 * log_tp2_raw) - (exp(log_sigma2) * z_corr[2]);
}

model{
  real tmp[2];
  int m;
  real logpi;
  real mu_1;
  real mu_2;
  real ls_1;
  real ls_2;
  
  for(i in 1:N_obs){
    m = model_obs[i];
    logpi = log(pi[m]);
    mu_1 = mu1[m];
    mu_2 = mu2[m];
    ls_1 = log_sigma1[m];
    ls_2 = log_sigma2[m];
    // numerator:   log( p * f1 * (1 - F2) + f2 * (1 - p * F1) )
    //            = log( exp(log(p) + log(f1) + log(1 - F2)) + exp(log(f2) + log(1 - exp(log(p) + log(F1)))) )
    tmp[1] = log_sum_exp(logpi + sev_logpdf(endtime_obs[i], mu_1, ls_1) +
               sev_logccdf(endtime_obs[i], mu_2, ls_2),
               sev_logpdf( endtime_obs[i], mu_2, ls_2) + 
               log1m_exp(logpi + sev_logcdf(endtime_obs[i], mu_1, ls_1)
               )
             );
    // denominator:  log((1 - p * F1) * (1 - F2))
    //            =  log(1 - p * F1) + log(1 - F2)
    tmp[2] = log1m_exp(logpi + sev_logcdf(starttime_obs[i], mu_1, ls_1)) + 
             sev_logccdf(starttime_obs[i], mu_2, ls_2);
             
    target += tmp[1] - tmp[2];
  }
  
  for(i in 1:N_cens){
    m = model_cens[i];
    logpi = log(pi[m]);
    mu_1 = mu1[m];
    mu_2 = mu2[m];
    ls_1 = log_sigma1[m];
    ls_2 = log_sigma2[m];
  
    // numerator:   log((1 - p * F1) * (1 - F2))
    //            =  log(1 - p * F1) + log(1 - F2)
    tmp[1] = log1m_exp(logpi + sev_logcdf(endtime_cens[i], mu_1, ls_1)) + 
             sev_logccdf(endtime_cens[i], mu_2, ls_2);
    // denominator:  log((1 - p * F1) * (1 - F2))
    //            =  log(1 - p * F1) + log(1 - F2)
    tmp[2] = log1m_exp(logpi + sev_logcdf(starttime_cens[i], mu_1, ls_1)) + 
             sev_logccdf(starttime_cens[i], mu_2, ls_2);
             
    target += tmp[1] - tmp[2];
  }
  
  //Likelihood terms for known Mode 2 (wearout) failures
  
    for(i in 1:N_obs2){
    m = model_obs2[i];
    mu_2 = mu2[m];
    ls_2 = log_sigma2[m];
    
    // numerator:   = log(f2) 
    tmp[1] = sev_logpdf(endtime_obs2[i], mu_2, ls_2);
    
    // denominator:  log(1 - F2)
    tmp[2] = sev_logccdf(starttime_obs2[i], mu_2, ls_2);
             
    target += tmp[1] - tmp[2];
  }
  
  // Likelihood terms for known Mode 2 (wearout) censored observations
  
    for(i in 1:N_cens2){
    m = model_cens2[i];
    mu_2 = mu2[m];
    ls_2 = log_sigma2[m];
  
    // numerator:   log(1 - F2)
    tmp[1] = sev_logccdf(endtime_cens2[i], mu_2, ls_2);
    // denominator:  log(1 - F2)
    tmp[2] = sev_logccdf(starttime_cens2[i], mu_2, ls_2);
             
    target += tmp[1] - tmp[2];
  }

}
