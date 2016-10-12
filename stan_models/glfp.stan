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
  int M;
  int N_obs;
  int N_cens;
  real endtime_obs[N_obs];
  real endtime_cens[N_cens];
  real starttime_obs[N_obs];
  real starttime_cens[N_cens];
  int<lower=1, upper=M> dm_obs[N_obs];
  int<lower=1, upper=M> dm_cens[N_cens];
  vector<lower=0, upper=1>[2] p;
}
transformed data{
  vector[2] z_corr;
  for(i in 1:2)
    z_corr[i] = log(-1.0 * log1m(p[i]));
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
  real<lower=0, upper=1> mu_pi;
  real log_phi_pi;
  real<lower=0, upper=1> pi[M];
}

transformed parameters{
  vector[M] mu1;
  vector[M] mu2;
  vector[M] log_sigma1;
  vector[M] log_sigma2;
  real<lower=0> alpha;
  real<lower=0> beta; 
  log_sigma1 = eta_ls1 + tau_ls1 * log_sigma1_raw;
  log_sigma2 = eta_ls2 + tau_ls2 * log_sigma2_raw;
  mu1 = (eta_ltp1 + tau_ltp1 * log_tp1_raw) - (exp(log_sigma1) * z_corr[1]);
  mu2 = (eta_ltp2 + tau_ltp2 * log_tp2_raw) - (exp(log_sigma2) * z_corr[2]);
  alpha = mu_pi * exp(log_phi_pi);
  beta = (1-mu_pi) * exp(log_phi_pi);
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
    m = dm_obs[i];
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
    m = dm_cens[i];
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
  
  log_phi_pi ~ normal(4, 1);
  
  log_tp1_raw ~ student_t(5, 0,1);
  log_tp2_raw ~ student_t(5, 0,1);
  log_sigma1_raw ~ normal(0,1);
  log_sigma2_raw ~ normal(0,1);
  
  pi ~ beta(alpha, beta);
  eta_ltp1 ~ normal(7, 1);
  eta_ls1  ~ normal(0, 1);
  eta_ltp2 ~ normal(9, 2);
  eta_ls2  ~ normal(0, 1);
  tau_ltp1 ~ normal(0, 1);
  tau_ltp2 ~ normal(0, 1);
  tau_ls1  ~ normal(0, .5);
  tau_ls2  ~ normal(0, .5);

}
