# Functions for use with remakeGenerator

prepare_data <- function(dm = NULL, lb_fails = 0, lb_late_fails = 0, lb_early_fails = 0){
  require(dplyr)
  #load data
  overview  <- readRDS("../BB_data/overview.rds")
  dat       <- readRDS("../BB_data/clean_unit_summaries.rds")
  #subsetting
  id        <- with(overview, which(f >= lb_fails & late_f >= lb_late_fails & early_f >= lb_early_fails))
  df <- with(subset(dat, model %in% id),
               data.frame(endtime = endtime,
                          starttime = starttime,
                          censored = censored,
                          model = as.integer(factor(model)))
  )
  if(!is.null(dm)){
    if(sum(df$model==dm) == 0) stop("There are no drive models matching the argument")
    df <- filter(df, model == dm)
  }
  # return(overview)
  return(df)
}

# run_mcmc <- function(dataset, chains = 4, iter = 2000, warmup = iter/2, p = c(.5,.2)){
#   require(rstan)
#   rstan_options(auto_write = TRUE)
#   options(mc.cores = parallel::detectCores())#format for Stan
#   stan_data <- with(dataset,
#                     list(M = length(unique(model)),
#                          N_obs = sum(!censored),
#                          N_cens = sum(censored),
#                          starttime_obs = log(starttime[!censored]+1),
#                          starttime_cens = log(starttime[censored]+1),
#                          endtime_obs = log(endtime[!censored]+1),
#                          endtime_cens = log(endtime[censored]+1),
#                          dm_obs = model[!censored],
#                          dm_cens = model[censored],
#                          p = p)
#   )
#   
#   m <- stan_model(file = "../stan_models/glfp_logodds_sigma_centered.stan")
#   s <- sampling(obj = m, data = stan_data, chains = chains, iter = iter, warmup = warmup, control = list(adapt_delta = .99), init=0)
#   return(s)
#   # return(stan_data)
# }

run_mcmc_logodds <- function(dataset, chains = 4, iter = 10, warmup = iter/2, p = c(.5,.2)){
  require(rstan)
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())#format for Stan
  stan_data <- with(dataset,
                    list(M = length(unique(model)),
                         N_obs = sum(!censored),
                         N_cens = sum(censored),
                         starttime_obs = log(starttime[!censored]+1),
                         starttime_cens = log(starttime[censored]+1),
                         endtime_obs = log(endtime[!censored]+1),
                         endtime_cens = log(endtime[censored]+1),
                         dm_obs = model[!censored],
                         dm_cens = model[censored],
                         p = p)
  )
  

  m <- stan_model(file = "../stan_models/glfp_lo_reduced_relaxed.stan")
  s <- sampling(obj = m, data = stan_data, chains = chains, iter = iter, warmup = warmup, control = list(adapt_delta = .999), init=0)
  return(s)
  # return(stan_data)
}

run_mcmc_logodds_nopool <- function(dataset, chains = 4, iter = 10, warmup = iter/2, p = c(.5,.2)){
  require(rstan)
  rstan_options(auto_write = TRUE)
  options(mc.cores = 4)#format for Stan
  stan_data <- with(dataset,
                    list(N_obs = sum(!censored),
                         N_cens = sum(censored),
                         starttime_obs = log(starttime[!censored]+1),
                         starttime_cens = log(starttime[censored]+1),
                         endtime_obs = log(endtime[!censored]+1),
                         endtime_cens = log(endtime[censored]+1),
                         p = p)
  )
  
  
  m <- stan_model(file = "../stan_models/glfp_no_pooling.stan")
  s <- sampling(obj = m, data = stan_data, chains = chains, iter = iter, warmup = warmup, control = list(adapt_delta = .999), init=0)
  return(s)
  # return(stan_data)
}

count_divergences <- function(fit, chain, inc_warmup) {
  sampler_params <- get_sampler_params(fit, inc_warmup=inc_warmup)
  sum(sapply(sampler_params, function(x) c(x[,'divergent__']))[,chain])
}

hist_treedepth <- function(fit, chain, inc_warmup) {
  sampler_params <- get_sampler_params(fit, inc_warmup=inc_warmup)
  hist(sapply(sampler_params, function(x) c(x[,'treedepth__']))[,chain], breaks=0:20, main="", xlab="Treedepth")
  abline(v=10, col=2, lty=1)
}
