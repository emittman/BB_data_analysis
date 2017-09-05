library(drake)

prepare_data <- function(min.fails = 0){
  require(dplyr)
  #load data
  overview  <- readRDS("../BB_data/overview.rds")
  dat       <- readRDS("../BB_data/clean_unit_summaries.rds")
  #subsetting
  id        <- with(overview, which(f >= min.fails))
  df <- with(subset(dat, model %in% id),
             data.frame(endtime = endtime,
                        starttime = starttime,
                        censored = censored,
                        model = as.integer(factor(model)))
  )
  # if(!is.null(dm)){
  #   if(sum(df$model==dm) == 0) stop("There are no drive models matching the argument")
  #   df <- filter(df, model == dm)
  # }
  # return(overview)
  return(df)
}

run_mcmc_generic <- function(dataset, chains = 4, iter = 10, warmup = iter/2, p = c(.5,.2), stanmodel=NULL){
  if(is.null(stanmodel)) stop("No stanmodel was supplied.")
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
  
  
  m <- stan_model(file = stanmodel)
  s <- sampling(obj = m, data = stan_data, chains = chains, iter = iter, warmup = warmup, control = list(adapt_delta = .99), init=0)
  return(s)
  # return(stan_data)
}

datasets <- plan(data = prepare_data(min.fails = 0))


methods <- plan(null_model = run_mcmc_generic(..datasets.., stanmodel="../stan_models/glfp_lo_null_model.stan"),
                strings_in_dots = "filenames")

analyses <- analyses(methods,datasets)

make(analyses, jobs=4)
