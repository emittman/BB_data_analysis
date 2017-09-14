library(remakeGenerator)

datasets = commands(
#  data_2_21 = prepare_data(lb_fails = 5, lb_late_fails = 1, lb_early_fails = 1),
  data_5_4 = prepare_data(lb_fails = 3, lb_late_fails = 0, lb_early_fails = 0)
)

analyses = analyses(
  commands = commands(
 #   samples_logodds_reduced_3_8 = run_mcmc_logodds(..dataset.., chains = 16, iter = 3000, warmup = 1500),
     s_tpvary = run_mcmc_generic(data_5_4, chains = 4, p = I(c(.8,.2)), iter = 1000, warmup = 800, I("glfp_tp2_vary.stan")),
     s_tpsigvary = run_mcmc_generic(data_5_4, chains = 4, p = I(c(.8,.2)), iter = 1000, warmup = 800, I("glfp_sig2_tp2_vary.stan"))
  ),
  datasets = datasets
)

targets = targets(datasets = datasets, analyses = analyses)

workflow(targets, sources = "../workflow/functions.R", packages = "rstan")