library(remakeGenerator)

datasets = commands(
#  data_2_21 = prepare_data(lb_fails = 5, lb_late_fails = 1, lb_early_fails = 1),
  data_3_10 = prepare_data(lb_fails = 3, lb_late_fails = 0, lb_early_fails = 0)
)

analyses = analyses(
  commands = commands(
 #   samples_logodds_reduced_3_8 = run_mcmc_logodds(..dataset.., chains = 16, iter = 3000, warmup = 1500),
    samples_double_variance_3_25 = run_mcmc_generic(dataset = data_3_10, chains = 16,
                                                    iter = 3000, warmup = 1500,
                                                    stanmodel = "../stan_models/glfp_lo_reduced_doublyrelaxed.stan")
  ),
  datasets = datasets
)

targets = targets(datasets = datasets, analyses = analyses)

workflow(targets, sources = "functions.R", packages = "rstan")