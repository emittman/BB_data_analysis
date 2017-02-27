library(remakeGenerator)

datasets = commands(
  data_2_21 = prepare_data(lb_fails = 5, lb_late_fails = 1, lb_early_fails = 1)
)

analyses = analyses(
  commands = commands(
    samples_logodds_reduced_2_28 = run_mcmc_logodds(..dataset.., chains = 16, iter = 3000, warmup = 1500)
  ),
  datasets = datasets
)

targets = targets(datasets = datasets, analyses = analyses)

workflow(targets, sources = "functions.R", packages = "rstan")