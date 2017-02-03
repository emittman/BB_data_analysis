library(remakeGenerator)

datasets = commands(
  data_2_3 = prepare_data(lb_fails = 5, lb_late_fails = 0, lb_early_fails = 0)
)

analyses = analyses(
  commands = commands(
    samples_logodds_2_3 = run_mcmc_logodds(..dataset.., chains = 16, iter = 3000, warmup = 1500)
  ),
  datasets = datasets
)

targets = targets(datasets = datasets, analyses = analyses)

workflow(targets, sources = "functions.R", packages = "rstan")