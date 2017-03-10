library(remakeGenerator)

datasets = commands(
  data_2_21 = prepare_data(lb_fails = 5, lb_late_fails = 1, lb_early_fails = 1),
  data_3_10 = prepare_data(lb_fails = 3, lb_late_fails = 0, lb_early_fails = 0)
)

analyses = analyses(
  commands = commands(
    samples_logodds_reduced_2_21 = run_mcmc_logodds(..dataset.., chains = 4, iter = 4000, warmup = 1000)
  ),
  datasets = datasets
)

targets = targets(datasets = datasets, analyses = analyses)

workflow(targets, sources = "functions.R", packages = "rstan")