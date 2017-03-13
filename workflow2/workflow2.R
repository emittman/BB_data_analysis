library(remakeGenerator)

datasets = commands(
  data_2_21 = prepare_data(lb_fails = 5, lb_late_fails = 1, lb_early_fails = 1)
)

analyses = analyses(
  commands = commands(
    samples_logodds_nopool_3_13 = run_mcmc_logodds_nopool(..dataset.., chains = 4, iter = 3000, warmup = 1500)
  ),
  datasets = datasets
)

targets = targets(datasets = datasets, analyses = analyses)

workflow(targets, sources = "../workflow/functions.R", packages = "rstan")