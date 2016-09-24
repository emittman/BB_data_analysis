library(remakeGenerator)

datasets = commands(
  data_8_1 = prepare_data(lb_fails = 8, lb_late_fails = 1, lb_early_fails = 0)
)

analyses = analyses(
  commands = commands(
    samples = run_mcmc(..dataset.., chains = 4, iter = 8000)
  ),
  datasets = datasets
)

targets = targets(datasets = datasets, analyses = analyses)

workflow(targets, sources = "functions.R", packages = "rstan")