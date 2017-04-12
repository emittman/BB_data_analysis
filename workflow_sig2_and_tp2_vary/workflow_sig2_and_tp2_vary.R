library(remakeGenerator)

datasets = commands(
  data_all = prepare_data(lb_fails = 3, lb_late_fails = 0, lb_early_fails = 0)
)

analyses = analyses(
  commands = commands(
    samples_lo_sig2_and_tp2_vary = run_mcmc_generic(..dataset.., chains = 4,
                                                    iter = 3000, warmup = 1500,
                                                    stanmodel=I("../stan_models/glfp_sig2_tp2_vary.stan"))
  ),
  datasets = datasets
)

targets = targets(datasets = datasets, analyses = analyses)

workflow(targets, sources = "../workflow/functions.R", packages = "rstan")