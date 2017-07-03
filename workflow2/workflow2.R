library(remakeGenerator)

datasets = commands(
  # data_dm1  = prepare_data(dm = 1, lb_fails = 5, lb_late_fails = 1, lb_early_fails = 1),
  # data_dm2  = prepare_data(dm = 2, lb_fails = 5, lb_late_fails = 1, lb_early_fails = 1),
  data_dm3  = prepare_data(dm = 3, lb_fails = 5, lb_late_fails = 1, lb_early_fails = 1),
  # data_dm4  = prepare_data(dm = 4, lb_fails = 5, lb_late_fails = 1, lb_early_fails = 1),
  # data_dm5  = prepare_data(dm = 5, lb_fails = 5, lb_late_fails = 1, lb_early_fails = 1),
  # data_dm6  = prepare_data(dm = 6, lb_fails = 5, lb_late_fails = 1, lb_early_fails = 1),
  # data_dm7  = prepare_data(dm = 7, lb_fails = 5, lb_late_fails = 1, lb_early_fails = 1),
  data_dm8  = prepare_data(dm = 8, lb_fails = 5, lb_late_fails = 1, lb_early_fails = 1),
  # data_dm9  = prepare_data(dm = 9, lb_fails = 5, lb_late_fails = 1, lb_early_fails = 1),
  data_dm10 = prepare_data(dm =10, lb_fails = 5, lb_late_fails = 1, lb_early_fails = 1),
  # data_dm11 = prepare_data(dm =11, lb_fails = 5, lb_late_fails = 1, lb_early_fails = 1),
  # data_dm12 = prepare_data(dm =12, lb_fails = 5, lb_late_fails = 1, lb_early_fails = 1),
  # data_dm13 = prepare_data(dm =13, lb_fails = 5, lb_late_fails = 1, lb_early_fails = 1),
  data_dm14 = prepare_data(dm =14, lb_fails = 5, lb_late_fails = 1, lb_early_fails = 1)#,
  # data_dm15 = prepare_data(dm =15, lb_fails = 5, lb_late_fails = 1, lb_early_fails = 1),
  # data_dm16 = prepare_data(dm =16, lb_fails = 5, lb_late_fails = 1, lb_early_fails = 1),
  # data_dm17 = prepare_data(dm =17, lb_fails = 5, lb_late_fails = 1, lb_early_fails = 1),
  # data_dm18 = prepare_data(dm =18, lb_fails = 5, lb_late_fails = 1, lb_early_fails = 1),
  # data_dm19 = prepare_data(dm =19, lb_fails = 5, lb_late_fails = 1, lb_early_fails = 1),
  # data_dm20 = prepare_data(dm =20, lb_fails = 5, lb_late_fails = 1, lb_early_fails = 1),
  # data_dm21 = prepare_data(dm =21, lb_fails = 5, lb_late_fails = 1, lb_early_fails = 1)
)

analyses = analyses(
  commands = commands(
    samples_nopool = run_mcmc_generic(..dataset.., chains = 4, iter = 5000, warmup = 1500,
                                      stanmodel = "../stan_models/glfp_no_pooling.stan")
  ),
  datasets = datasets
)

targets = targets(datasets = datasets, analyses = analyses)

workflow(targets, sources = "../workflow/functions.R", packages = "rstan")

