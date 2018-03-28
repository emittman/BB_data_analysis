library(drake)

dir.create("dbl_sd_sampling_output")

my_plan <- plan(
  data_3_25 = prepare_data(lb_fails = 3, lb_late_fails = 0, lb_early_fails = 0),
  samples_double_variance_3_25 = run_mcmc_generic(dataset = data_3_25, chains = 16,
    iter = 3000, warmup = 1500,
    stanmodel = "../stan_models/glfp_lo_reduced_doublyrelaxed.stan"),
  strings_in_dots = "literals"
)

check(my_plan)
# plot_graph(my_plan)

source('functions.R')

make(my_plan, packages = "rstan")