# library(remakeGenerator)
# 
# datasets = commands(
#   data_8_1 = prepare_data(lb_fails = 8, lb_late_fails = 1, lb_early_fails = 0)
# )
# 
# analyses = analyses(
#   commands = commands(
#     samples = run_mcmc(..dataset.., chains = 4, iter = 50)
#   ),
#   datasets = datasets
# )
# 
# targets = targets(datasets = datasets, analyses = analyses)
# 
# workflow(targets, sources = "functions.R", packages = "rstan")

source("functions.R")
data_8_1 = prepare_data(lb_fails = 8, lb_late_fails = 1, lb_early_fails = 0)
samples = run_mcmc(data_8_1, chains = 4, iter = 50)
