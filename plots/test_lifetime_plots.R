library(dplyr)
library(ggplot2)

#load data
df <- readRDS("BB_data/clean_unit_summaries.rds")
overvw <- readRDS("BB_data/overview.rds")
overvw$stan_id <- NA
overvw$stan_id[which(overvw$f>=3)] <- 1:44

#model 14
source("plotting_fns/lifetime_plot.R")

filter(df, model==overvw$model[which(overvw$stan_id==14)]) %>%
  lifetime_plot(n_to_show=800, lab="Drive-model 14", trans="log",
                xlimits=c(50,40000), xlabels=c(.01,.05,.25,1,2,4))

filter(df, model==overvw$model[which(overvw$stan_id==14)]) %>%
  lifetime_plot2(n_to_show=800, lab="Drive-model 14", trans="log",
                xlimits=c(50,40000), xlabels=c(.01,.05,.25,1,2,4))
