#setwd("workflow/")
dat <- readRDS("../BB_data/clean_unit_summaries.rds")
dat$model <- as.integer(dat$model)

library(plyr)
library(dplyr)
library(rstan)
overview <- ddply(dat, .(model), summarise,
                  n=length(model),
                  f=sum(failed>0),
                  early_f = sum(failed>0 & endtime<365*24*1),
                  late_f = sum(failed>0 & endtime>365*24*2))
# id <- unique(dat$model)
id <- with(overview, which(overview$early >= 1 & overview$late_f >= 1 & f >=5))
overview$stan_id <- NA
overview[id,]$stan_id <- 1:length(id)

s <- readRDS("../workflow/reduced_relaxed_21dms.rds")
samp <- extract(s)

source("../plotting_fns/KM_plot.R")

mods <- subset(dat, model %in% id)
max_id <- which.max(samp$lp__)
plot_list <-NULL
for(j in 1:length(id)){
  orig_id <- id[j]
  pi_j = exp(samp$log_pi[max_id,j])
  loc1_j = samp$mu1[max_id]
  loc2_j = samp$mu2[max_id,j]
  scl1_j = samp$sigma1[max_id]
  scl2_j = samp$sigma2[max_id,j]
  adj <- with(subset(mods, model == orig_id),get_tr_adj(min(starttime), pi_j, loc1_j, scl1_j, loc2_j, scl2_j))
  
  dat_tmp <- subset(mods, model == orig_id)

  kmb <- KM_with_band(data = dat_tmp, id =j, samp = samp, n_iter = 200, n = 50, quantiles = c(.05,.5,.95), tr_adj = adj,
                      xlimits = c(0,50000), ylimits = c(0,.9), fixed = T, linear_axes = T, verbose = F, model = "weibull")
  
  plot_list[[j]] <- kmb
}

gridExtra::grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]],
                        plot_list[[6]], plot_list[[7]], plot_list[[8]], plot_list[[9]], plot_list[[10]],
                        plot_list[[11]], plot_list[[12]], plot_list[[13]], plot_list[[14]], plot_list[[15]],
                        plot_list[[16]], plot_list[[17]], plot_list[[18]], plot_list[[19]], plot_list[[20]], 
                        plot_list[[21]], nrow=5, ncol=5)
