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
id <- with(overview, which(overview$early >= 0 & overview$late_f >= 0 & f >=3))
overview$stan_id <- NA
overview[id,]$stan_id <- 1:length(id)

# s <- readRDS("../workflow/samples_lor_only3fails.rds")
ssigtp <- readRDS("../workflow_sig2_and_tp2_vary/vary_s2_and_tp2_4_17.rds")
stp2 <- readRDS("../workflow_tp2_vary/samples_tp2_vary.rds")
snull <- readRDS("../workflow_null/samples_null_model_3_29.rds")
sfull <- readRDS("../workflow/samples_lor_only3fails.rds")
samps2tp <- extract(ssigtp)
samptp2 <- extract(stp2)
sampnull <- extract(snull)
sampfull <- extract(sfull)
source("../plotting_fns/KM_plot.R")
source("../plotting_fns/greenwood_errors.R")
mods <- subset(dat, model %in% id)
max_id <- which.max(sampfull$lp__)
full_list <- list()
s2tp_list <- list()
tp2_list <-list()
KM_list <- list()
xlimits <- c(100,50000)
ylimits <- c(.0001, .9)
null_band <- KM_band(id=1, n_iter= 100, samp=sampnull, xlim=xlimits, ylim=ylimits, quantiles=c(.05,.5,.95),
                     x_logscale=T, verbose=F, n = 100, pi.free=F, mu1.free=F, sigma1.free = F, mu2.free = F, sigma2.free = F,
                     colband="green", colline="darkblue")
for(j in 1:length(id)){
  orig_id <- id[j]
  pi_j = exp(sampfull$log_pi[max_id,j])
  loc1_j = sampfull$mu1[max_id]
  loc2_j = sampfull$mu2[max_id,j]
  scl1_j = sampfull$sigma1[max_id]
  scl2_j = sampfull$sigma2[max_id,j]
  adj <- with(subset(mods, model == orig_id),get_tr_adj(min(starttime), pi_j, loc1_j, scl1_j, loc2_j, scl2_j))
  
  dat_tmp <- subset(mods, model == orig_id)
  
  # KM_list[[j]] <- KM_plot(data = dat_tmp, model = "weibull", tr_adj = adj, title="", linear_axes = TRUE, fixed= TRUE, xlimits=xlimits, ylimits = ylimits,verbose = F)
  KM_list[[j]] <- KM_plot_NP(data=dat_tmp, model = "weibull", tr_adj=adj, title = NULL, linear_axes = FALSE, fixed=TRUE,
                             xlimits=xlimits, ylimits=ylimits, verbose=F, conf=.90)
  tp2_list[[j]] <- KM_band(id=j, samp=samptp2, xlim=xlimits, ylim=ylimits, n_iter=100, quantiles=c(.05,.5,.95),
                           x_logscale=T, verbose=F, n = 100, pi.free=F, mu1.free=F, sigma1.free = F, mu2.free = T, sigma2.free = F,
                           colband="red", colline="black")
  s2tp_list[[j]] <- KM_band(id=j, samp=samps2tp, xlim=xlimits, ylim=ylimits, n_iter=100, quantiles=c(.05,.5,.95),
                            x_logscale=T, verbose=F, n = 100, pi.free=F, mu1.free=F, sigma1.free = F, mu2.free = T, sigma2.free = T,
                            colband="blue", colline="black")
  
  full_list[[j]] <- KM_band(id=j, samp=sampfull, xlim=xlimits, ylim=ylimits, n_iter=100, quantiles=c(.05,.5,.95),
                            x_logscale=T, verbose=F, n = 100, pi.free=T, mu1.free=F, sigma1.free = F, mu2.free = T, sigma2.free = T,
                            colband="green")
}

j=26 #lack of fit on 4, compromise on 15, 37 hits the mark, all agree 26
KM_list[[j]] + #null_band[[1]] + null_band[[2]] +
  tp2_list[[j]][[1]] + tp2_list[[j]][[2]] +
  s2tp_list[[j]][[1]] + s2tp_list[[j]][[2]] + 
  full_list[[j]][[1]] + full_list[[j]][[2]] + null_band

# compare parameter estimates of log_tp2
results_list <- list(full = sfull, sig_tp = ssigtp, tp_only = stp2)
mean_df <- ldply(1:3, function(mod){
  means <- summary(results_list[[mod]])$summary[,"mean"]
  id <- grep("log_tp2", names(means))
  
  data.frame(model = names(results_list)[mod], means = means[id], dm = 1:44)
})

ggplot(mean_df, aes(x=means)) + geom_histogram(aes(y=..density..), bins=10) + facet_grid(model~.) + geom_density()

