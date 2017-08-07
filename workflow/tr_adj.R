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

samp <- extract(readRDS("../workflow/samples_lor_only3fails.rds"))
S <- length(samp$lp__)

source("../plotting_fns/KM_plot.R")

dat2 <- filter(dat, model %in% id)
tr_adj_df <- ddply(dat2, .(model), function(x){
  start = min(x$starttime)
  modelid = overview$stan_id[which(overview$model==x$model[1])]
  samp_tr = get_tr_adj(start, exp(as.numeric(samp$log_pi[,modelid])),
                       as.numeric(samp$mu1),
                       as.numeric(samp$sigma1),
                       as.numeric(samp$mu2[,modelid]), 
                       as.numeric(samp$sigma2[modelid]))
  adj = quantile(samp_tr, c(.25,.5,.75))
  data.frame(modelid=modelid, start=start, lower=adj[1], median=adj[2], upper=adj[3])
})

saveRDS(tr_adj_df, "../BB_data/tr_adj_tp2s2pi.rds")


dat3 <- filter(dat2, starttime >= .25*365*24)
tr_adj_df_quarteryr <- ddply(dat3, .(model), function(x){
  start = min(x$starttime)
  modelid = overview$stan_id[which(overview$model==x$model[1])]
  samp_tr = get_tr_adj(start, exp(as.numeric(samp$log_pi[,modelid])),
                       as.numeric(samp$mu1),
                       as.numeric(samp$sigma1),
                       as.numeric(samp$mu2[,modelid]), 
                       as.numeric(samp$sigma2[modelid]))
  adj = quantile(samp_tr, c(.25,.5,.75))
  data.frame(modelid=modelid, start=start, lower=adj[1], median=adj[2], upper=adj[3])
})
saveRDS(tr_adj_df_quarteryr, "../BB_data/tr_adj_fullmodel2.rds")
