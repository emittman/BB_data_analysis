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
id <- with(overview, which(overview$late_f >= 1 & f >=8))
overview$stan_id <- NA
overview[id,]$stan_id <- 1:length(id)

s <- readRDS("../workflow/samples1013.rds")
samp <- extract(s)

plot(s, pars="logit_pi_raw")
pairs(s, pars=c("eta_pi","tau_pi", "lp__"))

plot(s, pars="log_sigma1")
pairs(s, pars=c("eta_ls1", "tau_ls1", "lp__"))

plot(s, pars="log_tp1_raw")
pairs(s, pars=c("eta_ltp1","tau_ltp1", "lp__"))

plot(s, pars="log_sigma2")
pairs(s, pars=c("eta_ls2", "tau_ls2", "lp__"))

plot(s, pars="mu2")
pairs(s, pars=c("eta_ltp2","tau_ltp2", "lp__"))



filter(overview, stan_id==5)
filter(dat, model==id[5]) %>% ggplot(aes(x=starttime, y=endtime)) + geom_point()
pairs(s, pars=c("logit_pi_raw[5]", "log_tp1_raw[5]","log_tp2_raw[5]","log_sigma1_raw[5]", "log_sigma2_raw[5]"))
summary(s)$summary[c("mu1[5]","mu2[5]","log_sigma1[5]", "log_sigma2[5]", "log_pi[5]"),]

filter(overview, stan_id==14)
filter(dat, model==id[14]) %>% ggplot(aes(x=starttime, y=endtime)) + geom_point()
pairs(s, pars=c("logit_pi_raw[14]", "log_tp1_raw[14]","log_tp2_raw[14]","log_sigma1_raw[14]", "log_sigma2_raw[14]"))
summary(s)$summary[c("mu1[14]","mu2[14]","log_sigma1[14]", "log_sigma2[14]", "log_pi[14]"),]


source("../plotting_fns/KM_plot.R")

filter(dat, model==id[6]) %>%
  KM_plot("weibull")

mods <- subset(dat, model %in% id)
max_id <- which.max(samp$lp__)
plot_list <-NULL
for(j in 1:length(id)){
  orig_id <- id[j]
  pi_j = samp$pi[max_id,j]
  loc1_j = samp$mu1[max_id,j]
  loc2_j = samp$mu2[max_id,j]
  scl1_j = exp(samp$log_sigma1[max_id,j])
  scl2_j = exp(samp$log_sigma2[max_id,j])
  adj <- with(subset(mods, model == orig_id),get_tr_adj(min(starttime), pi_j, loc1_j, scl1_j, loc2_j, scl2_j))
  
  kmp <- subset(mods, model == orig_id) %>%
    KM_plot(model="weibull", tr_adj = adj, fixed=T)
  band_df <- data.frame(x = exp(seq.int(log(1000),
                                             log(50000),
                                             length.out=25))) %>%
    ddply(.(x), function(g){
      Fp <- sapply(1:1000, function(i) {
        1 -  (1 - samp$pi[i,j] * my_pweibull(g$x, samp$mu1[i,j], exp(samp$log_sigma1[i,j]))) *
          (1 - my_pweibull(g$x, samp$mu2[i,j], exp(samp$log_sigma2[i,j])))
      })
      q <- quantile(Fp, c(.05, .5, .95))
      data.frame(y = q[2], lower = max(.0001,q[1]), upper = q[3])
    })
  
  plot_list[[j]] <- kmp + geom_line(data = band_df, inherit.aes = FALSE, aes(x, y), lty=2) +
    geom_ribbon(data = band_df, inherit.aes = FALSE,
                aes(x=x, y=y, ymin=lower, ymax=upper), fill="red", alpha=.2) #+
    #geom_ribbon(data=check,aes(x=time,y=g,ymin=lower,ymax=upper),fill="green") #Add MLE Wald Bands; look funny though for mod 6
  
    
}