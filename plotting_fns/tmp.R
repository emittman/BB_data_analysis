dat <- readRDS("../BB_data/clean_unit_summaries.rds")
dat$model <- as.integer(dat$model)

library(plyr)
library(dplyr)
overview <- ddply(dat, .(model), summarise,
                  n=length(model),
                  f=sum(failed>0),
                  early_f = sum(failed>0 & end_time<365*24*1),
                  late_f = sum(failed>0 & end_time>365*24*2))
# id <- unique(dat$model)
id <- with(overview, which(overview$late_f >= 1 & f >=8))
overview$stan_id <- NA
overview[id,]$stan_id <- 1:length(id)

s <- readRDS("../workflow/samples_10_7.rds")
samp <- extract(s)

plot(s, pars="pi")
pairs(s, pars=c("mu_pi","log_phi_pi", "lp__"))

filter(overview, stan_id==5)
filter(dat, model==id[5]) %>% ggplot(aes(x=start_time, y=end_time)) + geom_point()
pairs(s, pars=c("pi[5]", "mu1[5]","mu2[5]","log_sigma1[5]", "log_sigma2[5]"))
summary(s)$summary[c("mu1[5]","mu2[5]","log_sigma1[5]", "log_sigma2[5]", "pi[5]"),]

filter(overview, stan_id==14)
filter(dat, model==id[14]) %>% ggplot(aes(x=starttime, y=endtime)) + geom_point()
pairs(s, pars=c("pi[14]", "mu1[14]","mu2[14]","log_sigma1[14]", "log_sigma2[14]"))
summary(s)$summary[c("mu1[14]","mu2[14]","log_sigma1[14]", "log_sigma2[14]", "pi[14]"),]


source("../plotting_fns/KM_plot.R")
names(dat)[c(5,6)] <- c("endtime","starttime")
dat$censored <- dat$failed == 0

filter(dat, model==id[1]) %>%
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
  adj <- with(subset(mods, model == orig_id),
              get_tr_adj(min(starttime), pi_j, loc1_j, scl1_j, loc2_j, scl2_j))
  
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
                aes(x=x, y=y, ymin=lower, ymax=upper), fill="red", alpha=.2)
  
    
}