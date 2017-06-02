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

s <- readRDS("../workflow/samples_lor_only3fails.rds")
samp <- extract(s)

source("../plotting_fns/KM_plot.R")

mods <- subset(dat, model %in% id)
max_id <- which.max(samp$lp__)
plot_list <-NULL

adjs <-readRDS("../BB_data/tr_adj_tp2s2pi.rds")$median

selected <- c(4,5,2,3,6,8,26,23)

for(i in 1:length(selected)){
  j <- selected[i]
  orig_id <- id[j]
  adj <- adjs[j]
  xlim <- c(100,50000)
  ylim <- c(.005,.5)
  dat_tmp <- subset(mods, model == orig_id)
  p <- KM_plot_NP(dat_tmp, "weibull", adj, as.character(j),FALSE, TRUE, xlimits=xlim, ylimits=ylim, conf = .9)
  p <- p + KM_band(num=as.character(j), id=j, samp=samp, n_iter = 1000, xlim = xlim, ylim = ylim, n = 100,
                   pi.free=T, mu1.free=F, mu2.free=T,sigma1.free=F,sigma2.free=T)+
    scale_y_continuous(trans="qsev", breaks=c(.01, .05, .15, .4), limits=c(.005,.5))+
    theme(legend.position="none",
          axis.title = element_blank())+
    scale_x_continuous(trans="log", breaks=c(.01, .1, 1, 5)*365*24, 
                       labels=as.character(c(.01,.1,1,5)), limits=xlim)
  
  plot_list[[i]] <- p
}

plot_grid(plot_list[[1]],plot_list[[2]], plot_list[[3]],plot_list[[4]],
          plot_list[[5]], plot_list[[6]], plot_list[[7]], plot_list[[8]], nrow=2)
