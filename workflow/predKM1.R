setwd("workflow/")
dat <- readRDS("../BB_data/clean_unit_summaries.rds")
dat$model <- as.integer(dat$model)

library(plyr)
library(dplyr)
library(rstan)
library(ggplot2)
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

sfull <- readRDS("../workflow/samples_lor_only3fails.rds")
tr_adj <- readRDS("../BB_data/tr_adj_tp2s2pi.rds")$median

sampfull <- extract(sfull)
source("../plotting_fns/KM_plot2.R")
source("../predictive/sample_from_predictive.R")


make_replicate_data <- function(orig, mcmc_samp, dm, samp_id){
  n <- nrow(orig)
  id_obs <- which(orig$censored==0)
  id_cens <- which(!(1:n %in% id_obs))
  t_tr <- orig$starttime
  t_cens <- rep(max(orig$endtime), n)
  t_cens[id_cens] <- orig$endtime[id_cens]
  
  replicates <- lapply(samp_id, function(s){
  y <- with(mcmc_samp, sample_glfp_rep(n, p=exp(log_pi[s,dm]),
                                      m1=mu1[s], s1=sigma1[s],
                                      m2=mu2[s,dm], s2=sigma2[s,dm]))
  
    is_cens <- y>t_cens
    y[is_cens] <- t_cens[is_cens]
  
    id_ok <- which(y>t_tr)
    n_new <- length(id_ok)
    print(cat(c("length ok:", n_new, "\n")))
    print(cat(c("number failed:", sum(!is_cens[id_ok]), "\n")))
    if(sum(!is_cens[id_ok]) == 0) return(NULL)
    out <-KM.survfit(data.frame(starttime=t_tr[id_ok], endtime=y[id_ok],
                          censored=is_cens[id_ok]), alpha = .1)
    print(cat(c("afterfit\n\n")))
    out
  })
  null_id <- which(sapply(replicates, is.null))
  if(length(null_id)>0){replicates <- replicates[-null_id]}
  replicates
}

get_plotlist <- function(dms){
  plist <- list()
  n <- length(dms)
  for(i in 1:n){
  dm <- dms[i]
  dat_dm <- filter(dat, model==overview$model[which(overview$stan_id==dm)])
  
  
  
  reps <- 19
  samp_id <- sample(length(sampfull$mu1), reps)
  replicates <- make_replicate_data(orig = dat_dm, mcmc_samp = sampfull, dm=dm, samp_id = samp_id)
  
  
  xlimits <- c(100,50000)
  ylimits <- c(.001, .9)
  
  bp <- baseKMplot.multiple(replicates, xlimits=xlimits, ylimits=ylimits,
                            color="black",linetype = "dashed", logscale = TRUE,
                            label="replicates", alpha=.4)
  
  km_obs <- KM.survfit(dat_dm)
  bp <- addKmToBaseplot(bp, km_obs, color="red",linetype = "solid", label = "observed")
  plist[[i]] <- plotFinally(bp, xbrks = c(100, 1000, 5000, 10000, 20000, 40000),
              ybrks = c(.001,.01,.1,.25,.5,.75,.9), years=TRUE) +
    theme_bw(base_size=14)+
    theme(panel.grid=element_blank(),
          legend.position = "none",
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) + ggtitle(as.character(dm))
  }
  plist
}
set.seed(10438)
ran_set <- sort(sample(44, 6))
plist <- get_plotlist(ran_set)

library(cowplot)
# length(plist)
pdf("../paper/fig/pospredictKM.pdf")
plot_grid(plist[[1]], plist[[2]],plist[[3]],plist[[4]],plist[[5]],plist[[6]], ncol=3)
dev.off()
# 
# plist2 <- get_plotlist(7:12)
# plot_grid(plist2[[7]], plist2[[8]],plist2[[9]],plist2[[10]],plist2[[11]],plist2[[12]], ncol=3)
# 
# plist3 <- get_plotlist(13:18)
# plot_grid(plist3[[13]], plist3[[14]], plist3[[15]], plist3[[16]], plist3[[17]], plist3[[18]], ncol=3)
# 
# plist4 <- get_plotlist(19:24)
# plot_grid(plist4[[19]], plist4[[20]], plist4[[21]], plist4[[22]], plist4[[23]], plist4[[24]], ncol=3)
# 
# plist5 <- get_plotlist(25:30)
# plot_grid(plist5[[25]], plist5[[26]], plist5[[27]], plist5[[28]], plist5[[29]], plist5[[30]], ncol=3)
# 
# plist6 <- get_plotlist(c(31:33,35:36))
# plot_grid(plist6[[31]], plist6[[32]], plist6[[33]], plist6[[35]], plist6[[36]], ncol=3)
# 
# plist7 <- get_plotlist(c(37:39))
# plist8 <- get_plotlist(c(40:41))
# plist9 <- get_plotlist(c(44))
# plot_grid(plist7[[37]],plist7[[38]],plist7[[39]],plist8[[40]],plist8[[41]],plist9[[44]], ncol=3)
