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
tr_adj <- readRDS("../BB_data/tr_adj_tp2s2pi.rds")$median


#ssigtp <- readRDS("MCMC_draws/vary_s2_and_tp2_4_17.rds")
#stp2 <- readRDS("MCMC_draws/samples_tp2_vary_new.rds")
#snull <- readRDS("MCMC_draws/samples_null_model_3_29.rds")
#sfull <- readRDS("MCMC_draws/samples_lor_only3fails.rds")



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
                     linetp = "longdash", colband="green", colline="darkblue")

for(j in 1:length(id)){
  orig_id <- id[j]
  pi_j = exp(sampfull$log_pi[max_id,j])
  loc1_j = sampfull$mu1[max_id]
  loc2_j = sampfull$mu2[max_id,j]
  scl1_j = sampfull$sigma1[max_id]
  scl2_j = sampfull$sigma2[max_id,j]
  adj <- tr_adj[j]   #New Median Truncation Values
  
  dat_tmp <- subset(mods, model == orig_id)
  
  #KM_list[[j]] <- KM_plot(data = dat_tmp, model = "weibull", tr_adj = adj, title="", linear_axes = FALSE, fixed= TRUE, xlimits=xlimits, ylimits = ylimits,verbose = F)
  
  KM_list[[j]] <- KM_plot_NP(data=dat_tmp, model = "weibull", tr_adj=adj, title = NULL, linear_axes = FALSE, fixed=TRUE,
                             xlimits=xlimits, ylimits=ylimits, verbose=F, conf=.90)
  tp2_list[[j]] <- KM_band(id=j, samp=samptp2, xlim=xlimits, ylim=ylimits, n_iter=100, quantiles=c(.05,.5,.95),
                           x_logscale=T, verbose=F, n = 100, pi.free=F, mu1.free=F, sigma1.free = F, mu2.free = T, sigma2.free = F,
                           colband="red", linetp="dotted", colline="black")
  s2tp_list[[j]] <- KM_band(id=j, samp=samps2tp, xlim=xlimits, ylim=ylimits, n_iter=100, quantiles=c(.05,.5,.95),
                            x_logscale=T, verbose=F, n = 100, pi.free=F, mu1.free=F, sigma1.free = F, mu2.free = T, sigma2.free = T,
                            colband="blue", linetp="dashed", colline="black")
  
  full_list[[j]] <- KM_band(id=j, samp=sampfull, xlim=xlimits, ylim=ylimits, n_iter=100, quantiles=c(.05,.5,.95),
                            x_logscale=T, verbose=F, n = 100, pi.free=T, mu1.free=F, sigma1.free = F, mu2.free = T, sigma2.free = T,
                            linetp="solid",colband="green")
}

j=26 #lack of fit on 4, compromise on 15, 37 hits the mark, all agree 26
KM_list[[j]] + #null_band[[2]] + null_band[[2]] +
  tp2_list[[j]][[1]] + tp2_list[[j]][[2]] +
  s2tp_list[[j]][[1]] + s2tp_list[[j]][[2]] + 
  full_list[[j]][[1]] + full_list[[j]][[2]] + null_band


#Models to Plot for Paper; Just show Bands for Full Model; Added LineType
#Possible Models to Show: 2, 9, 14, 23, 40, 45 

#Function to Grab 4 Drive Models and Make a List of the Plots
quad <- function(p1, p2, p3, p4){
  plts <- cbind(p1,p2,p3,p4)
  out <- list()
  for (i in 1:4){
    p1 <- plts[i]
    out[[i]] <- KM_list[[p1]] + null_list[[p1]][[2]] + tp2_list[[p1]][[2]] + s2tp_list[[p1]][[2]] + full_list[[p1]][[2]]
  }
  return((out))
}

#Test Function and Use Cow to Make Grid.
plot1 <- quad(2,9,14,40)
plot_grid(plot1[[1]], plot1[[2]], plot1[[3]], plot1[[4]], labels=c("2","9","14","40"), ncol = 2, nrow = 2)
  
# compare parameter estimates of log_tp2
results_list <- list(full = sfull, sig_tp = ssigtp, tp_only = stp2)
mean_df <- ldply(1:3, function(mod){
  means <- summary(results_list[[mod]])$summary[,"mean"]
  id <- grep("log_tp2", names(means))
  
  data.frame(model = names(results_list)[mod], means = means[id], dm = 1:44)
})

ggplot(mean_df, aes(x=means)) + geom_histogram(aes(y=..density..), bins=10) + facet_grid(model~.) + geom_density()


#compare pi's for full model (sampfull)
#Get Pi's Out with 95% Credible
#look at log_tp.05 quantile estimates
out.pi <- matrix(ncol=4, nrow=44)
for (i in 1:44){
  num=i
  seta <- paste0("log_pi[",num,"]",collapse="")
  etal <- summary(sfull)$summary[seta,"25%"]
  etam <- summary(sfull)$summary[seta,"50%"]
  etah <- summary(sfull)$summary[seta,"75%"]
  out.pi[i,1] <- exp(etal)
  out.pi[i,2] <- exp(etam)
  out.pi[i,3]<- exp(etah)
  out.pi[i,4]<- i
}

pi.dat <- as.data.frame(out.pi)
colnames(pi.dat) <- c("lower", "med","upper","model")

#Make Catepillar Plot
#Sort by Lower End Point Time
pi.dat$model <- factor(pi.dat$model,levels=pi.dat$model[order(pi.dat$lower)])  #Sort By Lower End Point


#Make Catepillar Plot for B10; Perhaps Sort by Sample Size?
p <- ggplot(pi.dat, aes(x=as.factor(model), y=med, ymin=lower, ymax=upper)) +
  geom_pointrange() +  
  coord_flip() + 
  xlab('Drive Model') + theme_bw() + 
  ylab(expression(pi))   +
  scale_x_discrete(breaks=seq(1,44,1))+
  scale_y_continuous(trans="logit", breaks=c(.01, .02, .05, .1, .2, .4, .6, .8, .9))
p





