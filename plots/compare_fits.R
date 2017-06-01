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
ssigtp <- readRDS("../workflow_sig2_and_tp2_vary/vary_s2_and_tp2_4_17.rds")
stp2 <- readRDS("../workflow_tp2_vary/samples_tp2_vary.rds")
snull <- readRDS("../workflow_null/samples_null_model_3_29.rds")
sfull <- readRDS("../workflow/samples_lor_only3fails.rds")
tr_adj <- readRDS("../BB_data/tr_adj_tp2s2pi.rds")$median


ssigtp <- readRDS("../MCMC_draws/vary_s2_and_tp2_4_17.rds")
stp2 <- readRDS("../MCMC_draws/samples_tp2_vary_new.rds")
snull <- readRDS("../MCMC_draws/samples_null_model_3_29.rds")
sfull <- readRDS("../MCMC_draws/samples_lor_only3fails.rds")



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

null_band <- KM_band("Model 1", id=1, n_iter= 100, samp=sampnull, xlim=xlimits, ylim=ylimits, quantiles=c(.05,.5,.95),
                     x_logscale=T, verbose=F, n = 100, pi.free=F, mu1.free=F, sigma1.free = F, mu2.free = F, sigma2.free = F,
                     linetp = "longdash", colband="green")

#for(j in 1:length(id)){
for(j in 1:44){
  orig_id <- id[j]
  pi_j = exp(sampfull$log_pi[max_id,j])
  loc1_j = sampfull$mu1[max_id]
  loc2_j = sampfull$mu2[max_id,j]
  scl1_j = sampfull$sigma1[max_id]
  scl2_j = sampfull$sigma2[max_id,j]
  adj <- tr_adj[j]   #New Median Truncation Values
  
  dat_tmp <- subset(mods, model == orig_id)
  
  KM_list[[j]] <- KM_plot(data = dat_tmp, model = "weibull", tr_adj = adj, title="", linear_axes = FALSE, fixed= TRUE, xlimits=xlimits, 
                          size=.2, ylimits = ylimits,verbose = F)
  
  #KM_list[[j]] <- KM_plot_NP(data=dat_tmp, model = "weibull", tr_adj=adj, title = NULL, linear_axes = FALSE, fixed=TRUE,
                            # xlimits=xlimits, ylimits=ylimits, verbose=F, conf=.90)
  
  tp2_list[[j]] <- KM_band("Model 2", id=j, samp=samptp2, xlim=xlimits, ylim=ylimits, n_iter=100, quantiles=c(.05,.5,.95),
                           x_logscale=T, verbose=F, n = 100, pi.free=F, mu1.free=F, sigma1.free = F, mu2.free = T, sigma2.free = F,
                           colband="red", linetp="dotted")
  
  s2tp_list[[j]] <- KM_band("Model 3", id=j, samp=samps2tp, xlim=xlimits, ylim=ylimits, n_iter=100, quantiles=c(.05,.5,.95),
                            x_logscale=T, verbose=F, n = 100, pi.free=F, mu1.free=F, sigma1.free = F, mu2.free = T, sigma2.free = T,
                            colband="blue", linetp="dashed")
  
  full_list[[j]] <- KM_band("Model 4",id=j, samp=sampfull, xlim=xlimits, ylim=ylimits, n_iter=100, quantiles=c(.05,.5,.95),
                            x_logscale=T, verbose=F, n = 100, pi.free=T, mu1.free=F, sigma1.free = F, mu2.free = T, sigma2.free = T,
                            linetp="solid")
}

j=26 #lack of fit on 4, compromise on 15, 37 hits the mark, all agree 26
KM_list[[j]] + null_band[[2]] + null_band[[2]] +
  tp2_list[[j]][[1]] + tp2_list[[j]][[2]] +
  s2tp_list[[j]][[1]] + s2tp_list[[j]][[2]] + 
  full_list[[j]][[1]] + full_list[[j]][[2]] + null_band

#Pick 4 Drive Models for Paper.  

#Make One Plot with Legend:
j=2
ptest = KM_list[[j]] + null_band[[2]] + tp2_list[[j]][[2]] + s2tp_list[[j]][[2]] + full_list[[j]][[2]]
p1 <- ptest + scale_colour_manual(values=c("green", "red","blue","orange"), name="Model") + 
  scale_linetype_manual(values=c("longdash","dotted","dashed","solid") , name="Model") + 
  theme(legend.key.width=unit(3,"line")) 

j=9
ptest = KM_list[[j]] + null_band[[2]] + tp2_list[[j]][[2]] + s2tp_list[[j]][[2]] + full_list[[j]][[2]]
p2 <- ptest + scale_colour_manual(values=c("green", "red","blue","orange"), name="Model") + 
  scale_linetype_manual(values=c("longdash","dotted","dashed","solid") , name="Model") + 
  theme(legend.key.width=unit(3,"line")) 

p2 <- p2 + theme(legend.position = "none")

j=14
ptest = KM_list[[j]] + null_band[[2]] + tp2_list[[j]][[2]] + s2tp_list[[j]][[2]] + full_list[[j]][[2]]
p3 <- ptest + scale_colour_manual(values=c("green", "red","blue","orange"), name="Model") + 
  scale_linetype_manual(values=c("longdash","dotted","dashed","solid") , name="Model") + 
  theme(legend.key.width=unit(3,"line")) 
p3 <- p3 + theme(legend.position = "none")

j=40
ptest = KM_list[[j]] + null_band[[2]] + tp2_list[[j]][[2]] + s2tp_list[[j]][[2]] + full_list[[j]][[2]]
p4 <- ptest + scale_colour_manual(values=c("green", "red","blue","orange"), name="Model") + 
  scale_linetype_manual(values=c("longdash","dotted","dashed","solid") , name="Model") + 
  theme(legend.key.width=unit(3,"line")) 
p4 <- p4 + theme(legend.position = "none")


plot_grid(p1, p2, p3, p4, labels=c("2","9","14","40"), ncol = 2, nrow = 2)
  
# compare parameter estimates of log_tp2
results_list <- list(full = sfull, sig_tp = ssigtp, tp_only = stp2)
mean_df <- ldply(1:3, function(mod){
  means <- summary(results_list[[mod]])$summary[,"mean"]
  id <- grep("log_tp2", names(means))
  
  data.frame(model = names(results_list)[mod], means = means[id], dm = 1:44)
})

ggplot(mean_df, aes(x=means)) + geom_histogram(aes(y=..density..), bins=10) + facet_grid(model~.) + geom_density()



#Plot B10
#Note: 50% CI for B10
b10 <- readRDS("../paper/b10_50ci.rds")
colnames(b10) <- c("model","lb","med","ub")
b10$lb <- b10$lb/(24*365)
b10$med <- b10$med/(24*365)
b10$ub <- b10$ub/(24*365)

#Sort by Lower End Point Time
b10$model <- factor(b10$model,levels=b10$model[order(b10$med)])  #Sort By Lower End Point


#Make Catepillar Plot for B10; Perhaps Sort by Sample Size?
b10.plot <- ggplot(b10, aes(x=model, y=med, ymin=lb, ymax=ub)) +
  geom_pointrange() +  
  coord_cartesian(ylim = c(0, 8)) + 
  xlab('Drive Model') + theme_bw() + 
  ylab("Years")   +
  scale_x_discrete(breaks=seq(1,44,1))
b10.plot





#compare pi's for full model (sampfull)
#Get Pi's Out with 50% Credible
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
#Sort by Median of B10

pi.dat$bmed <- b10$med
pi.dat$model <- factor(pi.dat$model,levels=pi.dat$model[order(pi.dat$bmed)])  #Sort By Lower End Point


#Make Catepillar Plot for B10; Perhaps Sort by Sample Size?
pi.plot <- ggplot(pi.dat, aes(x=as.factor(model), y=med, ymin=lower, ymax=upper)) +
  geom_pointrange() +  
  xlab('Drive Model') + theme_bw() + 
  ylab(expression(pi))   +
  scale_x_discrete(breaks=seq(1,44,1))+
  scale_y_continuous(trans="logit", breaks=c(.01, .02, .05, .1, .2, .4, .6, .8, .9))
pi.plot

#Make Cow Plot
plot_grid(b10.plot, pi.plot, ncol = 2, nrow = 1)



