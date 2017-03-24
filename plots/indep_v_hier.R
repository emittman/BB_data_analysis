library(backblaze)
library(plyr)
library(dplyr)
library(rstan)

source("../workflow/functions.R")
source("../plotting_fns/KM_plot.R")
source("../plotting_fns/greenwood_errors.R")

dat <- backblaze
dat$model <- as.integer(dat$model)
dat$censored <- !dat$failed
names(dat)[5:6] <- c("starttime","endtime")
overview <- ddply(dat, .(model), summarise,
                  n=length(model),
                  f=sum(failed>0),
                  early_f = sum(failed>0 & endtime<365*24*1),
                  late_f = sum(failed>0 & endtime>365*24*2))

# id <- unique(dat$model)
id <- with(overview, which(overview$early >= 1 & overview$late_f >= 1 & f >=5))
overview$stan_id <- NA
overview[id,]$stan_id <- 1:length(id)

# hierarchical model with informative priors
s <- readRDS("../workflow/samples_2_1.rds")
samp <- extract(s)

# hierarchical reduced model with relaxed priors
id2 <- with(overview, which(overview$early >= 0 & overview$late_f >= 0 & f >=3))
s_red <- readRDS("../workflow/samples_lor_only3fails.rds")
samp_red <- extract(s_red)

plot_count <- 0
plot_list <- NULL
plot_list2 <- NULL # compare credible intervals for pi
for(i in 1:21){
  #unpooled inference
  stemp  <- readRDS(paste(c("../workflow2/samples_nopool_data_dm",i,".rds"),collapse=""))
  #if divergent iterations present, go to next drive model
  if(count_divergences(stemp, 1:4, inc_warmup = FALSE)>0) next
  plot_count <- plot_count+1
  mod_id <- overview$model[which(overview$stan_id==i)]
  dat_tmp <- subset(dat, model==mod_id)
  samp_tmp <- extract(stemp)
  mxid_tmp <- which.max(samp_tmp$lp_)
  map_tmp <- with(samp_tmp, data.frame(p=exp(log_pi[mxid_tmp]), m1=mu1[mxid_tmp], m2=mu2[mxid_tmp], s1=sigma1[mxid_tmp], s2=sigma2[mxid_tmp]))
  adj_tmp <- with(map_tmp,get_tr_adj(min(dat_tmp$starttime), p, m1, s1, m2, s2))
  KM_np <- KM_plot_NP(data=dat_tmp, "weibull", tr_adj=adj_tmp, "Kaplan Meier with Greenwood Standard Errors",
                      linear_axes = F, fixed=T, xlimits = c(200, 6e4), ylimits = c(0.001,.8),
                      conf=.9, verbose=F)+theme(axis.text.x = element_text(angle = 60, hjust = 1))
  KM_up <- KM_with_band(paste(c("Unpooled"), collapse=""),
                        data=dat_tmp, id=1, samp=samp_tmp, n_iter=500, n=40,
                        quantiles = c(.1,.5,.9), tr_adj=adj_tmp, xlimits = c(200, 6e4),
                        ylimits = c(.001,.8), fixed=T, linear_axes = F, verbose = F,
                        model = "weibull")+theme(axis.text.x = element_text(angle = 60, hjust = 1))
  
  KM_pl <- KM_with_band("Pooled", data=dat_tmp, id=i, samp=samp, n_iter=500, n=40,
                        quantiles = c(.1,.5,.9), tr_adj=adj_tmp, xlimits = c(200, 6e4),
                        ylimits = c(.001,.8), fixed=T, linear_axes = F, verbose = F,
                        model = "weibull")+theme(axis.text.x = element_text(angle = 60, hjust = 1))
  
  KM_red <- KM_with_band("Pooled reduced", data=dat_tmp, id=which(id2 == id[i]), samp=samp_red, n_iter=500, n=40,
                        quantiles = c(.1,.5,.9), tr_adj=adj_tmp, xlimits = c(200, 6e4),
                        ylimits = c(.001,.8), fixed=T, linear_axes = F, verbose = F,
                        model = "weibull")+theme(axis.text.x = element_text(angle = 60, hjust = 1))
  
  plot_list[[plot_count]] <- gridExtra::grid.arrange(KM_np, KM_up, KM_pl, KM_red, ncol=2, top=paste(c("Drive model", i), collapse=""))
}

pdf("compare_converged_models.pdf")
for(p in plot_list){
  plot(p)
}
dev.off()
