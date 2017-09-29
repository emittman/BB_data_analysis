###########################
#' Compare KM plots for all drives simulaneously
#' to show differences between models
#'
#' Still TO DO:
#' linetype the drives to pick out individual lines across models
#' reconsider presentation (make 2 plots?)





# setwd("workflow/")
dat <- readRDS("../BB_data/clean_unit_summaries.rds")
dat$model <- as.integer(dat$model)

tr_adj <- readRDS("../BB_data/tr_adj_tp2s2pi.rds")$median

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

df1 <- dat %>% filter(model %in% id) %>%
  select(-start_date, -end_date) %>%
  mutate(model = overview$stan_id[match(model, overview$model)])

source("../plotting_fns/KM_plot2.R")

xlimits <- c(5, 50000)
ylimits <- c(0, .9)
x.breaks <- 0:5*1e4
y.breaks <- 0:4*.25

all_kms <- ldply(unique(df1$model), function(m){
  fit <- KM.survfit(filter(df1, model==m))
  fit <- tr.adj(fit, tr_adj = tr_adj[m])
  fit <- trimKMx(fit,xlimits)
  fit <- trimKMy(fit,ylimits)
  data.frame(model=m,fit)})

p <- ggplot(all_kms, aes(x=t, y=Ft, color=factor(model), group=factor(model))) + geom_step(size=1) +
  theme_bw(base_size = 14) + scale_colour_grey(start = 0, end = .9) +
  scale_x_continuous(name="Thousands of hours", breaks = x.breaks, labels=x.breaks/1000,
                     limits = xlimits) +
  scale_y_continuous(name="Fraction failing", breaks=y.breaks, limits=ylimits)+
  theme(legend.position = "none") + 
  ggtitle("Adjusted Kaplan-Meier")

null_s <- readRDS("../workflow_null/samples_null_model_3_29.rds")
null_samp <- with(extract(null_s), list(mu1=mu1, sigma1=sigma1, mu2=mu2, sigma2=sigma2, log_pi=log_pi)) 

b1 <- bandFromPSamp(null_samp, range=xlimits, length.out = 30, N = 1000, logscale = FALSE)
pnull <- ggplot(b1$band, aes(x=t, y=est)) + geom_line(size=1) +
  scale_x_continuous(breaks = x.breaks, labels=x.breaks/1000)+
  scale_y_continuous(breaks=y.breaks, limits=ylimits)+
  theme_bw(base_size = 14) + ggtitle("Model 1") +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

mod2_s <- readRDS("../workflow_tp2_vary/samples_tp2_vary.rds")

all_mod2 <- ldply(unique(df1$model), function(m){
  samp <- with(extract(mod2_s), list(mu1=mu1, sigma1=sigma1, mu2=mu2[,m], sigma2=sigma2, log_pi=log_pi)) 
  fit <- bandFromPSamp(samp, range=xlimits, length.out = 30, N = 1000, logscale=FALSE)
  fit <- bandTrimx(fit, xlimits)
  data.frame(model=m,fit$band)})

p_mod2 <- ggplot(all_mod2, aes(x=t,y=est, group=model)) +
  geom_line(size=1, alpha=.5) +
  scale_x_continuous(name="", breaks = x.breaks, labels=x.breaks/1000)+
  scale_y_continuous(name="", breaks=y.breaks, limits=ylimits)+
  theme_bw(base_size = 14) + ggtitle("Model 2")+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

mod3_s <- readRDS("../workflow_sig2_and_tp2_vary/vary_s2_and_tp2_4_17.rds")

all_mod3 <- ldply(unique(df1$model), function(m){
  samp <- with(extract(mod3_s), list(mu1=mu1, sigma1=sigma1, mu2=mu2[,m], sigma2=sigma2[,m], log_pi=log_pi)) 
  fit <- bandFromPSamp(samp, range=xlimits, length.out = 30, N = 1000, logscale=FALSE)
  fit <- bandTrimx(fit, xlimits)
  data.frame(model=m,fit$band)})

p_mod3 <- ggplot(all_mod3, aes(x=t,y=est, group=model)) +
  geom_line(size=1, alpha=.5) +
  scale_x_continuous(name="", breaks = x.breaks, labels=x.breaks/1000)+
  scale_y_continuous(name="", breaks=y.breaks, limits=ylimits)+
  theme_bw(base_size = 14) + ggtitle("Model 3")+
  theme(axis.title = element_blank())

mod4_s <- readRDS("../workflow/samples_lor_only3fails.rds")

all_mod4 <- ldply(unique(df1$model), function(m){
  samp <- with(extract(mod4_s), list(mu1=mu1, sigma1=sigma1, mu2=mu2[,m], sigma2=sigma2[,m], log_pi=log_pi[,m])) 
  fit <- bandFromPSamp(samp, range=xlimits, length.out = 30, N = 1000, logscale=FALSE)
  fit <- bandTrimx(fit, xlimits)
  data.frame(model=m,fit$band)})

p_mod4 <- ggplot(all_mod4, aes(x=t,y=est, group=model)) +
  geom_line(size=1, alpha=.5) +
  scale_x_continuous(name="", breaks = x.breaks, labels=x.breaks/1000)+
  scale_y_continuous(name="", breaks=y.breaks, limits=ylimits)+
  theme_bw(base_size = 14) + ggtitle("Model 4")+
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())


cowplot::plot_grid(p, cowplot::plot_grid(pnull,p_mod2,
                                p_mod3,p_mod4, nrow=2, align="hv"),
                   nrow=1, align="h")

# tp_s <- readRDS("../workflow_tp2_vary/samples_tp2_vary_new.rds")
# tp_samp <- extract(tp_s)
# tp_layers <- lapply(1:length(id), function(j){
#   KM_band(j, tp_samp, n_iter=1000, xlim=xlimits, ylim=ylimits,
#           pi.free=F, mu2.free=T, sigma2.free=F, x_logscale = F)[[2]]
# })
# p3 <- ggplot(data.frame(x=0,y=0), aes(x,y))
# for(j in 1:length(id)){
#   p3 <- p3 + tp_layers[[j]]
# }
# p3 <- p3 + theme_bw(base_size = 14)+xlab("Hours")+ylab("Fraction Failing")+
#   xlim(xlimits)+ylim(ylimits)
# 
# cowplot::plot_grid(p1 + no_legend,
#                    p2,
#                    p3, ncol=1)
#   
# tpsig_s <- readRDS("../workflow_sig2_and_tp2_vary/vary_s2_and_tp2_4_17.rds")
# tpsig_samp <- extract(tpsig_s)
# tpsig_layers <- lapply(1:length(id), function(j){
#   KM_band(j, tpsig_samp, n_iter=1000, xlim=xlimits, ylim=ylimits,
#           pi.free=F, mu2.free=T, sigma2.free=T, x_logscale = F)[[2]]
# })
# p4 <- ggplot(data.frame(x=0,y=0), aes(x,y))
# for(j in 1:length(id)){
#   p4 <- p4 + tpsig_layers[[j]]
# }
# p4 <- p4 + theme_bw(base_size = 14)+xlab("Hours")+ylab("Fraction Failing")+
#   xlim(xlimits)+ylim(ylimits)
# 
# cowplot::plot_grid(p1 + no_legend,
#                    p2,
#                    p3,
#                    p4, ncol=2)
# 
# tpsigpi <- readRDS("../workflow/samples_lor_only3fails.rds")
# tpsigpi_samp <- extract(tpsigpi)
# tpsigpi_layers <- lapply(1:length(id), function(j){
#   KM_band(j, tpsigpi_samp, n_iter=1000, xlim=xlimits, ylim=ylimits,
#           pi.free=T, mu2.free=T, sigma2.free=T, x_logscale = F)[[2]]
# })
# 
# p5 <- ggplot(data.frame(x=0,y=0), aes(x,y))
# for(j in 1:length(id)){
#   p5 <- p5 + tpsigpi_layers[[j]]
# }
# p5 <- p5 + theme_bw(base_size = 14)+xlab("Hours")+ylab("Fraction Failing")+
#   xlim(xlimits)+ylim(ylimits)
# 
# #decided to reorganize this (5-24-2017)
# # cowplot::plot_grid(p1 + no_legend, cowplot::plot_grid(p2,p3,p4,p5, ncol=2), ncol=2)
# 
# bytimeobs <- ddply(df1, .(model), function(x){
#   total_time = sum(x$endtime - x$starttime)
#   failures = sum(x$failed>0)
#   data.frame(total_time=total_time, failures=failures)
# })
# 
# p6 <- ggplot(bytimeobs, aes(factor(model), total_time)) + geom_bar(stat="identity") +
#   scale_y_continuous(trans="log", breaks=c(100,10000,1000000,1000000))+
#   ggtitle("Total operational hours")+ylab("")+xlab("")+
#   theme(axis.ticks.x = element_blank(),
#         axis.text.x = element_blank())
# p7 <- ggplot(bytimeobs, aes(factor(model), failures)) + geom_bar(stat="identity")+
#   scale_y_continuous(trans="log", breaks=c(1, 2 ,10, 100, 500, 1000, 2000)) +
#   xlab("") + ylab("") + ggtitle("Failures observed") +
#   theme(axis.ticks.x = element_blank(),
#         axis.text.x = element_blank())
# cowplot::plot_grid(p1 + no_legend, cowplot::plot_grid(p7,p6,ncol=1),ncol=2)             
