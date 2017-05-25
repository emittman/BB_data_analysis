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

source("../plotting_fns/KM_plot.R")

xlimits <- c(5, 50000)
ylimits <- c(.00001, .9)

x.breaks <- xbrks(xlimits[1], xlimits[2], N=10, prec=-2)
y.breaks <- ybrks(ylimits[1], ylimits[2], model="weibull")

x.scale <- scale_x_continuous(name="Hours", breaks=x.breaks, limits=xlimits)
y.scale <- scale_y_continuous(name="Fraction failing", breaks=y.breaks, limits=ylimits)

no_legend <- theme(legend.position="none")
x_log <- scale_x_continuous(trans="log")

p1 <- KM_plot_multi(df1, tr_adj = tr_adj, linear_axes = F,
                    xlimits = xlimits, ylimits=ylimits, grayscale = TRUE, size=1)

null_s <- readRDS("../workflow_null/samples_null_model_3_29.rds")
null_samp <- extract(null_s)
p2 <- ggplot(data.frame(x=0,y=0), aes(x,y)) +
  KM_band(1, null_samp, n_iter=1000, xlim=xlimits, ylim=ylimits,
          pi.free=F, mu2.free=F, sigma2.free=F, x_logscale = F)[[2]]+
  theme_bw(base_size = 14)+xlab("Hours")+ylab("Fraction Failing")+
  xlim(xlimits)+ylim(ylimits)


cowplot::plot_grid(p1 + no_legend,
                   p2, ncol=1)

tp_s <- readRDS("../workflow_tp2_vary/samples_tp2_vary_new.rds")
tp_samp <- extract(tp_s)
tp_layers <- lapply(1:length(id), function(j){
  KM_band(j, tp_samp, n_iter=1000, xlim=xlimits, ylim=ylimits,
          pi.free=F, mu2.free=T, sigma2.free=F, x_logscale = F)[[2]]
})
p3 <- ggplot(data.frame(x=0,y=0), aes(x,y))
for(j in 1:length(id)){
  p3 <- p3 + tp_layers[[j]]
}
p3 <- p3 + theme_bw(base_size = 14)+xlab("Hours")+ylab("Fraction Failing")+
  xlim(xlimits)+ylim(ylimits)

cowplot::plot_grid(p1 + no_legend,
                   p2,
                   p3, ncol=1)
  
tpsig_s <- readRDS("../workflow_sig2_and_tp2_vary/vary_s2_and_tp2_4_17.rds")
tpsig_samp <- extract(tpsig_s)
tpsig_layers <- lapply(1:length(id), function(j){
  KM_band(j, tpsig_samp, n_iter=1000, xlim=xlimits, ylim=ylimits,
          pi.free=F, mu2.free=T, sigma2.free=T, x_logscale = F)[[2]]
})
p4 <- ggplot(data.frame(x=0,y=0), aes(x,y))
for(j in 1:length(id)){
  p4 <- p4 + tpsig_layers[[j]]
}
p4 <- p4 + theme_bw(base_size = 14)+xlab("Hours")+ylab("Fraction Failing")+
  xlim(xlimits)+ylim(ylimits)

cowplot::plot_grid(p1 + no_legend,
                   p2,
                   p3,
                   p4, ncol=2)

tpsigpi <- readRDS("../workflow/samples_lor_only3fails.rds")
tpsigpi_samp <- extract(tpsigpi)
tpsigpi_layers <- lapply(1:length(id), function(j){
  KM_band(j, tpsigpi_samp, n_iter=1000, xlim=xlimits, ylim=ylimits,
          pi.free=T, mu2.free=T, sigma2.free=T, x_logscale = F)[[2]]
})

p5 <- ggplot(data.frame(x=0,y=0), aes(x,y))
for(j in 1:length(id)){
  p5 <- p5 + tpsigpi_layers[[j]]
}
p5 <- p5 + theme_bw(base_size = 14)+xlab("Hours")+ylab("Fraction Failing")+
  xlim(xlimits)+ylim(ylimits)

#decided to reorganize this (5-24-2017)
# cowplot::plot_grid(p1 + no_legend, cowplot::plot_grid(p2,p3,p4,p5, ncol=2), ncol=2)

bytimeobs <- ddply(df1, .(model), function(x){
  total_time = sum(x$endtime - x$starttime)
  failures = sum(x$failed>0)
  data.frame(total_time=total_time, failures=failures)
})

p6 <- ggplot(bytimeobs, aes(factor(model), total_time)) + geom_bar(stat="identity") +
  scale_y_continuous(trans="log", breaks=c(100,10000,1000000,1000000))+
  ggtitle("Total operational hours")+ylab("")+xlab("")+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
p7 <- ggplot(bytimeobs, aes(factor(model), failures)) + geom_bar(stat="identity")+
  scale_y_continuous(trans="log", breaks=c(1, 2 ,10, 100, 500, 1000, 2000)) +
  xlab("") + ylab("") + ggtitle("Failures observed") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
cowplot::plot_grid(p1 + no_legend, cowplot::plot_grid(p7,p6,ncol=1),ncol=2)             
