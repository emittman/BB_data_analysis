library(plyr)
library(dplyr)
library(rstan)
library(ggplot2)

setwd("plots/")
dat <- readRDS("../BB_data/clean_unit_summaries.rds")
dat$model <- as.integer(dat$model)

overview <- readRDS("../BB_data/overview.rds")
id <- with(overview, which(overview$early >= 0 & overview$late_f >= 0 & f >=3))
overview$stan_id <- NA
overview[id,]$stan_id <- 1:length(id)

source("../plotting_fns/KM_plot2.R")

ssigtp <- readRDS("../workflow_sig2_and_tp2_vary/vary_s2_and_tp2_4_17.rds")
stp2 <- readRDS("../workflow_tp2_vary/samples_tp2_vary.rds")
snull <- readRDS("../workflow_null/samples_null_model_3_29.rds")
sfull <- readRDS("../workflow/samples_lor_only3fails.rds")
tr_adj <- readRDS("../BB_data/tr_adj_tp2s2pi.rds")$median

samps2tp <- extract(ssigtp)
samptp2 <- extract(stp2)
sampnull <- extract(snull)
sampfull <- extract(sfull)

xlimits <- c(200, 30000)
ylimits <- c(.001, .8)
xbreaks <- 1:6*5000
ybreaks <- 1:8*0.1

#local function to produce ggplot with 4 curves for each drive-model
single_drive.all.4.Models <- function(dm, prob.log.scales){
  local_dat <- filter(dat, model==overview$model[which(overview$stan_id==dm)])
  
  km <- KM.survfit(local_dat, greenwood = FALSE)
  km <- tr.adj(km, tr_adj = tr_adj[dm])  
  if(!prob.log.scales){
    xlimits <- range(km$t) + c(-1000,1000)
    xlimits[1] <- max(xlimits[1], 0)
    ylimits <- range(km$Ft) + c(-.5, .2)
    ylimits[1] <- max(ylimits[1], 0)
    ylimits[2] <- min(ylimits[2], 1)
    xbreaks <- 0:ceiling(xlimits[2]/1000)*1000
    ybreaks <- 0:ceiling(ylimits[2]*20) /20
  } else{
    xlimits <- exp(range(log(km$t))+c(-.1,.1))
    ylimits <- psev(qsev(range(km$Ft)) + c(-1, 1))
    ylimits[1] <- min(ylimits,.001)
    xbreaks <- c(100, 500, 1000, 2000, 5000,10000, 20000, 40000, 80000, 160000)
    ybreaks <- signif(psev(seq.int(from=floor(qsev(ylimits[1])),
                                   ceiling(qsev(ylimits[2])), by=.5)),1)
  }
  if(!prob.log.scales){
    base_plot <- baseKMplot(km, xlimits=xlimits, ylimits=ylimits, color="black", linetype="solid",
                           logscale = FALSE, prob=FALSE, label = "adjusted Kaplan-Meier")
  } else{
    base_plot <- baseKMplot(km, xlimits=xlimits, ylimits=ylimits, color="black", linetype="solid",
                            logscale = TRUE, prob=TRUE, label = "adjusted Kaplan-Meier")
  }
  #Model 1
  samp1 <- with(sampnull, list(mu1=mu1, mu2=mu2, sigma1=sigma1, sigma2=sigma2, log_pi=log_pi))
  if(!prob.log.scales){
    band1 <- bandFromPSamp(samp1, range=xlimits, length.out = 30, N=1000, logscale = FALSE)
  } else{
    band1 <- bandFromPSamp(samp1, range=xlimits, length.out = 30, N=1000, logscale = TRUE)
  }
  
  #Model 2
  samp2 <- with(samptp2, list(mu1=mu1, mu2=mu2[,dm], sigma1=sigma1, sigma2=sigma2, log_pi=log_pi))
  if(!prob.log.scales){
    band2 <- bandFromPSamp(samp2, range=xlimits, length.out = 30, N=1000, logscale = FALSE)
  } else{
    band2 <- bandFromPSamp(samp2, range=xlimits, length.out = 30, N=1000, logscale = TRUE)
  }
  
  #Model 3
  samp3 <- with(samps2tp, list(mu1=mu1, mu2=mu2[,dm], sigma1=sigma1, sigma2=sigma2[,dm], log_pi=log_pi))
  if(!prob.log.scales){
    band3 <- bandFromPSamp(samp3, range=xlimits, length.out = 30, N=1000, logscale = FALSE)
  } else{
    band3 <- bandFromPSamp(samp3, range=xlimits, length.out = 30, N=1000, logscale = TRUE)
  }
  
  #Model 4
  samp4 <- with(sampfull, list(mu1=mu1, mu2=mu2[,dm], sigma1=sigma1, sigma2=sigma2[,dm], log_pi=log_pi[,dm]))
  if(!prob.log.scales){
    band4 <- bandFromPSamp(samp4, range=xlimits, length.out = 30, N=1000, logscale = FALSE)
  } else{
    band4 <- bandFromPSamp(samp4, range=xlimits, length.out = 30, N=1000, logscale = TRUE)
  }
  
  km_1 <- addBandToBaseplot(baseplot = base_plot, bandObj = band1,
                                color="red", linetype="3232", alpha=0, label="Model 1")
  km1_2 <- addBandToBaseplot(baseplot = km_1, bandObj = band2,
                            color="blue", linetype="6262", alpha=0, label="Model 2")
  km12_3 <- addBandToBaseplot(baseplot = km1_2, bandObj = band3,
                            color="darkgreen", linetype="b2b2", alpha=0, label="Model 3")
  km1234 <- addBandToBaseplot(baseplot = km12_3, bandObj = band4,
                            color="darkorange", linetype="3232b2b2", alpha=0, label="Model 4")
  
  p <- plotFinally(km1234, xbrks = xbreaks, ybrks = ybreaks, years=FALSE,
                         greenwood = FALSE, alpha = 0)
  p + theme_bw(base_size=14) + guides(fill=FALSE)
}

P <- single_drive.all.4.Models(10, FALSE) 
P2 <- single_drive.all.4.Models(3, TRUE)

plots.prob <- list()
for(dm in 1:44){
  plots.prob[[dm]] <- single_drive.all.4.Models(dm, TRUE)
}
save(plots.prob, file="plots-prob.RData")
library(cowplot)
lgnd <- get_legend(plots.prob[[1]])
lgnd_off <- theme(legend.position="none")
to.save=plot_grid(plot_grid(plots.prob[[2]] + ggtitle("2") + lgnd_off,
                    plots.prob[[9]] + ggtitle("9") + lgnd_off,
                    plots.prob[[14]] + ggtitle("14") + lgnd_off,
                    plots.prob[[40]] + ggtitle("40") + lgnd_off,ncol=2),
          lgnd, rel_widths = c(10,2), nrow=1)
ggsave("../paper/fig/single-drive-4-Models-ex.pdf", width=12, height=10)
# plots.prob <- lapply(1:44, function(dm) single_drive.all.4.Models(dm, FALSE))
# save(plots.nat, "ggplots-nat.RData")
