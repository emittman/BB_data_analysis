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
  source("../plotting_fns/lifetime_plot.R")
  
  xlabels=c(1,2,5,10, 20)
  ylabels=c(.001,.01,.1,.5,.9, .99, .999)
  xlimits=c(5000, 170000)
  ylimits=c(.001,.999)
  
  datdm <- filter(dat, model==overview$model[which(overview$stan_id==dm)])
  lt_plot <- lifetime_plot3(datdm, xlabels = xlabels, in_years=TRUE,
                            lab=paste(c("Drive-model ",dm),collapse=""), trans="log", xlimits = xlimits)
  
  kmdm <- KM.survfit(datdm, greenwood = FALSE)
  kmdm <- tr.adj(kmdm, tr_adj[dm])
  
  bp <- baseKMplot(fit=kmdm, xlimits=xlimits, ylimits=ylimits, color="black",
                   linetype = "solid", alpha = 1, logscale = TRUE, label="nonparametric", prob = .9)
  
  sampdm <- with(sampfull, list(mu1=mu1, sigma1=sigma1, log_pi=log_pi[,dm], mu2=mu2[,dm], sigma2=sigma2[,dm]))
  banddm <- bandFromPSamp(samples=sampdm, range=xlimits, length.out = 70, N=200, logscale = TRUE)
  bpp <- addBandToBaseplot(baseplot=bp, bandObj=banddm, color="black",
                           linetype="dashed", label="posterior median\n(90% interval)", alpha = .3)
  
  sampGlob <- with(sampfull, list(mu1=mu1, sigma1=sigma1, log_pi=eta_pi, mu2=eta_tp2 - qsev(.5)*exp(eta_s2), sigma2=exp(eta_s2)))
  bandGlob <- bandFromPSamp(samples=sampGlob, range=xlimits, length.out=50, N=300, logscale = TRUE)
  bandGlob$band <- bandGlob$band[which(bandGlob$band$est<ylimits[2]),]
  bppp <- addBandToBaseplot(baseplot = bpp, bandObj=bandGlob, color = "black", linetype = "dotted", alpha = 0, label = "global")
  
  combined <- plotFinally(plotList=bppp, xbrks=xlabels*24*365, ybrks=ylabels, years=TRUE, greenwood = FALSE) +
    guides(fill=FALSE)
  
  library(cowplot)
  # now extract the legends
  legendLT <- get_legend(lt_plot + theme(legend.text = element_text(size=7)))
  legendcomb <- get_legend(combined+ theme(legend.text = element_text(size=7)))
  
  # and replot suppressing the legend
  lt_plot <- lt_plot + theme(legend.position = "none")
  combined <- combined + theme(legend.position='none')
  
  # Now plots are aligned vertically with the legend to the right
  pdf(paste(c("../paper/fig/dm",dm,"-shrinkage.pdf"), collapse=""), width=6, height=6)
  ggdraw(plot_grid(plot_grid(lt_plot, combined, ncol=1, align='v'),
                   plot_grid(legendLT, legendcomb, ncol=1),
                   rel_widths=c(1, 0.2)))
  dev.off()

