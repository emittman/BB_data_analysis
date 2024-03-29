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

id <- with(overview, which(overview$early >= 0 & overview$late_f >= 0 & f >=3))
overview$stan_id <- NA
overview[id,]$stan_id <- 1:length(id)

sfull <- readRDS("../workflow/samples_lor_only3fails.rds")
tr_adj <- readRDS("../BB_data/tr_adj_tp2s2pi.rds")$median
tr_adj2 <- readRDS("../BB_data/tr_adj_fullmodel2.rds")

sampfull <- extract(sfull)
source("../plotting_fns/KM_plot2.R")
source("../predictive/sample_from_predictive.R")
source("../plotting_fns/lifetime_plot.R")

xlabels=c(.1,.5,2.5,10,50)
ylabels=c(.001,.005, .01,.05,.1)
xlimits=c(100, 50000)
ylimits=c(.001,.2)
font_size=14

dat4 <- filter(dat, model==overview$model[which(overview$stan_id==4)])
dat4 <- arrange(dat4, endtime)
dat4_tr <- dat4[which(dat4$starttime>800),]

set.seed(80717)

lt_plot <- lifetime_plot4(data = dat4, n_to_show=60, in_years=FALSE, lab="Drive-model 4",
                          trans="log", xlimits = xlimits, xlabels=xlabels,
                          font_size = font_size)
  # theme(axis.text.x = element_blank(),
        # axis.title.x = element_blank()
        # )
  

km4 <- KM.survfit(dat4)
km4 <- tr.adj(km4, tr_adj[4])
km4_tr <- KM.survfit(dat4_tr)
#calculate adjustment: estimate of proportion failing at first starttime
first_start <- min(dat4_tr$starttime)
draws <- with(sampfull, sapply(1:24000, function(iter){
  1 - (1 - exp(log_pi[iter,4]) * 
         my_pweibull(first_start, mu1[iter], sigma1[iter])) * 
    (1 - my_pweibull(first_start, mu2[iter,4], sigma2[iter,4]))}))

summary(draws)

km4_tr <- tr.adj(km4_tr, median(draws))

bp <- baseKMplot(fit=km4, xlimits=xlimits, ylimits=ylimits, color="black",
                 linetype = 4242, alpha = 1, logscale = TRUE, label="adjusted KM\nestimate (all drives)")

samp4 <- with(sampfull, list(mu1=mu1, sigma1=sigma1, log_pi=log_pi[,4], mu2=mu2[,4], sigma2=sigma2[,4]))
band4 <- bandFromPSamp(samples=samp4, range=xlimits, length.out = 50, N=300, logscale = TRUE)
bpp <- addBandToBaseplot(baseplot=bp, bandObj=band4, color="black",
                         linetype=11, label="\nposterior median\n(90% interval)", alpha = .3)

bppp <- addKmToBaseplot(baseplot = bpp, fitObj = km4_tr, color="black",size=1,
                        linetype="solid", label="\nadjusted KM\nestimate excluding\nnewest 75 drives")

combined <- plotFinally(plotList=bppp, xbrks=xlabels*1000, ybrks=ylabels, years=FALSE) +
  guides(fill=FALSE) + ggtitle("Drive-model 4, lifetime estimates") +theme_bw(base_size=14)

library(cowplot)

# now extract the legends
legendLT <- get_legend(lt_plot + theme(legend.text = element_text(size=12)))
legendcomb <- get_legend(combined+ theme(legend.text = element_text(size=12)))

# and replot suppressing the legend
lt_plot <- lt_plot + theme(legend.position = "none")
combined <- combined + theme(legend.position='none')

# Now plots are aligned vertically with the legend to the right
ggdraw(plot_grid(plot_grid(lt_plot, combined, ncol=1, align='v'),
                 plot_grid(legendLT, legendcomb, ncol=1),
                 rel_widths=c(1.3, 0.3)))
ggsave("../paper/fig/dm4-exception.png", width=9, height=11)

# ###################
# # removed, saved for later use
# bpp1 <- addKmToBaseplot(baseplot = bpp, fitObj = km4e1, color="black", linetype="solid", label="excl. 1")
# bpp12 <- addKmToBaseplot(baseplot = bpp1, fitObj = km4e2, color="black", linetype="solid", label="excl. 2")
# bpp123 <- addKmToBaseplot(baseplot = bpp12, fitObj = km4e3, color="black", linetype="solid", label="excl. 3")
# bpp1234 <- addKmToBaseplot(baseplot = bpp123, fitObj = km4e4, color="black", linetype="solid", label="excl. 4")
# bpp12345 <- addKmToBaseplot(baseplot = bpp1234, fitObj = km4e5, color="black", linetype="solid", label="excl. 5")
# 
# 
# 
# sampGlob <- with(sampfull, list(mu1=mu1, sigma1=sigma1, log_pi=eta_pi, mu2=eta_tp2 - qsev(.5)*exp(eta_s2), sigma2=exp(eta_s2)))
# bandGlob <- bandFromPSamp(samples=sampGlob, range=xlimits, length.out=50, N=300, logscale = TRUE)
# bppp <- addBandToBaseplot(baseplot = bpp, bandObj=bandGlob, color = "black", linetype = "dotted", alpha = 0, label = "global")
# 
