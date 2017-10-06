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

#Figure 14 (index 8 in independent fits)
indFit <- extract(readRDS("../workflow2/samples_nopool_data_dm8.rds"))
samp14 <- with(indFit, list(mu1=mu1, sigma1=sigma1, mu2=mu2, sigma2=sigma2, log_pi=log_pi))

#Get Summary Quantiles for Paper
quantile(1/samp14$sigma1, probs=c(.025, .5, .975))
q1 <- exp(samp14$mu1 + samp14$sigma1*log(-log(1-.5)))
quantile(q1, probs=c(.025, .5, .975))
q2 <- exp(samp14$mu2 + samp14$sigma2*log(-log(1-.2)))
quantile(q2, probs=c(.025, .5, .975))

#Linear scales
xlimits <- c(200, 30000)
ylimits <- c(.001, .8)
xbreaks <- 1:6*5000
ybreaks <- 1:8*0.1

km14 <- KM.survfit(filter(dat, model==overview$model[which(overview$stan_id==14)]),
                   greenwood = TRUE, alpha=.1)
base_p14 <- baseKMplot(km14, xlimits=xlimits, ylimits=ylimits, color="blue", linetype="solid",
                       logscale = FALSE, prob=FALSE, label = "Kaplan-Meier")
p14 <- plotFinally(base_p14, xbrks = xbreaks, ybrks = ybreaks, years=FALSE,
                   greenwood = TRUE, alpha = .1)
p14 + theme_bw(base_size=14) + theme(legend.position="none")
ggsave("../paper/fig/km14-natural.pdf", width = 6, height = 6)

#Add independent fit
band14_nat <- bandFromPSamp(samp14, range=xlimits, length.out = 40, N=1000, logscale = FALSE)
comb_p14 <- addBandToBaseplot(baseplot = base_p14, bandObj = band14_nat,
                              color="red", linetype="dashed", alpha=.1, label="GLFP posterior")
p14plus <- plotFinally(comb_p14, xbrks = xbreaks, ybrks = ybreaks, years=FALSE,
                       greenwood = TRUE, alpha = .1)
p14plus + theme_bw(base_size=14)
ggsave("../paper/fig/km14-natural-plus.pdf", width = 7.5, height=6)

#Probability + log scale
ybreaks<- c(.002,.01,.05,.1,.2,.4,.55)
xbreaks<- c(500, 1000, 2000, 4000, 8000, 15000, 30000)
base_p14_p <- baseKMplot(km14, xlimits=xlimits, ylimits=ylimits, color="blue", linetype="solid",
                         logscale = TRUE, prob=TRUE)
p14_p <- plotFinally(base_p14_p, xbrks = xbreaks, ybrks = ybreaks, years=FALSE,
                   greenwood = TRUE, alpha = .1)
p14_p + theme_bw(base_size=14) + theme(legend.position="none")
ggsave("../paper/fig/km14-prob.pdf", width = 6, height = 6)

band14_nat <- bandFromPSamp(samp14, range=xlimits, length.out = 40, N=1000, logscale = TRUE)
comb_p14_p <- addBandToBaseplot(baseplot = base_p14_p, bandObj = band14_nat,
                          color="red", linetype="dashed", alpha=.1, label="GLFP posterior")
p14_pplus <- plotFinally(comb_p14_p, xbrks = xbreaks, ybrks = ybreaks, years=FALSE,
                         greenwood = TRUE, alpha = .1)
p14_pplus + theme_bw(base_size=14)
ggsave("../paper/fig/km14-prob-plus.pdf", width=7.5, height=6)