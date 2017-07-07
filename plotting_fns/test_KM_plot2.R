#use new plotting functions
library(rstan)
library(plyr)
library(dplyr)

source("plotting_fns/KM_plot2.R")

s <- readRDS("workflow2/samples_nopool_data_dm8.rds") #number 8 after applying filter (f>=5, ef>=1, lf>=1)
samp <- extract(s)

dat14 <- readRDS("BB_data/clean_unit_summaries.rds") %>%
  filter(model == 20) #original model number; what we call 14

km14 <- KM.survfit(dat14, greenwood = TRUE)

xlimits <- range(km14$t)
ylimits <- range(km14$Ft)

p14base <- baseKMplot(km14, xlimits = xlimits, ylimits = ylimits,
           color="blue", linetype = 1)


band14 <- bandFromPSamp(samp, xlimits, length.out = 50, N = 1000, logscale = TRUE)

p14base <- addBandToBaseplot(p14base, band14, color="red", linetype=3, alpha=.2, label="Bayes 90%")

plotFinally(p14base, xbrks = c(.1, .5, 1, 2, 4)*365*24,
            ybrks = c(.01, .02, .05, .1, .2, .5), years=TRUE) +
  theme_bw()

plotFinally(p14base, xbrks = c(.1, .5, 1, 2, 4)*365*24,
            ybrks = c(.01, .02, .05, .1, .2, .5), years=TRUE, greenwood = TRUE) +
  theme_bw()

##################################################
sfull <- readRDS("workflow/samples_lor_only3fails.rds")
sampfull <- extract(sfull)
sampfull14 <- with(sampfull, data.frame(mu1 = mu1, sigma1 = sigma1,
                                        mu2 = mu2[,14], sigma2 = sigma2[,14],
                                        log_pi = log_pi[,14]))
band14full <- bandFromPSamp(sampfull14, xlimits, length.out = 50, N = 1000, logscale = TRUE)
p14base <- addBandToBaseplot(p14base, band14full, color="yellow", linetype=2, alpha=.2, label="Bayes hier. 90%")
plotFinally(p14base, xbrks = c(.1, .5, 1, 2, 4)*365*24,
            ybrks = c(.01, .02, .05, .1, .2, .5), years=TRUE) +
  guides(fill = FALSE) +
  theme_bw()

