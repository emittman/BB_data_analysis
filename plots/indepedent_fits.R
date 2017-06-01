setwd("workflow2/")
#3 8  14
s <- readRDS("samples_nopool_data_dm14.rds")
print(s)
library(plyr)
library(rstan)
library(dplyr)
library(ggplot2)
dat <- readRDS("../BB_data/clean_unit_summaries.rds")
ov <- readRDS("../BB_data/overview.rds")
mod_id <- filter(ov, f>=5, early_f>=1, late_f>=1) %>%
  select(model) %>%
  unlist
xlim = c(2000,50000)
ylim = c(.005,.8)
todo <- c(3,8,14)

#these corresponds to full_model ids: 
# match(mod_id[todo], filter(ov, f>=3)$model) 
#[1]  6 14 23

plots <- list()
bands <- list()

for(id in 1:3){
  dm = todo[id]
  subdat <- filter(dat, model == mod_id[dm])
  samp <- extract(readRDS(paste(c("samples_nopool_data_dm",dm,".rds"), collapse="")))
  plots[[id]] <- KM_plot_NP(subdat, "weibull", conf=.9, xlimits=xlim, ylimits=ylim)+
    theme_bw()+
    theme(legend.position="none",
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.grid = element_blank()
          )
  bands[[id]] <- KM_band("1", 1, samp, 100, xlim=xlim, ylim=ylim, pi.free=F, mu1.free=F, sigma1.free=F,
               mu2.free=F, sigma2.free=F)
}


pdf("../plots/ind_estimates.pdf")
library(cowplot)
plot_grid(plots[[1]]+bands[[1]],
          plots[[2]]+bands[[2]],
          plots[[3]]+bands[[3]], ncol=3)
dev.off()
