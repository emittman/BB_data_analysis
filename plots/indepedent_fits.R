setwd("workflow2/")
#3 8  14
s <- readRDS("samples_nopool_data_dm8.rds")
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
xlim = c(1000,55000)
ylim = c(.008,.8)
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
          axis.title.y = element_blank(),
          panel.grid = element_blank()
          )+
    scale_x_continuous(limits=xlim, breaks=c(1000, 2000, 5000, 10000, 20000, 40000), trans="log")+
    scale_y_continuous(limits=ylim, breaks=c(.01, .05, .1, .2, .5, .7, .9, .95, .99), trans="qsev")
  bands[[id]] <- KM_band("1", 1, samp, 2000, xlim=xlim, ylim=ylim, pi.free=F, mu1.free=F, sigma1.free=F,
               mu2.free=F, sigma2.free=F)
}

no_yaxis_nums <- theme(axis.text.y=element_blank())
library(cowplot)

pdf("../plots/ind_estimates.pdf", width=10, height=5)
plot_grid(plots[[1]]+bands[[1]]+ggtitle("6"),
          plots[[2]]+bands[[2]]+no_yaxis_nums+ggtitle("14"),
          plots[[3]]+bands[[3]]+no_yaxis_nums+ggtitle("23"), ncol=3)
dev.off()

plots[[1]] + geom_ribbon(data=bands[[1]]$band, inherit.aes = F,
                         aes(x=x, y=y, ymin=lower, ymax=upper), alpha=.2)
