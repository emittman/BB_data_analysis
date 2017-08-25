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

set.seed(82417)

make_replicate_data <- function(orig, mcmc_samp, dm, samp_id){
  n <- nrow(orig)
  id_obs <- which(orig$censored==0)
  id_cens <- which(!(1:n %in% id_obs))
  t_tr <- orig$starttime
  t_cens <- rep(max(orig$endtime), n)
  t_cens[id_cens] <- orig$endtime[id_cens]
  
  replicates <- lapply(samp_id, function(s){
  y <- with(mcmc_samp, sample_glfp_rep(n, p=exp(log_pi[s,dm]),
                                      m1=mu1[s], s1=sigma1[s],
                                      m2=mu2[s,dm], s2=sigma2[s,dm],
                                      left_tr=t_tr))
  
    is_cens <- y>t_cens
    y[is_cens] <- t_cens[is_cens]
  
    id_ok <- which(y>t_tr)
    n_new <- length(id_ok)
    print(cat(c("length ok:", n_new, "\n")))
    print(cat(c("number failed:", sum(!is_cens[id_ok]), "\n")))
    if(sum(!is_cens[id_ok]) == 0) return(NULL)
    out <-KM.survfit(data.frame(starttime=t_tr[id_ok], endtime=y[id_ok],
                          censored=is_cens[id_ok]), alpha = .1)
    print(cat(c("afterfit\n\n")))
    out
  })
  null_id <- which(sapply(replicates, is.null))
  if(length(null_id)>0){replicates <- replicates[-null_id]}
  replicates
}

reps <- 19

y_reps <- lapply(1:44, function(dm){
  dat_dm <- filter(dat, model==overview$model[which(overview$stan_id==dm)])
  samp_id <- sample(length(sampfull$mu1), reps)
  replicates <- make_replicate_data(orig = dat_dm, mcmc_samp = sampfull, dm=dm, samp_id = samp_id)
  replicates
})

saveRDS(y_reps,"../plots/y_rep.rds")

get_plotlist <- function(y_reps, dat){
  plist <- list()
  xaxis <- NULL
  yaxis <- NULL
  legend <- NULL
  n <- length(y_reps)
  for(i in 1:n){
    xlimits <- c(100,100000)
    ylimits <- c(.001, .9)
    xbrks <- 0:10 * 10000
    ybrks <- c(.001,.01,.1,.25,.5,.75,.9)
    bp <- baseKMplot.multiple(y_reps[[i]], xlimits=xlimits, ylimits=ylimits,
                              color="black",linetype = "dashed", logscale = FALSE,
                              label=expression(KM(y[rep])), alpha=.4)
    dat_dm <- filter(dat, model==overview$model[which(overview$stan_id==i)])
    km_obs <- KM.survfit(dat_dm)
    bp <- addKmToBaseplot(bp, km_obs, color="red",linetype = "solid",
                          label = expression(KM(y[obs])))
    pnew <- plotFinally(bp, xbrks=xbrks, ybrks=ybrks, years=FALSE) +
      theme_bw(base_size=14) + theme(axis.title.y=element_blank(),
                                     plot.margin = unit(c(0,0,0,0), "cm")) + 
      ggtitle(as.character(i))
  
    if(i==1) legend <- cowplot::get_legend(pnew)
    
    pnew <- pnew + theme(legend.position = "none")
    # if(!(i %in% firstcol)){
    #   pnew <- pnew + theme(axis.text.y=element_blank())
    # } 
    # if(!(i %in% lastrow)){
    #   pnew <- pnew + theme(axis.title.x=element_blank(),
    #                        axis.text.x=element_blank(),
    #                        plot.margin = unit(c(0,0,-1,0), "cm"))
    # }
    plist[[i]] <-pnew
  }
  list(plots=plist, xaxis=xaxis, yaxis=yaxis, legend=legend)
}
set.seed(10438)
# ran_set <- sort(sample(44, 6))
all <- 1:44

plist <- get_plotlist(y_reps, dat)

plist <- get_plotlist(all, firstcol = 1+(0:8)*5, lastrow = 41:44)
library(cowplot)
pdf("../paper/fig/post-pred-KM-all2.pdf", height=14, width=9)
plot_grid(plotlist = plist$plots, ncol=5, align="hv")
dev.off()

legnd <- get_legend(plist$plots[[1]] + theme(legend.position="right"))

no.x.title <- lapply(plist$plots, function(p) p + theme(axis.title = element_blank()))
for(j in 1:6){
  if(j<6){
    plot_grid(plot_grid(plotlist=plist$plots[1:8 + (j-1)*8], ncol=2, align="h"), plist$legend,
              nrow=1, rel_widths = c(1,.2))
  } else{
    plot_grid(plot_grid(plotlist=plist$plots[1:4 + (j-1)*8], ncol=2, align="h"), plist$legend,
                     nrow=1, rel_widths = c(1,.2))
  }
  ggsave(paste("../paper/fig/ppcheck",j,".pdf"), width = 8.5, height=11)
}

set.seed(12122017)
smplID <- sort(sample(44, 4))
plot_grid(plot_grid(plotlist=plist$plots[c(smplID)], ncol=2, align="hv"), plist$legend,
          nrow=1, rel_widths = c(1,.2))

ggsave("../paper/fig/ppcheck-sample.pdf", width=11, height=8)
#debug
# p <- get_plotlist(42)
# 
# # length(plist)
# pdf("../paper/fig/pospredictKM.pdf")
# plot_grid(plist[[1]], plist[[2]],plist[[3]],plist[[4]],plist[[5]],plist[[6]], ncol=3)
# dev.off()
# 
# plist2 <- get_plotlist(7:12)
# plot_grid(plist2[[7]], plist2[[8]],plist2[[9]],plist2[[10]],plist2[[11]],plist2[[12]], ncol=3)
# 
# plist3 <- get_plotlist(13:18)
# plot_grid(plist3[[13]], plist3[[14]], plist3[[15]], plist3[[16]], plist3[[17]], plist3[[18]], ncol=3)
# 
# plist4 <- get_plotlist(19:24)
# plot_grid(plist4[[19]], plist4[[20]], plist4[[21]], plist4[[22]], plist4[[23]], plist4[[24]], ncol=3)
# 
# plist5 <- get_plotlist(25:30)
# plot_grid(plist5[[25]], plist5[[26]], plist5[[27]], plist5[[28]], plist5[[29]], plist5[[30]], ncol=3)
# 
# plist6 <- get_plotlist(c(31:33,35:36))
# plot_grid(plist6[[31]], plist6[[32]], plist6[[33]], plist6[[35]], plist6[[36]], ncol=3)
# 
# plist7 <- get_plotlist(c(37:39))
# plist8 <- get_plotlist(c(40:41))
# plist9 <- get_plotlist(c(44))
# plot_grid(plist7[[37]],plist7[[38]],plist7[[39]],plist8[[40]],plist8[[41]],plist9[[44]], ncol=3)
