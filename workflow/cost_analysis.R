library(rstan)
# setwd("workflow/")
samp <- extract(readRDS("samples_lor_only3fails.rds"))
source("../predictive/sample_from_predictive.R")
# tr_adj <- readRDS("../BB_data/tr_adj_tp2s2pi.rds")
V <- function(x){
  #devaluation curve
  exp(-.4*x)
}

library(plyr)
library(dplyr)

res <- ldply(1:44, function(dm){
  p <- exp(samp$log_pi[,dm])
  m1 <- samp$mu1
  s1 <- samp$sigma1
  m2 <- samp$mu2[,dm]
  s2 <- samp$sigma2[,dm]
  newtime <- sample_glfp(n=length(samp$lp__), p, m1, s1, m2, s2)/24/365
  meantime <- mean(newtime)
  loss <- V(newtime)
  q <- quantile(newtime, c(.05,.5,.95))
  q2 <- quantile(loss, c(.05,.5,.95))
  data.frame(model = dm, mean_time = meantime, mean_loss = mean(loss),
             lower = q[1], median = q[2], upper = q[3],
             lower2 = q2[1], median2 = q2[2], upper2 = q2[3])
})

ggplot(res, aes(x=mean_time)) + geom_histogram()

library(ggplot2)
res$model <- factor(res$model, levels=res$model[order(res$mean_loss, decreasing=FALSE)])
p1 <- res %>%
  filter(mean_loss<1/8) %>%
  ggplot(aes(x=model, y=mean_time, ymin=lower, ymax=upper)) + geom_pointrange() +
  # geom_point(aes(y=mean_time), shape=2, color=2)+
  scale_y_continuous(trans="log10", breaks=c(1,2,5,10,20,50),limits=c(1,65)) +
  theme_bw()+
  theme(axis.title.x = element_blank())+
        # axis.ticks.x = element_blank(),
        # axis.text.x = element_blank()) +
  ggtitle("Predictive Time to Failure") + ylab("")

p2 <- res %>%
  filter(mean_loss<1/8) %>%
  ggplot(aes(x=model, y=mean_loss)) + geom_point(color="red")

p3 <- res %>%
  filter(mean_loss<1/8) %>%
  ggplot(aes(x=model, y=mean_loss, ymin=lower2, ymax=upper2)) +
  geom_pointrange()+
  geom_hline(yintercept=.11, color="red") +
  # geom_point(aes(y=mean_loss), shape=2, color=2)+
  ylab("") + theme_bw()+ xlab("top drive-models") +
  scale_y_continuous(trans="sqrt", breaks=c(.01,.05,.10,.2, .5))+
  ggtitle("Predictive VAF")

c1 <- cowplot::plot_grid(p1,p3,ncol=1)

p4 <- res %>%
  ggplot(aes(x=model, y=mean_loss, ymin=lower2, ymax=upper2)) + geom_pointrange(alpha=.5)+
  # geom_pointrange(aes(y=median2))+
  theme_bw()+
  scale_y_continuous()+ geom_hline(yintercept=.11, color="red") +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_y_continuous(trans="sqrt", breaks=c(.01,.05,.10,.2, .5,1))+
  
  ylab("Predicted Value at Failure")
cowplot::plot_grid(p4, c1)

dat$model_name[match(overview[overview$stan_id %in% c(2,3,4,5),]$model, dat$model)]
