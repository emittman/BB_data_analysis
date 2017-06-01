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
  q <- quantile(newtime, c(.25,.5,.75))
  q2 <- quantile(loss, c(.25,.5,.75))
  data.frame(model = dm, mean_time = meantime, mean_loss = mean(loss),
             lower = q[1], median = q[2], upper = q[3],
             lower2 = q2[1], median2 = q2[2], upper2 = q2[3])
})

ggplot(res, aes(x=mean_time)) + geom_histogram()

library(ggplot2)
res$model <- factor(res$model, levels=res$model[order(res$median, decreasing=TRUE)])
p1 <- res %>%
  filter(median>9) %>%
  ggplot(aes(x=model, y=median, ymin=lower, ymax=upper)) + geom_pointrange() +
  scale_y_continuous() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  ggtitle("Predicted Time to Failure") + ylab("")

p2 <- res %>%
  filter(median>9) %>%
  ggplot(aes(x=model, y=annualized_loss*10)) + geom_point(color="red")

p3 <- res %>%
  filter(median>9) %>%
  ggplot(aes(x=model, y=annualized_loss, ymin=lower2, ymax=upper2)) +
  geom_pointrange()+
  ylab("") +
  scale_y_continuous(trans="sqrt", breaks=c(.01,.05,.10,.15,.25))+
  ggtitle("Annualized adjusted loss")

c1 <- cowplot::plot_grid(p1,p3,ncol=1)

p4 <- res %>%
  ggplot(aes(x=model, y=median2, ymin=lower2, ymax=upper2)) + geom_pointrange(linetype=2) +
  scale_y_continuous()+#geom_hline(yintercept=9, color="red") +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  ylab("Predicted TTF for new unit")
cowplot::plot_grid(p4, c1)

dat$model_name[match(overview[overview$stan_id %in% c(2,3,4,5),]$model, dat$model)]