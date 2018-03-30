s <- readRDS("workflow/samples_lor_only3fails.rds")
sdbl <- readRDS("workflow/samples_doublevar.rds")

library(rstan)

original <- extract(s)
`inflated variance` <- extract(sdbl)

plot(s, pars=c("eta_tp2","eta_s2","eta_pi"))
plot(sdbl, pars=c("eta_tp2","eta_s2","eta_pi"))

library(plyr)

plot_dat <- ldply(c("original","inflated variance"), function(mod){
  ldply(c("eta_tp2","eta_s2","eta_pi"), function(par){
    data.frame(model = mod, param = par, t(quantile(get(mod)[[par]], c(.05,.5,.95))))
  })
})

names(plot_dat)[3:5] <- c("lower","median","upper")


library(ggplot2)

ggplot(plot_dat, aes(x=model, ymin=lower, ymax=upper, y=median))+
  geom_pointrange(aes(color=param)) + facet_wrap(~param, scales="free") +
  theme_bw()
