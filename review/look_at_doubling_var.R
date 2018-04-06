s <- readRDS("workflow/samples_lor_only3fails.rds")
sdbl <- readRDS("workflow/samples_doublevar.rds")

library(rstan)

`original posterior`<- extract(s)
`inflated variance posterior` <- extract(sdbl)

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

normal.priors <- data.frame(param = rep(c("eta_tp2","eta_s2","eta_pi"),each=2),
                            Distribution = rep(c("original prior","inflated variance prior"), 3),
                            mu = rep(c(9,0,-3), each=2),
                            sigma = c(2,4,2,4,1,2))

df_plot_prior <- ddply(normal.priors, .(param,Distribution), function(x){
  data.frame(param=x$param, Distribution = x$Distribution, value=rnorm(10000, x$mu,x$sigma))
})

df_plot_poster <- ldply(c("original posterior","inflated variance posterior"), function(m){
  ldply(c("eta_tp2","eta_s2","eta_pi"), function(p){
    data.frame(value=get(m)[[p]], Distribution=m, param=p)
  })
})
df_plot = rbind(df_plot_prior,df_plot_poster)
df_plot$Distribution <- factor(df_plot$Distribution, levels=c("original prior","inflated variance prior","original posterior","inflated variance posterior"))
ggplot(df_plot, aes(value)) + geom_density(aes(fill=Distribution, color=Distribution), alpha=.5, lty=2) + facet_wrap(~param, scale="free",ncol=1)+
  theme_bw(base_size=14) + xlab("")
