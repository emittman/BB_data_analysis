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
  data.frame(param=x$param, Distribution = x$Distribution, value=rnorm(100000, x$mu,x$sigma))
})

df_plot_poster <- ldply(c("original posterior","inflated variance posterior"), function(m){
  ldply(c("eta_tp2","eta_s2","eta_pi"), function(p){
    data.frame(value=get(m)[[p]], Distribution=m, param=p)
  })
})
df_plot = rbind(df_plot_prior,df_plot_poster)
df_plot$Distribution <- factor(df_plot$Distribution, levels=c("original prior","inflated variance prior","original posterior","inflated variance posterior"))


p1 <- filter(df_plot, param == "eta_pi") %>%
ggplot(aes(value)) + geom_density(aes(fill=Distribution, color=Distribution), alpha=.5, lty=2) + 
  xlim(c(-6,1))+
  # facet_wrap(~param, scale="free",ncol=1)+
  theme_bw(base_size=14) + xlab("") + ylab("") + ggtitle(expression(eta[pi])) +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5)) #+
#  guides(color="none", fill="none")

p2 <- filter(df_plot, param == "eta_tp2") %>%
  ggplot(aes(value)) + geom_density(aes(fill=Distribution, color=Distribution), alpha=.5, lty=2) + 
  xlim(c(6,13))+
  # facet_wrap(~param, scale="free",ncol=1)+
  theme_bw(base_size=14) + xlab("") + ylab("") + ggtitle(expression(eta[t[p[2]]])) +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  guides(color="none", fill="none")

p3 <- filter(df_plot, param == "eta_s2") %>%
  ggplot(aes(value)) + geom_density(aes(fill=Distribution, color=Distribution), alpha=.5, lty=2) + 
  xlim(c(-4,4))+
  # facet_wrap(~param, scale="free",ncol=1)+
  theme_bw(base_size=14) + xlab("") + ylab("") + ggtitle(expression(eta[sigma[2]])) +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  guides(color="none", fill="none")

grobs <- ggplotGrob(p1)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
p1 <- p1 + guides(color="none", fill="none")
library(cowplot)

pdf(file="paper/fig/double_var.pdf", width=11, height=6)
plot_grid(p1,p2,p3,legend, ncol=2)
dev.off()
