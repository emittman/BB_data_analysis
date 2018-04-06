library(ggplot2)

## Priors for Model 1  ##

#t, .5, 1 ~ Normal(8, 4)
t1 <- as.data.frame(exp(rnorm(100000, 8, sd=4)))
colnames(t1) <- c("y")
t1den <- ggplot(t1, aes(x=y)) + 
  geom_density() + xlim(100, 100000) + theme_bw() + 
  xlab(bquote(''~t[0.5]~''))
t1den

#t, .2, 2 ~ Normal(10, 4)
t2 <- as.data.frame(exp(rnorm(100000, 10, sd=4)))
colnames(t2) <- c("y")
t2den <- ggplot(t2, aes(x=y)) + 
  geom_density() + xlim(100, 100000) + theme_bw() + 
  xlab(bquote(''~t[0.2]~'')) + theme(axis.title.y=element_blank())
t2den

#InvLogitPi ~ Normal(-3, 2)
pi <- as.data.frame(rnorm(100000, -3, sd=2))
pi.tran <- exp(pi)/(1+exp(pi))
colnames(pi.tran) <- c("y")
piden <- ggplot(pi.tran, aes(x=y)) + 
  geom_density() + xlim(0, 1) + theme_bw() + 
  xlab(expression(pi))
piden

#Sigma Terms (Log Normal (0, 2.5))
sig <- as.data.frame(rlnorm(100000, 0, sd=2.5))
colnames(sig) <- c("y")
sigden <- ggplot(sig, aes(x=y)) + 
  geom_density() + xlim(0, 100) + theme_bw() + 
  xlab(expression(sigma[1] ~ "," ~ sigma[2])) +  theme(axis.title.y=element_blank())
sigden


#Note: x-axis is truncated so the plots don't just look like a spike; also could plot normal priors for log distn.
cowplot::plot_grid(t1den, t2den, piden, sigden, align = "v")


## Priors for Model 2 ##
#t, .5, 1 ~ Normal(8, 4)
t1m2 <- as.data.frame(exp(rnorm(100000, 7, sd=2)))
colnames(t1m2) <- c("y")
t1denm2 <- ggplot(t1m2, aes(x=y)) + 
  geom_density() + xlim(100, 100000) + theme_bw() + 
  xlab(bquote(''~t[0.5]~''))
t1denm2

#Sigma1 Term (Log Normal (0, 1))
sig1 <- as.data.frame(rlnorm(100000, 0, sd=1))
colnames(sig1) <- c("y")
sigdenm2 <- ggplot(sig1, aes(x=y)) + 
  geom_density() + xlim(0, 100) + theme_bw() + 
  xlab(expression(sigma[1])) +  theme(axis.title.y=element_blank())
sigdenm2

cowplot::plot_grid(t1denm2, sigdenm2, align = "v")

## Priors for Model 3 & 4 ##

#InvLogitPi ~ Normal(-3, 1)
etapi <- as.data.frame(rnorm(100000, -3, sd = 1))
etapi.tran <- exp(etapi)/(1 + exp(etapi))
colnames(etapi.tran) <- c("y")
etapiden <- ggplot(etapi.tran, aes(x=y)) + 
  geom_density() + xlim(0, 1) + theme_bw() + 
  xlab(expression(eta[pi]))
etapiden

#eta_sigma2 ~ Normal(0, 2)
etasig2 <- as.data.frame(exp(rnorm(10000000, 0, sd = 2)))
colnames(etasig2) <- c("y")
etasig2den <- ggplot(etasig2, aes(x=y)) + 
  geom_density() + xlim(100, 5000) + theme_bw() + 
  xlab(expression(eta[sigma[2]])) + theme(axis.title.y=element_blank())
etasig2den

#eta_tp2 ~ Normal(9, 2)
etat2 <- as.data.frame(exp(rnorm(1000000, 9, sd = 2)))
colnames(etat2) <- c("y")
etat2den <- ggplot(etat2, aes(x=y)) + 
  geom_density() + xlim(100, 100000) + theme_bw() + 
  xlab(expression(eta[t[0.2]]))
etat2den

#eta_tp2 ~ Normal(9, 2)
tau <- as.data.frame(rhalft(100000, scale = 1, nu = 1))
colnames(tau) <- c("y")
tauden <- ggplot(tau, aes(x=y)) + 
  geom_density() + xlim(0, 50) + theme_bw() + 
  xlab(expression(tau[pi] ~ "," ~ tau[sigma[2]] ~ "," ~ tau[0.20])) + theme(axis.title.y=element_blank())
tauden

cowplot::plot_grid(etapiden, etasig2den, etat2den, tauden, align = "v")



