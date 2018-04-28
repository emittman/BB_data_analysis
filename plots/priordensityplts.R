library(ggplot2)
library(LaplacesDemon)

## Priors for Model 1  ##
formatter1000 <- function(x){ 
  x/1000 
}


#t, .5, 1 ~ Normal(8, 4)
t1 <- as.data.frame(exp(rnorm(100000, 8, sd=4)))
colnames(t1) <- c("y")
t1den <- ggplot(t1, aes(x=y)) + 
  geom_density()  + theme_bw() + 
  xlab(expression(t[0.5])) + scale_x_continuous(labels = formatter1000, limits = c(0, 100000)) + ylab("Density") + 
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
t1den

#t, .2, 2 ~ Normal(10, 4)
t2 <- as.data.frame(exp(rnorm(100000, 10, sd=4)))
colnames(t2) <- c("y")
t2den <- ggplot(t2, aes(x=y)) + 
  geom_density() + theme_bw() + scale_x_continuous(labels = formatter1000, limits = c(0, 100000))  + 
  xlab(expression(t[0.2]))  + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
t2den

#InvLogitPi ~ Normal(-3, 2)
pi <- as.data.frame(rnorm(100000, -3, sd=2))
pi.tran <- exp(pi)/(1+exp(pi))
colnames(pi.tran) <- c("y")
piden <- ggplot(pi.tran, aes(x=y)) + 
  geom_density() + xlim(0, 1) + theme_bw() +   theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + ylab("Density") + 
  xlab(expression(pi))
piden

#Sigma Terms (Log Normal (0, 2.5))
sig <- as.data.frame(rlnorm(100000, 0, sd=2.5))
colnames(sig) <- c("y")
sigden <- ggplot(sig, aes(x=y)) + 
  geom_density() + xlim(0, 10) + theme_bw() + 
  xlab(expression(sigma[1] ~ "," ~ sigma[2])) + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
sigden


#Note: x-axis is truncated so the plots don't just look like a spike
cowplot::plot_grid(t1den, t2den, piden, sigden, align = "h")


## Priors for Model 2 ##
#t, .5, 1 ~ Normal(7, 2)
t1m2 <- as.data.frame(exp(rnorm(100000, 7, sd=2)))
colnames(t1m2) <- c("y")
t1denm2 <- ggplot(t1m2, aes(x=y)) + 
  geom_density() + theme_bw() + scale_x_continuous(labels = formatter1000, limits = c(0, 100000))  + 
  xlab(expression(t[0.5]))  + theme(axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.text.y = element_blank())
t1denm2


#Sigma 1 Term ~ Log Normal (0, 1)
sig1 <- as.data.frame(rlnorm(100000, 0, sd=1))
colnames(sig1) <- c("y")
sigdenm2 <- ggplot(sig1, aes(x=y)) + 
  geom_density() + xlim(0, 10) + theme_bw() + ylab("Density") + 
  xlab(expression(sigma[1])) +  theme(axis.ticks.y = element_blank(), axis.text.y=element_blank())
sigdenm2

#t_p2_g ~ Lognormal(eta_tp2, tau_tp2) eta_tp2 ~ Normal(9, 2), tau_sig2 ~ Ca(0, 1)
tau_tp2 <- rhalft(1000000, scale = 1, nu = 1)
eta_tp2 <- rnorm(1000000, mean = 9, sd = 2)

out <- numeric(1000000)
for(i in 1:1000000){
  out[i] <- rnorm(1, mean = eta_tp2[i], sd = tau_tp2[i])
}

tp2g.tran <- as.data.frame(exp(out))

colnames(tp2g.tran) <- c("y")
tp2.g.den <- ggplot(tp2g.tran, aes(x=y)) + 
  geom_density()  + theme_bw() + scale_x_continuous(labels = formatter1000, limits = c(0, 100000)) + 
  xlab(expression(t[0.2][g])) + theme(axis.title.y=element_blank(), axis.ticks.y = element_blank(), axis.text.y=element_blank()) + 
  ylab("Density")
tp2.g.den



## Priors for Model 3 ##

#sigma2_g ~ Lognormal(eta_sig2, tau_sig2) eta_sig2 ~ Normal(0, 2), tau_sig2 ~ Ca(0, 1).  Note: This is Truncated T(0, 1)
tau_sig2 <- rhalft(1000000, scale = 1, nu = 1)
eta_sig2 <- rnorm(1000000, mean = 0, sd = 2)

out <- numeric(1000000)
for(i in 1:1000000){
  out[i] <- rnorm(1, mean = eta_sig2[i], sd = tau_sig2[i])
}

sig2g.tran <- as.data.frame(exp(out))

colnames(sig2g.tran) <- c("y")
sig2.tran.truncate <- as.data.frame(sig2g.tran[sig2g.tran$y < 1,])
colnames(sig2.tran.truncate) <- c("y")
sig2.g.den <- ggplot(sig2.tran.truncate, aes(x=y)) + 
  geom_density() + xlim(0, 1) + theme_bw() + 
  xlab(expression(sigma[2][g])) + theme(axis.ticks.y = element_blank(), axis.text.y=element_blank()) + 
  ylab("Density")
sig2.g.den


## Priors for Model 4 ##

#pi_g  Logit-Normal(eta_pi, tau_pi); eta_pi ~ N(-3,1), tau_pi ~ Ca(0, 1)
tau_pi <- rhalft(1000000, scale = 1, nu = 1)
eta_pi <- rnorm(1000000, mean = -3, sd = 1)

out <- numeric(1000000)
for(i in 1:1000000){
out[i] <- rnorm(1, mean = eta_pi[i], sd = tau_pi[i])
}

pi.tran <- as.data.frame((exp(out))/(1 + exp(out)))

colnames(pi.tran) <- c("y")
pi.g.den <- ggplot(pi.tran, aes(x=y)) + 
  geom_density() + xlim(0, 1) + theme_bw() + 
  xlab(expression(pi[g])) + theme(axis.title.y=element_blank(), axis.ticks.y = element_blank(), axis.text.y=element_blank()) + 
  ylab("Density")
pi.g.den


#Put Priors for Models 2, 3, and 4 together.

cowplot::plot_grid(sigdenm2, t1denm2, tp2.g.den, sig2.g.den, pi.g.den, align = "h")



