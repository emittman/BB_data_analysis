library(rstan)

sfull <- readRDS("workflow/samples_lor_only3fails.rds")
fullsamp <- extract(sfull)

#3, 8, 14 correspond to 6, 14, 23 in full analys.
s6 <- readRDS("workflow2/samples_nopool_data_dm3.rds")
samp6 <- extract(s6)

hist(samp6$log_tp2)
hist(fullsamp$log_tp2_raw[,6] * fullsamp$tau_tp2 + fullsamp$eta_tp2, add=TRUE)
plot(density(samp6$log_tp2))
lines(density(fullsamp$log_tp2_raw[,6] * fullsamp$tau_tp2 + fullsamp$eta_tp2), lty=2)
curve(dnorm(x, 10, 4), col=2, add=TRUE)

s14 <- readRDS("workflow2/samples_nopool_data_dm8.rds")
samp14 <- extract(s14)


str(fullsamp)
str(samp14)

hist(samp14$log_tp2)
hist(fullsamp$log_tp2_raw[,14] * fullsamp$tau_tp2 + fullsamp$eta_tp2, add=TRUE)
plot(density(samp14$log_tp2))
lines(density(fullsamp$log_tp2_raw[,14] * fullsamp$tau_tp2 + fullsamp$eta_tp2), lty=2)
curve(dnorm(x, 10, 4), col=2, add=TRUE)

s23 <- readRDS("workflow2/samples_nopool_data_dm14.rds")
samp23 <- extract(s23)


str(fullsamp)
str(samp23)

hist(samp23$log_tp2)
hist(fullsamp$log_tp2_raw[,23] * fullsamp$tau_tp2 + fullsamp$eta_tp2, add=TRUE)
plot(density(samp23$log_tp2))
lines(density(fullsamp$log_tp2_raw[,23] * fullsamp$tau_tp2 + fullsamp$eta_tp2), lty=2)
curve(dnorm(x, 10, 4), col=2, add=TRUE)


plot(density(samp6$sigma2))
lines(density(fullsamp$sigma2[,6]), lty=2)
curve(dlnorm(x, 0, 2.5), col=2, add=TRUE)

plot(density(samp6$logit_pi))
lines(density(fullsamp$logit_pi_raw[,6]*fullsamp$tau_pi + fullsamp$eta_pi), lty=2)
curve(dnorm(x, -3, 2), col=2, add=TRUE)

plot(density(samp14$sigma2))
lines(density(fullsamp$sigma2[,14]), lty=2)
curve(dlnorm(x, 0, 2.5), col=2, add=TRUE)

plot(density(samp14$logit_pi))
lines(density(fullsamp$logit_pi_raw[,14]*fullsamp$tau_pi + fullsamp$eta_pi), lty=2)
curve(dnorm(x, -3, 2), col=2, add=TRUE)
