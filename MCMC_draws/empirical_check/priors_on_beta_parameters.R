df <- data.frame(lphi = rnorm(1000, 0, 5),
                 mu = runif(1000))

df$a <- with(df, mu*exp(lphi))
df$b <- with(df, (1-mu)*exp(lphi))

library(ggplot2)

p <- ggplot(data = data.frame(x=c(0.001,.999)), aes(x))

for(i in 1:100){
  p <- p + stat_function(fun = dbeta, args = list(shape1=df$a[i], shape2=df$b[i]), alpha=.2)
}

p + ylim(c(0,1.25))

f <- function(x) x*(1-x) / (exp(2)+1)
curve(sqrt(f(x)), from=0, to=.25)

#log-odds delta where pi ~ beta
b <- rbeta(10000, .4, 2.85)
lo <- log(b/(1-b))
par(mfrow=c(2,1))
hist(b, 100)
hist(lo, 100)
mu = median(lo)
sigma = sd(lo)
hist(lo, 100, prob=T)
curve(dnorm(x, mu, sigma), add=T, lty=2)
ggplot(data = data.frame(x=b), aes(x=x))+
    geom_histogram(aes(y=..density..))+
    geom_segment(data=data.frame(x=0, y=-.5, xend=qbeta(.999, .4, 2.85), yend=-.5),
                 aes(x=x,y=y,xend=xend,yend=yend))
ggplot(data = data.frame(x=lo), aes(x=x))+
    geom_histogram(aes(y=..density..))+
    geom_segment(data=data.frame(x=qnorm(.001, mu, sigma), xend=qnorm(.999, mu, sigma), y=-.01, yend=-.01),
                                 aes(x=x,y=y,xend=xend,yend=yend))
