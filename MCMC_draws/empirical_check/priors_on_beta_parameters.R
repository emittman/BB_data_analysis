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
