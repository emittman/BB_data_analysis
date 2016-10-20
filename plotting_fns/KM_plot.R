###################################################################################
#### Suite of functions for plotting K-M estimate of cdf with truncation correction
#### overlaying point estimate of cdf with optional pointwise credible/confidence band
####
### Expects data columns
# starttime(numeric > 0) endtime(numeric > starttime) censored(0 or 1)

qsev = function(x, location = 0, scale = 1) log(-log(1-x))
psev = function(x, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE){
  upper <- exp(-exp(x))
  if(lower.tail){
    if(log.p)
      return(log(1 - upper))
    else
      return(1-upper)
  } else {
    if(log.p)
      return(log(upper))
    else
      return(upper)
  }
}

pgev <- function(x) exp(-exp(-x))
qgev <- function(p) -log(-log(p))

## Transformations for Weibull plot
qsev_trans <- function(){
  trans <- qsev
  inv <- psev
  scales::trans_new("qsev", transform = trans, inverse = inv)
}

## Transformations for Lognormal plot
qnorm_trans <- function(){
  trans <- qnorm
  inv <- pnorm
  scales::trans_new("qnorm", transform = trans, inverse = inv)
}

## Transformations for Frechet plot
qgev_trans <- function(){
  trans <- qgev
  inv <- pgev
  scales::trans_new("qgev", transform = trans, inverse = inv)
}

## Generate sensible breaks for Weibull plot
xbrks <- function(min, max, N = 8, prec = -3) {
  n <- log2(max/min)
  round(2^seq.int(0, n, length.out = N)*min, prec)
}

ybrks <- function(min_p, max_p, model){
  switch(model,
         weibull = {invcdf <- qsev
         cdf <- psev},
         lnorm   = {invcdf <- qnorm
         cdf <- pnorm},
         frechet = {invcdf <- qgev
         cdf <- pgev})
  
  
  quantiles <- seq.int(invcdf(min_p)-.5, invcdf(max_p)+.5, length.out = 7)
  out <- round(cdf(quantiles), 4)
  
  out <- out + c(.0001, rep(0,5), -.0001)
}

my_pweibull <- function(x, location, scale, ...){
  pweibull(x, 1/scale, exp(location), ...)
}

get_tr_adj <- function(starttime, pi, loc1, scl1, loc2, scl2){
  1 - (1 - pi * my_pweibull(starttime, loc1, scl1) * (1 - my_pweibull(starttime, loc2, scl2)))
}

## ggplot object containing Weibull probability plot with point estimate and optional credible band
KM_plot <- function(data, model, tr_adj = 0, title = NULL, linear_axes = FALSE, fixed=FALSE){
  require(ggplot2)
  require(plyr)
  require(dplyr)
  
  switch(model,
         weibull = {   cdf <- my_pweibull
         ytrans <- "qsev"},
         lnorm   = {   cdf <- plnorm
         ytrans <- "qnorm"},
         frechet = {   cdf <- pfrechet
         ytrans <- "qgev"},
         {stop("`model' must be one of lnorm, weibull, frechet")})
  
  # df <- with(data, data.frame(time = endtime[!censored]))
  df <- data.frame(t = sort(unique(data$endtime[!(data$censored)]))) %>%
    ddply(.(t), summarise,
          n = sum(data$starttime < t & data$endtime >=t),
          d = sum(data$starttime < t & data$endtime == t & data$censored == 0)) %>%
    mutate(p = (n-d)/n)
  df$St <- cumprod(df$p)
  df$Ft <- (1 - df$St) * (1 - tr_adj) + tr_adj
  df <- df[which(df$Ft < 1),]
  
  p <- df %>%
    ggplot(aes(t, Ft)) + geom_point()
  
  if(!linear_axes){
    if(!fixed){
      p <- p +
        scale_x_continuous(trans = "log", breaks = xbrks(min(df$t), max(df$t))) +
        scale_y_continuous(trans = ytrans, breaks = ybrks(min(df$Ft), max(df$Ft), model="weibull"))
      
    } else{
      p <- p +
        scale_x_continuous(trans = "log", breaks = xbrks(1000, 50000, prec=0), limits=c(1000, 50000)) +
        scale_y_continuous(trans = ytrans, breaks = ybrks(.0001, .9999, model="weibull"), limits=c(.0001, .9999))
    }
  }
  p <- p +
    theme_bw(base_size = 14) + ggtitle(title)+xlab("Hours")+ylab("Fraction Failing")
  p
}

# # No truncation example
# 
# df <- data.frame(endtime=rweibull(1000, 2))
# df$starttime <- 0
# df$cens <- 0
# 
# band_df <- data.frame(time = seq(min(df$endtime), max(df$endtime), length.out=100)) %>%
#                         mutate(lower=pweibull(time, 2, 4),
#                                upper = pweibull(time, 2, .5))
# 
# p <- KM_plot(df, c(2,1), "Baloney plot", band_df, T)
# p
# 
# # Truncation example
# 
# df <- data.frame(endtime=rweibull(1000, 2))
# df$starttime <- .2
# df <- df[df$endtime>.2, ]
# df$cens <- 0
# 
# band_df <- data.frame(time = seq(min(df$endtime), max(df$endtime), length.out=100)) %>%
#   mutate(lower=pweibull(time, 2, 4),
#          upper = pweibull(time, 2, .5))
# 
# p <- KM_plot(df, c(2,1), "Some Truncated Malarkey", band_df, T)
# p
