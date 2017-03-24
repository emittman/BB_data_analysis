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
KM_plot <- function(data, model, tr_adj = 0, title = NULL,
                    linear_axes = FALSE, fixed=FALSE, xlimits=c(1000, 50000), ylimits=c(.0001, .9999), verbose=F){
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
  
  if(verbose){
    print(df)
  }
  
  p <- df %>%
    ggplot(aes(t, Ft)) + geom_point()
  
  if(!linear_axes){
    if(!fixed){
      p <- p +
        scale_x_continuous(trans = "log", breaks = xbrks(min(df$t), max(df$t))) +
        scale_y_continuous(trans = ytrans, breaks = ybrks(min(df$Ft), max(df$Ft), model="weibull"))
      
    } else{
      p <- p +
        scale_x_continuous(trans = "log", breaks = xbrks(xlimits[1], xlimits[2], prec=0), limits=xlimits) +
        scale_y_continuous(trans = ytrans, breaks = ybrks(ylimits[1], ylimits[2], model="weibull"), limits=ylimits)
    }
  } else{
    p <- p + scale_x_continuous(limits=xlimits) + scale_y_continuous(limits=ylimits)
  }
  p <- p +
    theme_bw(base_size = 14) + ggtitle(title)+xlab("Hours")+ylab("Fraction Failing")
  p
}

KM_band <- function(id, samp, n_iter=NULL, xlim, ylim, quantiles=c(.05, .5, .95), x_logscale=T, n=30, verbose=F){
  total_iter <- dim(samp[[1]])[1]
  iter <- floor(total_iter/n_iter) * 1:n_iter
  
  if(x_logscale){
    x <- exp(seq(log(xlim[1]), log(xlim[2]), length.out=n))
  } else{
    x <- seq(xlim[1], xlim[2], length.out=n)
  }
  
  band <- data.frame(x=x) %>%
    ddply(.(x), function(g){
      Fp <- sapply(iter, function(i) {
        1 -  (1 - exp(samp$log_pi[i+total_iter*(id-1)]) * my_pweibull(g$x, samp$mu1[id], samp$sigma1[id])) *
          (1 - my_pweibull(g$x, samp$mu2[i+total_iter*(id-1)], samp$sigma2[i+total_iter*(id-1)]))
      })
      q <- quantile(Fp, quantiles)
      data.frame(y = q[2], lower = max(ylim[1],q[1]), upper = min(ylim[2],q[3]))
    })
  if(verbose) print(band)
  geom_ribbon(data = band, inherit.aes = FALSE, aes(x=x, ymin=lower, ymax=upper), fill="red", alpha=.2)
}

KM_with_band <- function(title = NULL, data, id, samp, n_iter, n, quantiles, tr_adj, xlimits, ylimits, fixed=T, linear_axes=F, verbose=F, model="weibull"){
  p <- KM_plot(data = data, model = model, tr_adj = tr_adj, title=title, fixed=fixed, linear_axes = linear_axes, xlimits=xlimits, ylimits=ylimits, verbose)
  q <- KM_band(id, samp, n_iter, xlimits, ylimits, quantiles, !linear_axes, n, verbose)
  p + q
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
