KM_plot_NP <- function(data, model, tr_adj = 0, title = NULL,
                    linear_axes = FALSE, fixed=FALSE, xlimits=c(1000, 50000),
                    ylimits=c(.0001, .9999), conf=.95, verbose=F){
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
  
  zcrit <- qnorm(1-(1-conf)/2)
  # df <- with(data, data.frame(time = endtime[!censored]))
  df <- data.frame(t = sort(unique(data$endtime[!(data$censored)]))) %>%
    ddply(.(t), summarise,
          n = as.numeric(sum(data$starttime < t & data$endtime >=t)),
          d = as.numeric(sum(data$starttime < t & data$endtime == t & data$censored == 0))) %>%
    mutate(p = (n-d)/n)
  df$St <- cumprod(df$p)
  df$Vt <- with(df, St * cumsum(d/(n*(n-d))))
  df$lowS <- with(df, St - zcrit*sqrt(Vt))
  df$uprS <- with(df, St + zcrit*sqrt(Vt))
  
  df$Ft <- (1 - df$St) * (1 - tr_adj) + tr_adj
  df$lowFt <- pmax((1 - df$uprS) * (1 - tr_adj) + tr_adj,1e-5)
  df$uprFt <- pmin((1 - df$lowS) * (1 - tr_adj) + tr_adj,1 - 1e-5)
  df <- df[which(df$Ft < 1),]
  
  
  if(verbose){
    print(df)
  }
  
  p <- df %>%
    ggplot(aes(t, Ft)) + geom_point() +
    geom_line(aes(y=lowFt), color="green", lty=2)+
    geom_line(aes(y=uprFt), color="green", lty=2)
    
  
  if(!linear_axes){
    if(!fixed){
      p <- p +
        scale_x_continuous(trans = "log", breaks = xbrks(min(df$t), max(df$t)), limits=xlimits) +
        scale_y_continuous(trans = ytrans, breaks = ybrks(min(df$Ft), max(df$Ft), model="weibull"), limits=ylimits)
      
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