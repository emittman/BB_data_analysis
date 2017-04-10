KM_hazard <- function(data, linear_axes = FALSE, xlimits=c(1000, 50000),
                    ylimits=c(.0001, .001), conf=.95, verbose=F){
  require(ggplot2)
  require(plyr)
  require(dplyr)
  
zcrit <- qnorm(1-(1-conf)/2)
  # df <- with(data, data.frame(time = endtime[!censored]))
  df <- data.frame(t = sort(unique(data$endtime[!(data$censored)]))) %>%
    ddply(.(t), summarise,
          n = as.numeric(sum(data$starttime < t & data$endtime >=t)),
          d = as.numeric(sum(data$starttime < t & data$endtime == t & data$censored == 0))) %>%
    mutate(p = 1 -(n-d)/n)
  
  
  if(verbose){
    print(df)
  }
  
  p <- df %>%
    ggplot(aes(t, p)) + geom_line(color="orange") + 
    scale_x_continuous(trans="log", limits=xlimits) + 
    scale_y_continuous(trans="log", limits=ylimits)

  # if(!linear_axes){
  #   p <- p +
  #     scale_x_continuous(trans = "log")+#, breaks = xbrks(xlimits[1], xlimits[2], prec=0), limits=xlimits) +
  #     scale_y_continuous(trans = "log")#, breaks = xbrks(min(df$p), max(df$p), prec=4), limits=ylimits)
  # } else{
  #   p <- p + scale_x_continuous(limits=xlimits) + scale_y_continuous(limits=ylimits)
  # }
  # p <- p +
  #   theme_bw(base_size = 14) + ggtitle(title)+xlab("Hours")+ylab("Hazard Rate")
  # p
}