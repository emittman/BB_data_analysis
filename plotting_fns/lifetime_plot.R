lifetime_plot <- function(data, n_to_show=NULL, in_years=TRUE, lab = "", trans="identity",
                          xlabels=c(.01,.05,.25,1,2,5), xlimits=c(100, 50000)){
  if(!in_years){
    xlabels <- signif(xlabels * 24*365, 0)
  } 
  if(!is.null(n_to_show)){
    N = nrow(data)
    data <- sample_n(data, size=n_to_show)
  }
  data <- data %>%
    mutate(ID=as.integer(factor(serial_number, levels=serial_number[order(starttime)])),
           starttime=pmax(starttime, xlimits[1]))
  
  p <- filter(data, censored) %>%
    ggplot(aes(x=ID, xend=ID, y=starttime, yend=endtime)) +
    geom_segment(arrow = arrow(length=unit(.15, "cm"), type = "open"), alpha=.5)+
    # geom_segment(alpha=.5)+
    # geom_point(aes(y=endtime), color="blue", alpha=.5) +
    geom_segment(data = filter(data, !censored),
                 aes(x=ID, xend=ID, y=starttime, yend=endtime), alpha=.5)+
    geom_point(data = filter(data, !censored),
               aes(x=ID,y=endtime), shape=4, size=3) +
    # geom_point(data = filter(data, !censored),
               # aes(x=ID, y=endtime), color="red", alpha=.5) +
    coord_flip() +
    theme_bw() +
    theme(axis.title.y= element_blank(),
          axis.ticks.y= element_blank(),
          axis.text.y = element_blank())
  
  transform <- switch(trans,
                      "log" = scales::log_trans,
                      "identity" = scales::identity_trans)
  
  if(in_years){
  p <-   p + scale_y_continuous(name="time(years)",limits = xlimits,
                           breaks=xlabels*365*24,
                           labels=xlabels,
                           trans=trans)
  } else{
  p <-   p + scale_y_continuous(name="time(years)",limits = xlimits,
                           breaks=xlabels,
                           labels=xlabels,
                           trans=trans)
  }
    
  
  if(is.null(n_to_show)) return(p+ggtitle(paste(mod)))
  p + ggtitle(paste(lab, "(sample size = ", n_to_show, ", N = ", N, ")"))
}

lifetime_plot2 <- function(data, n_to_show=NULL, in_years=TRUE, lab = "", trans="identity",
                       xlabels=c(.01,.05,.25,1,2,5), xlimits=c(100, 50000)){
  if(!in_years){
    xlabels <- signif(xlabels * 24*365, 0)
  } 
  if(!is.null(n_to_show)){
    N = nrow(data)
    data <- sample_n(data, size=n_to_show)
  }
  n_cens <- sum(data$censored)
  n_failed <- sum(!data$censored)
  data <- data %>% arrange(censored,starttime) %>%
    mutate(ID=c(1:n_failed, 1:n_cens),
           status=factor(ifelse(censored, "censored","failed")),
           starttime=pmax(starttime, xlimits[1]))
  
  p <- ggplot(data, aes(x=ID, xend=ID, y=starttime, yend=endtime)) +
    geom_segment(alpha=.5, color="white")+
    facet_grid(.~status)+
    coord_flip() +
    theme_dark() +
    theme(axis.title.y= element_blank(),
          axis.ticks.y= element_blank(),
          axis.text.y = element_blank())
  
  transform <- switch(trans,
                      "log" = scales::log_trans,
                      "identity" = scales::identity_trans)
  
  if(in_years){
    p <-   p + scale_y_continuous(name="time(years)",limits = xlimits,
                                  breaks=xlabels*365*24,
                                  labels=xlabels,
                                  trans=trans)
  } else{
    p <-   p + scale_y_continuous(name="time(years)",limits = xlimits,
                                  breaks=xlabels,
                                  labels=xlabels,
                                  trans=trans)
  }
  
  
  if(is.null(n_to_show)) return(p+ggtitle(paste(mod)))
  p + ggtitle(paste(lab, "(sample size = ", n_to_show, ", N = ", N, ")"))
}

lifetime_plot3 <- function(data, n_to_show=NULL, in_years=TRUE, lab = "", trans="identity",
                           xlabels=c(.01,.05,.25,1,2,5), xlimits=c(100, 50000), font_size=10){
  if(!in_years){
    xlabels <- signif(xlabels * 24*365, 0)
  } 
  if(!is.null(n_to_show)){
    N = nrow(data) #referenced in plot main title
    data <- sample_n(data, size=n_to_show)
  }
  n <- nrow(data)
  data <- data %>% arrange(starttime) %>%
    mutate(ID=1:n)

  data_c <- filter(data, censored==TRUE) %>%
    mutate(status="censored",
           starttime=endtime,
           endtime=xlimits[2])
  
  # data_f <- filter(data, censored==FALSE) %>%
  #   mutate(starttime=endtime,
  #          endtime=xlimits[2],
  #          status="failed")
  
  data_o <- mutate(data,
                   starttime=pmax(starttime, xlimits[1]),
                   status="observed")
  
  data_new <- rbind(data_c,data_o)
  data_new$status <- factor(data_new$status)
  
  p <- ggplot(data_new, aes(x=ID, xend=ID, y=starttime, yend=endtime, color=status, linetype=status, alpha=status)) +
    geom_segment()+
    coord_flip(expand=F) +
    theme_bw(base_size=font_size) +
    theme(axis.title.y= element_blank(),
          axis.ticks.y= element_blank(),
          axis.text.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank()) +
    scale_x_continuous(limits=c(0, n+1))
  
  transform <- switch(trans,
                      "log" = scales::log_trans,
                      "identity" = scales::identity_trans)
  
  if(in_years){
    p <-   p + scale_y_continuous(name="time(years)", limits=xlimits,
                                  breaks=xlabels*365*24,
                                  labels=xlabels,
                                  trans=trans)
  } else{
    p <-   p + scale_y_continuous(name="Thousands of hours",limits = xlimits,
                                  breaks=xlabels*1000,
                                  labels=xlabels,
                                  trans=trans)
  }
  p <- p + scale_color_manual(values = c("censored"="black", "observed"="black"))+
    scale_linetype_manual(values = c("censored"="dashed", "observed"="solid"))+
    scale_alpha_manual(values = c("censored"="1", "observed"="1"))
  
  if(is.null(n_to_show)) return(p+ggtitle(paste(lab)))
  p + ggtitle(paste(lab, "(sample size = ", n_to_show, ", N = ", N, ")"))
}
