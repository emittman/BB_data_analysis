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

## Transformations for Weibull plot
qsev_trans <- function(){
  trans <- qsev
  inv <- psev
  scales::trans_new("qsev", transform = trans, inverse = inv)
}

my_pweibull <- function(x, location, scale, ...){
  pweibull(x, 1/scale, exp(location), ...)
}

## Compute Kaplan-Meier estimate for possibly left-trunctated,
# right-censored univariate lifetime data (s3 class "myKM")
KM.survfit <- function(life.data, greenwood = FALSE, alpha = .05){
  if(!all(c("starttime","endtime","censored") %in% names(life.data))){
    stop("Check names of life.data; must be a data.frame with 'starttime','endtime','censored'")
  }
  require(plyr)
  require(dplyr)
  df <- data.frame(t = sort(unique(life.data$endtime[!(life.data$censored)]))) %>%
    ddply(.(t), summarise,
          n = sum(life.data$starttime < t & life.data$endtime >=t),
          d = sum(life.data$starttime < t & life.data$endtime == t &
                    life.data$censored == 0)) %>%
    mutate(p = (n-d)/n)
  df$St <- cumprod(df$p)
  df$Ft <- 1 - df$St
  df <- df[which(df$Ft < 1),]
  if(greenwood){
    zcrit <- qnorm(1-alpha)
    df <- mutate(df,
                 Vt = St * cumsum(d/(n*(n-d))),
                 lower = pmax(0, Ft - zcrit*sqrt(Vt)),
                 upper = pmin(1, Ft + zcrit*sqrt(Vt)))
    
  }
  attributes(df)$greenwood <- greenwood
  class(df) <- c("data.frame", "myKM")
  df
}

#Given a KM truncation adjustment estimate (Ft at min(starttime)), adjust myKM object
tr.adj <- function(fit, tr_adj){
  stopifnot(tr_adj>=0)
  if(!("myKM" %in% class(fit))) stop("fit must be output from KM.survfit")
  fit$Ft <- fit$Ft * (1 - tr_adj) + tr_adj
  fit
}

#Given xlimits (for plotting) adjust myKM obj
trimKMx <- function(fit, xlimits){
  stopifnot(length(xlimits)==2, xlimits[1]<xlimits[2])
  if(!("myKM" %in% class(fit))) stop("fit must be output from KM.survfit")
  require(dplyr)
  id <- which(fit$t<xlimits[1])
  if(length(id)>1){
    fit <- fit[-id[-length(id)],]
  }
  fit$t <- pmax(fit$t, xlimits[1])
  
  id <- which(fit$t>xlimits[2])
  if(length(id)>1){
    fit <- fit[-id[-1],]
  }
  fit$t <- pmin(fit$t, xlimits[2])
  fit
}

#Given ylimits (for plotting) adjust myKM obj
trimKMy <- function(fit, ylimits){
  stopifnot(length(ylimits)==2, ylimits[1]<ylimits[2])
  if(!("myKM" %in% class(fit))) stop("fit must be output from KM.survfit")
  id <- which(fit$Ft<ylimits[1])
  if(length(id>0)){
    fit <- fit[-id,]
  }
  id <- which(fit$Ft>ylimits[2])
  if(length(id)>1){
    fit <- fit[-id[-1],]
  }
  fit$Ft <- pmin(fit$Ft, ylimits[2])
  if("Vt" %in% names(fit)){
    fit$lower <- pmax(fit$lower, ylimits[1])
    fit$upper <- pmin(fit$upper, ylimits[2])
  }
  fit
}

baseKMplot.multiple <- function(fit_list, xlimits=NULL, ylimits=NULL, color="black", linetype=1, logscale=NULL, prob=TRUE, label="nonparametric", alpha=1){
  require(ggplot2)
  classes <- sapply(fit_list, function(li) "myKM" %in% class(li))
  if(!(all(classes))) stop("fits must be outputs from KM.survfit")
  if(!is.null(xlimits)){
    fit_list <- lapply(fit_list, function(li) trimKMx(li, xlimits))
  }
  if(!is.null(ylimits)){
    fit_list <- lapply(fit_list, function(li) trimKMy(li, ylimits))
  }
  
  #check to see if any KMs are now empty
  empty_id <- which(sapply(fit_list, nrow)==0)
  if(length(empty_id)>0) {fit_list <- fit_list[ -empty_id]}
  
  fit_list <- lapply(1:length(fit_list), function(i) cbind(fit_list[[i]], rep=i))
  df.combined <- do.call(rbind, fit_list)
  class(df.combined) <- c("data.frame", "myKM")
  df.combined$grp <- as.character(1)
  p <- ggplot(df.combined, aes(x=t, y=Ft, color=grp, linetype=grp)) + geom_step(aes(group=rep), alpha=alpha)
  
  #for scale_color_manual
  scm <- list(colorvals=c("1"=color), linetypevals=c("1"=linetype), labels=label)
  
  out <- list(xlimits=xlimits, ylimits=ylimits, base=p, scm=scm, logscale=logscale, prob=prob)
  attributes(out)$greenwood <- attr(df.combined, "green")
  class(out) <- "myPlotList"
  return(out)
}

#Make base KMplot
baseKMplot <- function(fit, xlimits=NULL, ylimits=NULL, color="black", linetype=1, logscale=NULL, prob=TRUE, label="nonparametric", alpha=1){
  require(ggplot2)
  if(!("myKM" %in% class(fit))) stop("fit must be output from KM.survfit")
  if(!is.null(xlimits)){
    fit <- trimKMx(fit, xlimits)
  }
  if(!is.null(ylimits)){
    fit <- trimKMy(fit, ylimits)
  }
  
  fit$grp <- as.character(1)
  p <- ggplot(fit, aes(x=t, y=Ft, color=grp, linetype=grp)) + geom_step()

  #for scale_color_manual
  scm <- list(colorvals=c("1"=color), linetypevals=c("1"=linetype), labels=label)
  
  out <- list(xlimits=xlimits, ylimits=ylimits, base=p, scm=scm, logscale=logscale, prob=prob)
  attributes(out)$greenwood <- attr(fit, "green")
  class(out) <- "myPlotList"
  return(out)
}

#Compute data.frame containing pointwise posterior estimate of lifetime for GLFP model
bandFromPSamp <- function(samples, range, length.out=25, N=NULL, mean=F, quantiles=c(.05,.95), logscale=TRUE){
  #check inputs
  stopifnot(all(c("mu1","mu2","sigma1","sigma2","log_pi") %in% names(samples)))
  S <- length(samples$mu1)
  stopifnot(length(range)==2, length(quantiles)==2, length.out>1)
  if(is.null(N)){ N <- floor(length(samples$mu1)/4)}
  stopifnot(N<=S)
  
  #grid of times
  if(logscale){
    t <- exp(seq(log(range[1]), log(range[2]), length.out=length.out))
  } else{
    t <- seq(range[1], range[2], length.out=length.out)
  }
  #representative sample of posterior samples
  id <- sample(1:S, N)
  vals <- sapply(t, function(x){
    p <- sapply(id, function(i) {
      1 -  (1 - exp(samples$log_pi[i]) * my_pweibull(x, samples$mu1[i], samples$sigma1[i])) *
        (1 - my_pweibull(x, samples$mu2[i], samples$sigma2[i]))
    })
  })
  if(mean) {
    summaries <- apply(vals, 2, function(col){
    q <- quantile(col, c(quantiles[1], quantiles[2]), names=FALSE)
    mean <- mean(col)
    return(c(lower=q[1], est=mean, upper=q[2]))
    })
  } else{
    summaries <- apply(vals, 2, function(col){
      q <- quantile(col, c(quantiles[1], .5, quantiles[2]), names=FALSE)
      return(c(lower=q[1], est=q[2], upper=q[3]))
    })
  }
  df <- cbind(data.frame(t=t), as.data.frame(t(summaries)))
  out <- list(band=df, logscale=logscale)
  class(out) <- c("myBand")
  out
}

bandTrimy <- function(bandObj, ylimits){
  stopifnot("myBand" %in% class(bandObj))
  #these are removed
  id <- which(bandObj$band$upper<ylimits[1] | bandObj$band$lower > ylimits[2])
  if(length(id)){
    bandObj$band <- bandObj$band[-id,]
  }
  #trim y
  bandObj$band$lower<- pmax(ylimits[1], bandObj$band$lower)
  bandObj$band$est <- pmax(ylimits[1], bandObj$band$est)
  bandObj$band$est <- pmin(ylimits[2], bandObj$band$est)
  bandObj$band$upper<- pmin(ylimits[2], bandObj$band$upper)
  bandObj
}

bandTrimx <- function(bandObj, xlimits){
  stopifnot("myBand" %in% class(bandObj))
  id <- which(bandObj$band$t > xlimits[2] | bandObj$band$t < xlimits[1])
  if(length(id)>0) {bandObj$band <- bandObj$band[-id,]}
  bandObj
}

addKmToBaseplot <- function(baseplot, fitObj, color, linetype, size=1, label=""){
  stopifnot(class(baseplot) == "myPlotList")
  stopifnot("myKM" %in% class(fitObj))
  if(!is.null(baseplot$xlimits)){
    fitObj <- trimKMx(fitObj, baseplot$xlimits)
  }
  if(!is.null(baseplot$ylimits)){
    fitObj <- trimKMy(fitObj, baseplot$ylimits)
  }
  grp_num <- as.character(length(baseplot$scm$colorvals)+1)
  fitObj$grp <- grp_num
  baseplot$base <- baseplot$base +
    geom_step(data=fitObj, inherit.aes=FALSE, aes(x=t, y=Ft, color=grp, linetype=grp), size=size)
  
  concat_names <- c(names(baseplot$scm$colorvals), grp_num)
  baseplot$scm$colorvals <- with(baseplot$scm, c(colorvals, color))
  baseplot$scm$linetypevals <- with(baseplot$scm, c(linetypevals, linetype))
  names(baseplot$scm$colorvals) <- concat_names
  names(baseplot$scm$linetypevals) <- concat_names
  baseplot$scm$labels <- c(baseplot$scm$labels, label)
  
  baseplot
}

addBandToBaseplot <- function(baseplot, bandObj, color, linetype, alpha=.5, label=""){
  stopifnot(class(baseplot) == "myPlotList")
  stopifnot(is.null(baseplot$logscale) || baseplot$logscale == bandObj$logscale)
  if(!is.null(baseplot$xlimits)){
    bandObj <- bandTrimx(bandObj, baseplot$xlimits)
  }
  if(!is.null(baseplot$ylimits)){
    bandObj <- bandTrimy(bandObj, baseplot$ylimits)
  }
  
  grp_num <- as.character(length(baseplot$scm$colorvals)+1)
  bandObj$band$grp <- grp_num
  
  baseplot$base <- baseplot$base + 
    geom_ribbon(data=bandObj$band, inherit.aes=FALSE, aes(x=t, ymin=lower, ymax=upper, fill=grp), alpha=alpha)+
    geom_line(data=bandObj$band, inherit.aes=FALSE, aes(x=t, y=est, color=grp, linetype=grp))
  
  concat_names <- c(names(baseplot$scm$colorvals), grp_num)
  baseplot$scm$colorvals <- with(baseplot$scm, c(colorvals, color))
  baseplot$scm$linetypevals <- with(baseplot$scm, c(linetypevals, linetype))
  names(baseplot$scm$colorvals) <- concat_names
  names(baseplot$scm$linetypevals) <- concat_names
  baseplot$scm$labels <- c(baseplot$scm$labels, label)

  baseplot$logscale <- bandObj$logscale
  
  baseplot
}

plotFinally <- function(plotList, xbrks, ybrks, years=FALSE, greenwood=FALSE, alpha=.2){
  stopifnot(!greenwood || attr(plotList, "greenwood"))
  
  p <- plotList$base
  
  if(greenwood){
    p <- p + geom_ribbon(inherit.aes = F, aes(x=t, ymin=lower, ymax=upper, fill=grp), alpha=alpha)
  }
  
  
  if(plotList$logscale & years){
    p <- p + scale_x_continuous(name="Time (Years)", trans="log",
                                            breaks = xbrks, labels=signif(xbrks/365/24, 0), limits = plotList$xlimits)
  } else if(plotList$logscale & !years) {
    p <- p + scale_x_continuous(name="Thousands of hours", trans="log", labels = xbrks/1000,
                                            breaks = xbrks, limits = plotList$xlimits)
  } else if(years){
    p <- p + scale_x_continuous(name="Time (Years)", breaks = xbrks,
                                            labels=signif(xbrks/365/24, 0), limits = plotList$xlimits)
  } else{
    p <- p + scale_x_continuous(name="Thousands of hours", breaks = xbrks, labels=xbrks/1000
                                            , limits = plotList$xlimits)
  }
  
  if(plotList$prob){
    p <- p + scale_y_continuous(name="Proportion failing", trans="qsev", breaks = ybrks,
                           limits = plotList$ylimits)
  } else{
    p <- p + scale_y_continuous(name="Proportion failing", breaks = ybrks, limits = plotList$ylimits)
  }
  
  p + scale_color_manual(name="", values = plotList$scm$colorvals, labels = plotList$scm$labels)+
  scale_linetype_manual(name="", values = plotList$scm$linetypevals, labels = plotList$scm$labels)+
    coord_cartesian(xlim = plotList$scm$xlimits, ylim = plotList$scm$ylimits, expand = F)+
    scale_fill_manual(name="", values=plotList$scm$colorvals, labels=plotList$scm$labels)
}
