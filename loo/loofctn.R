#LOO Model Comparison#
library(rstan)

s <- readRDS("MCMC_draws/samples_lor_only3fails.rds")
samp <- extract(s)


#Make Matrix of Draws to Pass to Function


d <- readRDS("BB_data/clean_unit_summaries.rds")
d$model <- as.integer(d$model)
d$starttime <- log(d$starttime + 1)
d$endtime <- log(d$endtime + 1)

#Functions for log sev CDF and log sev PDF
sev_logpdf<- function(y, mu, sigma){
    z = (y - mu) / sigma
    return (-log(sigma) + z - exp(z))
}

sev_logcdf <- function(y, mu, sigma){
  out = (-exp((y - mu) / sigma))
  return(log(1-exp(out)))
}

sev_logccdf <- function(y, mu, sigma){
  return (-exp((y - mu) / sigma))
}


#Function to Take an Observation and Return a Vector of the Likelihood over all Posterior Parameters

llfun <- function(i, data, draws){
  dat <- data[i,]
  model <- dat$model
  ll <- vector(length=10)
  for (j in 1:10){
  #Get Parameter Estimates  
  mu1 <- draws$mu1[j]
  sigma1 <- draws$sigma1[j]
  log_pi <- draws$log_pi[j]
  mu2 <- draws$mu2[j,model]
  sigma2 <- draws$sigma2[j,model]
  if (dat$failed==0){
    #log(p * f1 * (1 - F2) + f2 * (1 - p * F1))
    likenum  <-  log(exp(log_pi + sev_logpdf(dat$endtime, mu1, sigma1) +
                               sev_logccdf(dat$endtime, mu2, sigma2)) + 
                             exp(sev_logpdf(dat$endtime, mu2, sigma2) + 
                               log(1-exp(log_pi + sev_logcdf(dat$endtime, mu1, sigma1)))))
    
    #log(1 - p * F1) + log(1 - F2)
    likedem <- log(1-exp(log_pi + sev_logcdf(dat$starttime, mu1, sigma1))) + 
      sev_logccdf(dat$starttime, mu2, sigma2)
    l <- likenum - likedem                
  }
  
  else{
    
    #log(1 - p * F1) + log(1 - F2)
    likenum = log(1-exp(log_pi + sev_logcdf(dat$endtime, mu1, sigma1))) + 
      sev_logccdf(dat$endtime, mu2, sigma2)
    
    #log((1 - p * F1) * (1 - F2))
    likedem = log(1-exp(log_pi + sev_logcdf(dat$starttime, mu1, sigma1))) + 
      sev_logccdf(dat$starttime, mu2, sigma2)
    
    l <- likenum - likedem 
  }
  ll[j] <- l
  }
  return(ll)
}

test <- llfun(1, data = d,draws = samp)
    
