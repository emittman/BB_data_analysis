#LOO Model Comparison#
library(rstan)
library(loo)

s <- readRDS("MCMC_draws/samples_lor_only3fails.rds")
samp <- extract(s)




#Data for Loo
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
#Need to Add Arguments for Different Models

llfun <- function(i, data, draws){
  #dat <- data[i,]
  model <- data$model
  ll <- vector(length=24000)
  for (j in 1:24000){
    #Get Parameter Estimates  
    mu1 <- draws$mu1[j]
    sigma1 <- draws$sigma1[j]
    log_pi <- draws$log_pi[j]
    mu2 <- draws$mu2[j,model]
    sigma2 <- draws$sigma2[j,model]
        if (data$failed==0){
            #log(p * f1 * (1 - F2) + f2 * (1 - p * F1))
            likenum  <-  log(exp(log_pi + sev_logpdf(data$endtime, mu1, sigma1) +
                               sev_logccdf(data$endtime, mu2, sigma2)) + 
                             exp(sev_logpdf(data$endtime, mu2, sigma2) + 
                               log(1-exp(log_pi + sev_logcdf(data$endtime, mu1, sigma1)))))
    
            #log(1 - p * F1) + log(1 - F2)
            likedem <- log(1-exp(log_pi + sev_logcdf(data$starttime, mu1, sigma1))) + 
            sev_logccdf(data$starttime, mu2, sigma2)
            l <- likenum - likedem                
          }
        else{
        #log(1 - p * F1) + log(1 - F2)
        likenum = log(1-exp(log_pi + sev_logcdf(data$endtime, mu1, sigma1))) + 
        sev_logccdf(data$endtime, mu2, sigma2)
  
        #log((1 - p * F1) * (1 - F2))
        likedem = log(1-exp(log_pi + sev_logcdf(data$starttime, mu1, sigma1))) + 
        sev_logccdf(data$starttime, mu2, sigma2)
    
        l <- likenum - likedem 
          }
      ll[j] <- l
    }
  return(ll)
}

#test <- llfun(1, data = d,draws = samp)

N <- nrow(d)
log_like_mat <- sapply(1:N, function(i) llfun(i,d[i,,drop=FALSE], samp))
loo_with_mat <- loo(log_like_mat)
    
