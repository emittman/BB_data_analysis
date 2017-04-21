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

llfun <- function(i, data, draws, which.model){
  # which.model is logical vector (pi.free, mu1.free, sigma1.free, mu2.free, sigma2.free)
  #dat <- data[i,]
  niter <- length(draws$lp__)
  col_idx <- n_iter*(data$model-1)
  model <- data$model
  idata <- data[i,]
  if(data$failed==0){
    likenum <- with(which.model,
                  log(exp(draws$logpi[1:n_iter + col_idx*pi.free] +
                  sev_logpdf(idata$endtime, draws$mu1[1:n_iter + col_idx*mu1.free], draws$sigma1[1:n_iter + col_idx*sigma1.free]) +
                  sev_logccdf(idata$endtime, draws$mu2[1:n_iter + col_idx*mu2.free], draws$sigma2[1:n_iter + col_idx*sigma2.free])) +
                  exp(sev_logpdf(idata$endtime, draws$mu2[1:n_iter + col_idx*mu2.free], draws$sigma2[1:n_iter + col_idx*sigma2.free])+
                  log(1 - exp(log_pi + sev_logcdf(idata$endtime, draws$mu1[1:n_iter + col_idx*mu1.free], draws$sigma1[1:n_iter + col_idx*sigma1.free]))))))
    likedem <- with(which.model, log(1)
                    )
  } else{
    
  }
  return(likenum-likedem)
}
  
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


#test <- llfun(1, data = d,draws = samp)

N <- nrow(d)
log_like_mat <- sapply(1:N, function(i) llfun(i,d[i,,drop=FALSE], samp))

#this should be
#loo_output <- loo(llfun, args = list(data=d, N=N, S=S, draws=samp, which.model=c(T, F, F, T, T)))
loo_with_mat <- loo(log_like_mat)
    
