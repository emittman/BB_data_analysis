#LOO Model Comparison#
library(rstan)
library(loo)
setwd("C:/Users/Colin/Documents/GitHub/BB_data_analysis/MCMC_draws")


#Differnt Posterior Draws for Each Model
s1 <- readRDS("samples_null_model_3_29.rds") # Model 1  (all fixed)
s2 <- readRDS("samples_tp2_vary_new.rds")    #Model 2 (mu2 free)
s3 <- readRDS("vary_s2_and_tp2_4_17.rds") # Model 3 (mu2 free, sigma2 free, pi fixed)
s4 <- readRDS("samples_lor_only3fails.rds") # Is this Model 4 (mu2 free, sigma2 free, pi free)



#Data for Loo: Need To Make Sure We Grab the Same Data Used to Fit the Model as We are Not Using All The Models!
source("../workflow/functions.R")

data_all = prepare_data(lb_fails = 3, lb_late_fails = 0, lb_early_fails = 0) #Same Data Set for All Models

data_all$starttime <- log(data_all$starttime + 1)  
data_all$endtime <- log(data_all$endtime + 1)

#Function to Take an Observation and Return a Vector of the Likelihood over all Posterior Parameters
#Need to Add Arguments for Different Models

llfun <- function(i, data, draws){
  niter <- length(draws$lp__)
  #idata <- data[i,]
  idata <- data
  col_idx <- niter*(idata$model-1)  
 
  #Functions for log sev CDF and log sev PDF
  #Seems like We need this Within the Function for loo to recognize themm
  sev_logpdf <- function(y, mu, sigma){
    z = (y - mu) / sigma
    out= (-log(sigma) + z - exp(z))
    return(out)
  }
  
  sev_logcdf <- function(y, mu, sigma){
    out = (-exp((y - mu) / sigma))
    return(log(1-exp(out)))
  }
  
  sev_logccdf <- function(y, mu, sigma){
    return (-exp((y - mu) / sigma))
  }
  
  if(idata$censored==TRUE){
    likenum <- with(idata,
                  log(exp(draws$log_pi[1:niter + col_idx*pi.free] +
                  sev_logpdf(endtime, draws$mu1[1:niter + col_idx*mu1.free], draws$sigma1[1:niter + col_idx*sigma1.free]) +
                  sev_logccdf(endtime, draws$mu2[1:niter + col_idx*mu2.free], draws$sigma2[1:niter + col_idx*sigma2.free])) +
                  exp(sev_logpdf(endtime, draws$mu2[1:niter + col_idx*mu2.free], draws$sigma2[1:niter + col_idx*sigma2.free]) +
                  log(1 - exp(draws$log_pi[1:niter + col_idx*pi.free] + sev_logcdf(endtime, draws$mu1[1:niter + col_idx*mu1.free], draws$sigma1[1:niter + col_idx*sigma1.free])))))
                    )
    likedem <- with(idata, 
                    log(1-exp(draws$log_pi[1:niter + col_idx*pi.free] + sev_logcdf(starttime,draws$mu1[1:niter + col_idx*mu1.free], draws$sigma1[1:niter + col_idx*sigma1.free]))) + 
                      sev_logccdf(starttime, draws$mu2[1:niter + col_idx*mu2.free], draws$sigma2[1:niter + col_idx*sigma2.free])
                      )
  } 
  else{
    likenum <- with(idata, 
                    log(1-exp(draws$log_pi[1:niter + col_idx*pi.free] + sev_logcdf(endtime,draws$mu1[1:niter + col_idx*mu1.free], draws$sigma1[1:niter + col_idx*sigma1.free]))) + 
                    sev_logccdf(endtime,draws$mu2[1:niter + col_idx*mu2.free], draws$sigma2[1:niter + col_idx*sigma2.free])
                    )
    
    likedem <- with(idata, 
                    log(1-exp(draws$log_pi[1:niter + col_idx*pi.free] + sev_logcdf(starttime, draws$mu1[1:niter + col_idx*mu1.free], draws$sigma1[1:niter + col_idx*sigma1.free]))) + 
                    sev_logccdf(starttime,draws$mu2[1:niter + col_idx*mu2.free], draws$sigma2[1:niter + col_idx*sigma2.free])
                    )
    }
  return(likenum-likedem)
}




#Make Matrix for Loo Based On Bayes Model with Extra Flag Columns
data.model.flags <- function(d,pi.free=F, mu1.free=F, sigma1.free=F, mu2.free=F, sigma2.free=F){
  n <- nrow(d)
  d$pi.free <- rep(pi.free,n)
  d$mu1.free <- rep(mu1.free,n)
  d$sigma1.free <- rep(sigma1.free,n)
  d$mu2.free <- rep(mu2.free,n)
  d$sigma2.free <- rep(sigma2.free,n)
  return(d)
}


#I'm Sure There is a Way to Write This in a Loop or Function

samp = extract(s1)
d.augment <- data.model.flags(data_all, pi.free=F, mu1.free=F, sigma1.free = F, mu2.free = F, sigma2.free = F)
N <- nrow(d.augment)
S <- nrow(samp$lp__)
loo_output_1 <- loo(llfun, args = list(data=d.augment, N=N, S=S, draws=samp)) 

samp = extract(s2)
d.augment <- data.model.flags(data_all, pi.free=F, mu1.free=F, sigma1.free = F, mu2.free = T, sigma2.free = F)
N <- nrow(d.augment)
S <- nrow(samp$lp__)
loo_output_2 <- loo(llfun, args = list(data=d.augment, N=N, S=S, draws=samp)) 

samp = extract(s3)
d.augment <- data.model.flags(data_all, pi.free=F, mu1.free=F, sigma1.free = F, mu2.free = T, sigma2.free = T)
N <- nrow(d.augment)
S <- nrow(samp$lp__)
loo_output_3 <- loo(llfun, args = list(data=d.augment, N=N, S=S, draws=samp))

samp = extract(s4)
d.augment <- data.model.flags(data_all, pi.free=T, mu1.free=F, sigma1.free = F, mu2.free = T, sigma2.free = T)
N <- nrow(d.augment)
S <- nrow(samp$lp__)
loo_output_4 <- loo(llfun, args = list(data=d.augment, N=N, S=S, draws=samp))


#loo_matrix <- sapply(1:2000, function(i) llfun(i, d.augment[i,, drop=FALSE], samp))




