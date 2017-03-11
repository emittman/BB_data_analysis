setwd("constrained_MLE")
library(rstan)
library(plyr)
library(dplyr)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
mixmle <- stan_model(file = "MLE_gflp.stan")

#Get Data Ready.  Note: Two Years is 17520 hours; 1 Year is 8760
#For overview, the definitions are: early is less than 1 year, late is greater than 2 years.

prepare_data <- function(lb_fails = 0, lb_late_fails = 0, lb_early_fails = 0, infant=0){
  #load data
  overview  <- readRDS("../BB_data/overview.rds")
  dat       <- readRDS("../BB_data/clean_unit_summaries.rds")
  #subsetting
  id        <- with(overview, which(f >= lb_fails & late_f >= lb_late_fails & early_f >= lb_early_fails))
  dat$model <- as.integer(dat$model)
  #dat$mode2<-ifelse(dat$end_time>wear, 1, 0)
  dat$mode1<-ifelse(dat$endtime<infant, 1, 0)
  df <- with(subset(dat, model %in% id),
             data.frame(endtime = endtime,
                        starttime = starttime,
                        censored = failed == 0,
                        mode1 = mode1 == 1,
                        model = as.integer(factor(model)))
  )
  # return(overview)
  return(df)
}


dat <- prepare_data(5,1,1,8760)  # 5 total failures, 1 early, 1 late, Less than 1 Year is Mode 1 (infant)

#Stan MLE Function for One Drive Brand Model

mle_stan <- function(brand,start){
  
  df <- subset(dat,model==brand)
  #Get Data Ready for Stan
  stan_dats <- with(df,
                  list(N_obs = sum(!(censored) & !(mode1)),
                       N_cens = sum((censored)),
                       N_obs1 = sum(!(censored) & (mode1)),
                       starttime_obs = log(starttime[!(censored) & !(mode1)]+1),
                       starttime_cens = log(starttime[(censored)]+1),
                       endtime_obs = log(endtime[!(censored) & !(mode1)]+1),
                       endtime_cens = log(endtime[(censored)]+1),
                       starttime_obs1 = as.array(log(starttime[!(censored) & (mode1)]+1)),
                       endtime_obs1 = as.array(log(endtime[!(censored) & (mode1)]+1)),
                       p = c(.5, .2)))
  
  o <- optimizing(object = mixmle, data = stan_dats,algorithm="BFGS", init = start, hessian=TRUE)
  return(o)
}
  

#Initial Values Based on Bayes Posteriors
inits <- list(log_tp1=7,
              log_tp2=9,
              log_sigma1=.25,
              log_sigma2=-.5,
              pi=.02)

#Get MLE from Stan and Standard Errors in Terms of mu1, mu2, sigma1, sigma2, pi
get_mle_stan <- function(model,initial){
  mle <- mle_stan(model,initial) #this step gets data, too.
  est <- mle$par
  cov <- solve(-mle$hessian) #est. covariance for (log_tp1, log_tp2,log_sigma1,log_sigma2,pi)
  sev1 <- log(-log(1-.5)) #estimating at .5 quantile for early failure
  sev2 <- log(-log(1-.2)) #estimating at .2 quantile for wearout
  D <- matrix(c(1,0,0,0,0,0,1,0,0,0,-exp(est[3])*sev1,0,exp(est[3]),0,0,0,-exp(est[4])*sev2,0,exp(est[4]),0,0,0,0,0,(1/(1+exp(-est[5])))*((exp(-est[5]))/(1+exp(-est[5])))),nrow=5,ncol=5)
  cov2 <- D %*% cov %*% t(D)   #est. covariance for (mu1,mu2, sigma1,sigma2,pi);
  out <- list(est=est, cov=cov, cov2=cov2) #note: point estimates are still in original scale
  attr(out,"model") <- model
  return(out)
}

mod2 <- get_mle_stan(19,inits)
  
#Return Parameters and CI for MLE, i diag of cov2 (mu1,mu2, sigma1,sigma2,pi).
bandsmle <- function(out,alpha, pm, i){
  lb = out$est[pm] - qnorm(ci)*out$cov2[i,i]
  mean = out$est[pm] 
  ub = out$est[pm] + qnorm(ci)*out$cov2[i,i]
  return(cbind(lb,mean,ub))
}

#Example, Get 95% CI for mu1 
bandsmle(test,.975,"mu1",1)  


#Write Function to Get CI for Quantile?

#Wald Bands for F(t); Section 8.4.3 Meeker
  
#SEV for standardized time
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

#Derivative of SEV
  dsev = function(x, location = 0, scale = 1){
    pdf <- exp(x-exp(x))
    return(pdf)
  }
  

    

#Delta Method Here to Get Wald Bands on the GFLP CDF
  WaldBand_mle <- function(mle, alpha, begin, end, length.out=100){
    b <- data.frame(time = seq(begin, end, length.out=length.out)) %>%
      ddply(.(time), summarise,
            z1     = (log(time) - mle$est["mu1"]) / exp(mle$est["log_sigma1"]),
            z2     = (log(time) - mle$est["mu2"]) / exp(mle$est["log_sigma2"]),
            g=1-((1-mle$est["pi"]*psev(z1,lower.tail = TRUE))*(psev(z2,lower.tail=FALSE))),
            gm1=-(mle$est["pi"])*dsev(z1)*(1/exp(mle$est["log_sigma1"]))+(mle$est["pi"]*dsev(z1)*psev(z2,lower.tail=TRUE)*(1/exp(mle$est["log_sigma1"]))),
            gm2=-dsev(z2)*(1/exp(mle$est["log_sigma2"]))+((mle$est["pi"])*psev(z1,lower.tail=TRUE)*dsev(z2)*(1/exp(mle$est["log_sigma2"]))),
            gs1=-(mle$est["pi"])*psev(z2,lower.tail = TRUE)*dsev(z1)*((-log(time)+mle$est["mu1"])/((exp(mle$est["log_sigma1"]))^2))+(mle$est["pi"]*dsev(z1)*((-log(time)+mle$est["mu1"])/((exp(mle$est["log_sigma1"]))^2))),
            gs2=dsev(z2)*((-log(time)+mle$est["mu2"])/((exp(mle$est["log_sigma2"]))^2))-(mle$est["pi"]*psev(z1,lower.tail=TRUE)*dsev(z2)*((-log(time)+mle$est["mu2"])/((exp(mle$est["log_sigma2"]))^2))),
            gp=psev(z1)-(psev(z1)*psev(z2)),
            see=(as.matrix(cbind(gm1,gm2,gs1,gs2,gp)) %*% mle$cov2 %*% t((as.matrix(cbind(gm1,gm2,gs1,gs2,gp))))) ^ 0.5,
            dF=(mle$est["pi"]*dsev(z1)*(1/exp(mle$est["log_sigma1"])))*(psev(z2,lower.tail=FALSE))+(dsev(z2)*(1/exp(mle$est["log_sigma2"]))*(1-mle$est["pi"]*psev(z1,lower.tail = TRUE))))  #Pg 169 (Hong & Meeker) df/dlog(t)= 1/sigma1*f1*(1-F2)*p+f2*(1/sigma2)*(1-p*F1)
    
    
    b$lower = with(b, g - qnorm(1-alpha/2) * see)
    b$upper = with(b, g + qnorm(1-alpha/2) * see)
    b$yl=with(b,exp(log(time)-((qnorm(1-alpha/2) * see)/dF)))
    b$yu=with(b,exp(log(time)+((qnorm(1-alpha/2) * see)/dF)))
    
    b
  }

#Meeker Bands that Respect Parameter Space
  
#Lets make Function for Standardized GFLP CDF: F-hat(exp(yl),Fhat(exp(yu)))
gfp<-function(time,mle){
  z1     = (log(time) - mle$est["mu1"]) / exp(mle$est["log_sigma1"])
  z2     = (log(time) - mle$est["mu2"]) / exp(mle$est["log_sigma2"])
  g=1-((1-mle$est["pi"]*psev(z1,lower.tail = TRUE))*(psev(z2,lower.tail=FALSE)))
  return(g)
}
  
x1=gfp((check$yl),db6)
x2=gfp((check$yu),db6)
                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                          

#Test on Model 6
df6=dat[dat$model==6,]  
check=WaldBand_mle(db6,alpha=.05,begin = 1000, end = 50000)
plot(check$time,check$g,ylab="CDF",type="l")
lines(check$time,check$lower,lty=3)
lines(check$time,check$upper,lty=3)



###Ignore Below ###  
  
  
#Check Analytical Derivatives to Numerical Ones#
library(numDeriv)
gfp<-function(time,mu1,mu2,p,sig1,sig2){
  z1     = (log(time) - mu1) / sig1
  z2     = (log(time) - mu2) / sig2
  g=1-((1-p*psev(z1,lower.tail = TRUE))*(psev(z2,lower.tail=FALSE)))
  return(g)
}

mu1=6
mu2=12
sig1=1
sig2=2
p=.06
time=5


f <- function(x) {
1-((1-p*psev((x - mu1) / sig1,lower.tail = TRUE))*(psev((x - mu2) / sig2,lower.tail=FALSE)))
}

grad(f,c(6,12,1,2,.06,1.609438))
check1=-.06*dsev((log(time)-mu1)/sig1)*(1/sig1)+p*dsev((log(time)-mu1)/sig1)*psev(log(time)-mu2/sig2)*(1/sig1)
check2=-dsev((log(time)-mu2)/sig2)*(1/sig2)+p*psev((log(time)-mu1)/sig1)*dsev(log(time)-mu2/sig2)*(1/sig2)
check3=-p*psev(log(time)-mu2/sig2,lower.tail = TRUE)*dsev(log(time)-mu1/sig1)*((-log(time)+mu1)/((sig1)^2))+(p*dsev(log(time)-mu1/sig1)*((-log(time)+mu1)/(sig1^2)))
check4=dsev(z2)*((-log(time)+mu2)/(sig2^2))-(p*psev(z1,lower.tail=TRUE)*dsev(z2)*((-log(time)+mu2)/((sig2)^2)))
check5=psev(log(time)-mu1/sig1)-(psev(log(time)-mu1/sig1)*psev(log(time)-mu2/sig2))
check6=p*dsev(z1)*(1-psev(z2))*(1/sig1)+(dsev(z2)*(1-p*psev(z1))*(1/sig2))



#Compare MLEs to Bayes Point Estimtes for 10/13
mles <- readRDS("mle_infant1year.rds")

#Extract Bayes Parameters
s <- readRDS("../MCMC_draws/samples_logodds_reduced_2_26.rds")
samp <- extract(s)

#Get Parameters for a Specific Model
summary(s)$summary[c("mu1[1]","mu2[1]","sigma1[1]", "sigma2[1]", "log_pi[1]"),]

bayesout <- matrix(ncol=5, nrow=22)
for (i in 1:22){
    num=i
    spi <- paste0("pi[",num,"]",collapse="")
    smu1 <- paste0("mu1[",num,"]",collapse="")
    smu2 <- paste0("mu2[",num,"]",collapse="")
    slogsig1 <- paste0("log_sigma1[",num,"]",collapse="")
    slogsig2 <- paste0("log_sigma2[",num,"]",collapse="")
    pi <- summary(s)$summary[spi,"50%"]
    mu1 <- summary(s)$summary[smu1,"50%"]
    mu2 <- summary(s)$summary[smu2,"50%"]
    logsig1 <- summary(s)$summary[slogsig1,"50%"]
    logsig2 <- summary(s)$summary[slogsig2,"50%"]
    

  bayesout[i,1]<-mu1
  bayesout[i,2]<-mu2 
  bayesout[i,3]<-logsig1
  bayesout[i,4]<-logsig2
  bayesout[i,5]<-pi
}

bayesout.dat<-as.data.frame(bayesout)
colnames(bayesout.dat) <- c("mu1_b","mu2_b","log_sigma1_b","log_sigma2_b","pi_b")
colnames(out)<- c("mu1_f","mu2_f","log_sigma1_f","log_sigma2_f","pi_f","LL","Grid")
out<-as.data.frame(out)

#Combine Bayes and MLE estimates
combo_parm <- data.frame("mu1_b" = bayesout.dat$mu1_b, "mu_1f" = out$mu1_f,
                         "mu2_b"=bayesout.dat$mu2_b,"mu2_f"=out$mu2_f,
                         "log_sigma1_b"=bayesout.dat$log_sigma1_b,"log_sigma1_f"=out$log_sigma1_f,
                         "log_sigma2_b"=bayesout.dat$log_sigma2_b,"log_sigma2_f"=out$log_sigma2_f,
                         "pi_b"=bayesout.dat$pi_b,"pi_f"=out$pi_f)

saveRDS(combo_parm,"bayesmle.rds")



