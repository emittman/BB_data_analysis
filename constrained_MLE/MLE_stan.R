setwd("constrained_MLE")
library(rstan)
library(plyr)
library(dplyr)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
mixmle <- stan_model(file = "MLE_gflp.stan")

#Get Data Ready.  Note: Two Years is 17520 hours; 1 Year is 8760
#For OverView early is less than 1 year, late is greater than 2 years.

prepare_data <- function(lb_fails = 0, lb_late_fails = 0, lb_early_fails = 0,infant=0){
  #load data
  overview  <- readRDS("../BB_data/overview.rds")
  dat       <- readRDS("../BB_data/clean_unit_summaries.rds")
  #subsetting
  id        <- with(overview, which(f >= lb_fails & late_f >= lb_late_fails & early_f >= lb_early_fails))
  dat$model <- as.integer(dat$model)
  #dat$mode2<-ifelse(dat$end_time>wear, 1, 0)
  dat$mode1<-ifelse(dat$end_time<infant, 1, 0)
  df <- with(subset(dat, model %in% id),
             data.frame(endtime = end_time,
                        starttime = start_time,
                        censored = failed == 0,
                        mode1 = mode1 == 1,
                        model = as.integer(factor(model)))
  )
  # return(overview)
  return(df)
}

#8 total failures, 1 early, and < 1 year is infant. 
dat<-prepare_data(8,1,0,8760)

#Stan MLE Function for One Drive Brand Model

mle_stan<-function(brand,start){
  
  df<-subset(dat,model==brand)
  
  #Get Data Ready for Stan
  stan_dats <- with(df,
                  list(N_obs = sum(!(censored) & !(mode1)),
                       N_cens = sum((censored)),
                       N_obs1 = sum(!(censored) & (mode1)),
                       starttime_obs = log(starttime[!(censored) & !(mode1)]+1),
                       starttime_cens = log(starttime[(censored)]+1),
                       endtime_obs = log(endtime[!(censored) & !(mode1)]+1),
                       endtime_cens = log(endtime[(censored)]+1),
                       starttime_obs1 = log(starttime[!(censored) & (mode1)]+1),
                       endtime_obs1 = log(endtime[!(censored) & (mode1)]+1),
                       p = c(.5, .1)))
  
  o <- optimizing(object = mixmle, data = stan_dats,algorithm="BFGS", init = start, hessian=TRUE)
  return(o)
  
}
  
inits <- list(log_tp1=4,
              log_tp2=9,
              log_sigma1=.25,
              log_sigma2=-.5,
              pi=.01)


#All Models Seem OK except for 18, which has only 1 early failure
out<-matrix(nrow=22,ncol=7)
colnames(out)<-c('log_tp1','log_tp2','log_sigma1','log_sigma2','pi','mu1','mu2')
for (i in 1:22){
test=mle_stan(i,inits)
out[i,]<-test$par
}

#Grid Search: 243 combinations
grid<-expand.grid(log_tp1=c(2,4,6), log_tp2=c(9,11,13),
                  log_sigma1=c(.05,.15,.25),log_sigma2=c(-.8,-.6,-.2),pi=c(.01,.05,.10)) 

#Store Results
mlegrid<-matrix(nrow=243,ncol=6) #Matrix
colnames(mlegrid)<-c('log_tp1','log_tp2','log_sigma1','log_sigma2','pi','likelihood')

#Run Search for Model 8, which is pretty stable
  m<-8
  for (i in 1:243){
  fit<-mle_stan(m,grid[i,])
  mlegrid[i,1:5]<-fit$par[1:5]
  mlegrid[i,6]<-fit$value
}


#Run Grid Search on Models with Lots of Failures: 6,7,13,20,29,49,
grid_out<-list()
subfail<-c(6,7,13,20,29,49)
for (i in subfail){
  est<-grid_mle(i)
  grid_out[i]<-est
}

m6_grid<-grid_mle(6)
m7_grid<-grid_mle(7)
m13_grid<-grid_mle(13)
m29_grid<-grid_mle(29)

#Model 20 is pretty stable; let's get MLE and Standard Errors for mu and sigma.
inits <- list(log_tp1=4,
              log_tp2=10,
              log_sigma1=.25,
              log_sigma2=-.5,
              pi=.1
)

  get_mle_stan<-function(model,initial){
  mle<-mle_stan(model,initial) #this step gets data, too.
  est <- mle$par
  cov <- solve(-mle$hessian) #est. covariance for (log_tp1, log_tp2,log_sigma1,log_sigma2,pi)
  sev1<-log(-log(1-.5)) #estimating at .5 quantile for early failure
  sev2<-log(-log(1-.2)) #estimating at .2 quantile for wearout
  D<-matrix(c(1,0,0,0,0,0,1,0,0,0,-exp(est[3])*sev1,0,exp(est[3]),0,0,0,-exp(est[4])*sev2,0,exp(est[4]),0,0,0,0,0,1),nrow=5,ncol=5)
  cov2 <- D %*% cov %*% t(D)   #est. covariance for (mu1,mu2, sigma1,sigma2,pi);
  out <- list(est=est, cov=cov, cov2=cov2) #note: point estimates are still in original scale
  attr(out,"model") <- model
  return(out)
  }
  
#Get MLE for Model 20  
db20<-get_mle_stan(29,inits)
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
  

    

#Delta Method Here to Get Wald Bands on the CDF
  WaldBand_mle <- function(mle, alpha, begin, end, length.out=100){
    b <- data.frame(time = seq(begin, end, length.out=length.out)) %>%
      ddply(.(time), summarise,
            z1     = (log(time) - mle$est["mu1"]) / exp(mle$est["log_sigma1"]),
            z2     = (log(time) - mle$est["mu2"]) / exp(mle$est["log_sigma2"]),
            g=1-((1-mle$est["pi"]*psev(z1,lower.tail = TRUE))*(psev(z2,lower.tail=FALSE))),
            gm1=(mle$est["pi"])*psev(z2,lower.tail = FALSE)*dsev(z1)*(1/exp(mle$est["log_sigma1"])),
            gm2=(1-(mle$est["pi"]*psev(z1,lower.tail=TRUE)))*dsev(z2)*(-1/exp(mle$est["log_sigma2"])),
            gs1=(mle$est["pi"])*psev(z2,lower.tail = FALSE)*dsev(z1)*((log(time)-mle$est["mu1"])/((exp(mle$est["log_sigma1"]))^2)),
            gs2=(1-mle$est["pi"]*psev(z1,lower.tail = TRUE))*dsev(z2)*((log(time)-mle$est["mu2"])/((exp(mle$est["log_sigma2"]))^2)),
            gp=-psev(z1)+(psev(z1)*psev(z2)),
            se_g=(as.matrix(cbind(gm1,gm2,gs1,gs2,gp)) %*% mle$cov2 %*% t((as.matrix(cbind(gm1,gm2,gs1,gs2,gp))))) ^ 0.5)
            #dF<-(mle$est["pi]*dsev(z1)*(1/exp(mle$est["log_sigma1"])))-(mle$est["pi]*dsev(z1)*psev(z2)*(1/exp(mle$est["log_sigma1"])))+dsev(z2)*(1/exp(mle$est["log_sigma2"]))-mle$est["pi]*dsev(z2)*(1/exp(mle$est["log_sigma2"]))*psev(z1)
    
    
    b$lower = with(b, g - qnorm(1-alpha/2) * se_g)
    b$upper = with(b, g + qnorm(1-alpha/2) * se_g)
    #b$yl=with(b,log(time)-(qnorm(1-alpha/2) * se_g)/dF)
    #b$yu=with(b,log(time)+(qnorm(1-alpha/2) * se_g)/dF)
    
    b
  }

#Meeker Bands that Respect Parameter Space

                                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                          

#Test on Model 20
df20=dat[dat$model==20,]  
check=WaldBand_mle(db20,alpha=.05,begin = min(df20$end_time[df20$failed==1])*.9, end = max(df20$end_time[df20$failed==1])*1.1)
plot(check$g,ylab="CDF")
lines(check$lower)
lines(check$upper)

                                                                                                           ###Ignore Below ###  
  
  
#Get MLE and SE for p; looks like SE is not on constrained space, but p is.
p<-NULL
db<-NULL
se<-NULL
for (j in id){
  fit<-mle_stan(j,inits)
  er=tryCatch(sqrt(diag(solve(-fit$hessian))),error=function(e) NA)
  p<-append(p,fit$par[13])
  db<-append(db,j)
  se<-append(se,er[13])
}  
d = data.frame(db, p,se)

#Tranform Interval with Inv Logit
ivlog<-function(x){
  return(exp(x)/(1+exp(x)))
}

#Get 95% CI 
pu=sapply(d[,2],function(x) log(x/(1-x)))
lb=pu-1.96*d$se
ub=pu+1.96*d$se
d$lb<-ivlog(lb)
d$ub<-ivlog(ub)

#Extract log sigma 1 and 2
#Get MLE and SE for p; looks like SE is not on constrained space, but p is.
betamle<-NULL
for (j in id){
  fit<-mle_stan(j,inits)
  er=tryCatch(sqrt(diag(solve(-fit$hessian))),error=function(e) NA)
  b1<-1/exp(fit$par[16])
  b2<-1/exp(fit$par[17])
  se1=er[11]
  se2=er[12]
  betamle <- rbind(betamle,data.frame(model=j, b1=b1,b2=b2,se1=se1,se2=se2))
  rownames(betamle) <- NULL
}  



