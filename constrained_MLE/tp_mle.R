##  Hong Meeker Plug in Method for tp CI, Page 172 Hong & Meeker ##

# First Get the Value of the Quantile for Fixed MLE estimated Parameters ##

#Weibull Function
my_pweibull <- function(x, location, scale, ...){
  pweibull(x, 1/scale, exp(location), ...)
}

#Function to get time from p
quant <- function(loc1, loc2, scl1, scl2, pi, b, t) {
  1 - (1 - pi * my_pweibull(t, loc1, scl1)) * (1 - my_pweibull(t, loc2, scl2)) - b
}

#Function to Get Quantiles of Interest for MLE. Takes MLE output and Quantile.
mletp <- function(mle, b){
  Qp = uniroot(function(t) quant(mle$est["mu1"],mle$est["mu2"],exp(mle$est["log_sigma1"]), exp(mle$est["log_sigma2"]), mle$est["pi"], b, t),c(0,300000))$root
  return(Qp)
}


# SEV for standardized time
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

#Get Implicit Derivative of GFLP
impli_deriv <- function(mle, tp){
    z1   = (log(tp) - mle$est["mu1"]) / exp(mle$est["log_sigma1"])
    z2   = (log(tp) - mle$est["mu2"]) / exp(mle$est["log_sigma2"])
    gF1  = -psev(z2,lower.tail=FALSE)*mle$est["pi"]
    gF2  = -(1 - mle$est["pi"]*psev(z1,lower.tail=TRUE))
    F1yp = dsev(z1)*(1/exp(mle$est["log_sigma1"]))
    F2yp = dsev(z2)*(1/exp(mle$est["log_sigma2"]))
    F1m1 = -dsev(z1)*(1/exp(mle$est["log_sigma1"]))
    F1s1 = dsev(z1)*((-log(tp) + mle$est["mu1"])/((exp(mle$est["log_sigma1"]))^2))
    F2m2 = -dsev(z2)*(1/exp(mle$est["log_sigma2"]))
    F2s2 = dsev(z2)*((-log(tp) + mle$est["mu2"])/((exp(mle$est["log_sigma2"]))^2))
    denom = gF1*F1yp + gF2*F2yp
    l1 = -gF1*F1m1/denom
    l2 = -gF2*F2m2/denom
    ls1 = -gF1*F1s1/denom
    ls2 = -gF2*F2s2/denom
    lp = 0
    out <- cbind(l1,l2,ls1,ls2,lp)
    rownames(out) <- NULL
    return(out)
}

#Get Standard Error for y=log(tp)
se_tp <- function(impl, mle){
  out <- impl%*%mle$cov2%*%t(impl)
  return(sqrt(out))
}

#Get CI for tp using Delta Plug in Method
tp_bounds <- function(mle, tp, se, alpha){
  lb <- exp(log(tp) - qnorm(1-alpha/2)*se)
  med <- tp
  ub  <- exp(log(tp) + qnorm(1-alpha/2)*se)
  out <- cbind(lb,med,ub)
  return(out)
}

## Test the Series of Functions for Model 2 ##


mod2 <- get_mle_stan(2,inits)  #get MLE estimate
tp_mod2 <- mletp(mod2,.10)     #use MLE point estimates for tp
idev_mod2 <- impli_deriv(mod2, tp_mod2) #get implicit derivative of GFLP at log(tp)
se_mod2 <- se_tp(idev_mod2,mod2) #se for log(tp)
mod2_bounds <- tp_bounds(mod2, tp_mod2, se_mod2, .05) #get 95% tp CI
mod2_bounds


