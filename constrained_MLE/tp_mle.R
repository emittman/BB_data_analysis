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
    F1mu1 = -dsev(z1)*(1/exp(mle$est["log_sigma1"]))
    F1s1 = dsev(z1)*((-log(tp) + mle$est["mu1"])/((exp(mle$est["log_sigma1"]))^2))
    F2mu2 = -dsev(z2)*(1/exp(mle$est["log_sigma2"]))
    F2s2 = dsev(z2)*((-log(tp) + mle$est["mu2"])/((exp(mle$est["log_sigma2"]))^2))
    denom = gF1*F1yp + gF2*F2yp
    l1 = (-gF1*F1mu1)/denom
    l2 = (-gF2*F2mu2)/denom
    ls1 = (-gF1*F1s1)/denom
    ls2 = (-gF2*F2s2)/denom
    lp = 0 #note: p is not a parameter in the individual failure distributions
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
  mle.est <- tp
  ub  <- exp(log(tp) + qnorm(1-alpha/2)*se)
  out <- cbind(lb, mle.est, ub)
  return(out)
}

## Run for All Models ##

tp_all <- matrix(NA, ncol=4, nrow=21)

for (i in 1:21){
  mod <- i
  mlefit <- get_mle_stan(mod, inits)  #get MLE estimate
  tp_mle <- mletp(mlefit,.10)     #use MLE point estimates for tp
  idev <- impli_deriv(mlefit, tp_mle) #get implicit derivative of GFLP at log(tp)
  se_ylog <- se_tp(idev, mlefit) #se for log(tp)
  tp.bounds <- tp_bounds(mlefit, tp_mle, se_ylog, .05) #get 95% tp CI
  tp_all[i,2:4] <- tp.bounds
  tp_all[i, 1] <- mod
}




