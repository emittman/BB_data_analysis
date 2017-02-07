## Load In Samples ##

s <- readRDS("workflow/samples_2_1.rds")
samp <- extract(s)

## Get Quantiles from GFLP Model ##

#Weibull Function
my_pweibull <- function(x, location, scale, ...){
  pweibull(x, 1/scale, exp(location), ...)
}

#Function to get time from p
quant <- function(loc1, loc2, scl1, scl2, pi, b, t) {
  1 - (1 - pi * my_pweibull(t, loc1, scl1)) * (1 - my_pweibull(t, loc2, scl2)) - b
  }

#Function to Get Quantiles of Interest for All Posterior Draws.  Takes Model Argument and Quantile.
bfcn <- function(model, b){
    j <- model 
    Qp <- sapply(1:24000, function(i){
          uniroot(function(t) quant(samp$mu1[i,j],samp$mu2[i,j],samp$sigma1[i,j],samp$sigma2[i,j], exp(samp$log_pi[i,j]), b, t),c(0,300000))$root
  })
  q <- quantile(Qp, c(.025,.5,.975))
  data.frame(model=j,lower=q[1],med=q[2],upper=q[3])
}

#Run this Function for All Models.  Note:  Should be a quick way to combine this with first function.
df <- NULL
for (j in 1:21){
  df <- rbind(df,bfcn(j,.1))
  rownames(df) <- NULL
}


