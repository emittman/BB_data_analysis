## Load In Samples ##

#s <- readRDS("samples_logodds_reduced_2_26.rds")
s <- readRDS("samples_2_1.rds")
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
bfcn <- function(model, b, samp, quantiles = c(.025, .5, .975)){
    j <- model 
    Qp <- sapply(1:dim(samp$lp__), function(i){
          uniroot(function(t) quant(samp$mu1[i,j],samp$mu2[i,j],samp$sigma1[i,j],samp$sigma2[i,j], exp(samp$log_pi[i,j]), b, t),c(0,300000))$root
  })
  q <- quantile(Qp, quantiles)
  df <- data.frame(model=j,t(q))
  rownames(df) <- NULL
  names(df)[-1] <- sapply(quantiles, function(q) as.character(round(q,3)))
  df
}

#reduced version
bfcn2 <- function(model, b, samp, quantiles = c(.025, .5, .975)){
  
  j <- model 
  Qp <- sapply(1:dim(samp$lp__), function(i){
    uniroot(function(t) quant(samp$mu1[i],samp$mu2[i,j],samp$sigma1[i],samp$sigma2[i,j], exp(samp$log_pi[i,j]), b, t),c(0,300000))$root
  })
  q <- quantile(Qp, quantiles)
  df <- data.frame(model=j,t(q))
  rownames(df) <- NULL
  names(df)[-1] <- sapply(quantiles, function(q) as.character(round(q,3)))
  df
}

#Run this Function for All Models.  Note:  Should be a quick way to combine this with first function.
# df_reduced <- NULL
# for (j in 1:21){
#   df_reduced <- rbind(df_reduced,bfcn2(j,.1))
# }

df <- NULL
for (j in 1:21){
  df <- rbind(df,bfcn(j,.1, samp))
}

xtable(cbind(df, df_reduced))