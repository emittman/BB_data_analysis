library(plyr)
library(dplyr)
library(bshazard)

dat <- readRDS("../BB_data/clean_unit_summaries.rds")

overview <- ddply(dat, .(model), summarise,
                  n=length(model),
                  f=sum(failed>0),
                  early_f = sum(failed>0 & endtime<365*24*1),
                  late_f = sum(failed>0 & endtime>365*24*2))

dat$model <- as.integer(dat$model)
id <- with(overview, which(overview$early >= 5 & overview$late_f >= 10))
overview$stan_id <- NA
overview[id,]$stan_id <- 1:length(id)

fits <- list()

for(i in 1:length(id)){
Survtmp <- with(filter(dat, model==id[i]), Surv(starttime, endtime, pmin(failed,1)))

fits[[i]] <- bshazard(Survtmp~1, data = dat)
}
length(fits)
pdf("hazards.pdf")
opar <- par(mfrow=c(3,3))
  for(i in 1:length(id)){
    plot(fits[[i]], main=paste(c("Drive model %d", id[i]), collapse=""))
  }
dev.off()

#all data
Survall <- with(filter(dat, model == 7), Surv(starttime, endtime, pmin(failed, 1)))
fitall <- bshazard(Survall~1, data=dat, phi=7)
plot(fitall)

#dummy data

my_rweibull <- function(n, mu, sigma){
  rweibull(n, shape=1/sigma, scale=exp(mu))
}

rgflp <- function(n, p, mu1, sigma1, mu2, sigma2, censortime=NULL){
  def <- rbinom(n, size=1, prob=p)
  out <- ifelse(!def, my_rweibull(n, mu2, sigma2),
                pmin(my_rweibull(n, mu1, sigma1), my_rweibull(n, mu2, sigma2)))
  cens <- 0
  if(!is.null(censortime)){
    cens <- out>censortime
    out <- pmin(censortime, out)
  }
  round(out)
  return(data.frame(endtime=out, cens=cens))
}

par(mfrow=c(3,3))
for(i in 1:9){
  dummy <- rgflp(5000, .05, 8, 2, 11, .8, 60000)
  dummy$died <- dummy$cens != 1
  dummy$start <- 0
  
  fit <- bshazard(Surv(dummy$start, dummy$endtime, dummy$died)~1, data=dummy)
  # with(fit, plot(time, log(hazard)))
  plot(fit)
}
  ?Surv

?bshazard
datsmall <- dat[1:10,]
sur <- with(datsmall, Surv(endtime))
plot(bshazard(endtime~1, datsmall))
