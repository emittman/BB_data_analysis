setwd("review")

ov <- readRDS("../BB_data/overview.rds")
full_id <- with(ov, which(f>=3))
up_id <- with(ov, which(f>=5 & early_f>=1 & late_f>=1))
ov$full_id <- NA
ov$up_id <- NA
ov$full_id[full_id] <- 1:length(full_id)
ov$up_id[up_id] <- 1:length(up_id)

s3 <- readRDS("../workflow2/samples_review_dm3.rds") #6
s8 <- readRDS("../workflow2/samples_review_dm8.rds") #14
s10 <- readRDS("../workflow2/samples_review_dm10.rds") #16
s14 <- readRDS("../workflow2/samples_review_dm14.rds") #23
s18 <- readRDS("../workflow2/samples_review_dm18.rds") #31


full <- readRDS("../workflow/samples_lor_only3fails.rds")
samp_full <- extract(full)

psev = function(x, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE){
  z = (x - location)/scale
  upper <- exp(-exp(z))
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
solve4B10 <- Vectorize(function(log_pi, mu1, sigma1, mu2, sigma2){
  #solve .1 = 1 - (1 - pi * pweibull(x, mu1, sigma1))(1 - pweibull(x, mu2, sigma2))
  f <- function(x) {
      .9 - (1 - exp(log_pi)*psev(x, location=mu1, scale=sigma1)) *
                        (1 - psev(x, location=mu2, scale=sigma2))
    }
  uniroot(f, c(-100, 1000))$root
})

solve4B10(-4.4, 8.5, 1.2, 12.2, .4)

unpooled_df <- ldply(c("s3","s8","s10","s14","s18"), function(i){
  samp <- extract(get(i))
  B10 <- with(samp, solve4B10(log_pi, mu1, sigma1, mu2, sigma2))
  B10_interval <- quantile(B10, c(.05,.5,.95))
  data.frame(set = i, lower_B10 = exp(B10_interval[1]), median_B10 = exp(B10_interval[2]),
             upper_B10 = exp(B10_interval[3]))
})
unpooled_df$set <- c(6, 14, 16, 23, 31)
pooled_df <- ldply(c(6, 14, 16, 23, 31), function(i){
  B10 <- with(samp_full, solve4B10(log_pi[,i], mu1, sigma1, mu2[,i], sigma2[,i]))
  B10_interval <- quantile(B10, c(.05,.5,.95))
  data.frame(set = i, lower_B10 = exp(B10_interval[1]), median_B10 = exp(B10_interval[2]),
             upper_B10 = exp(B10_interval[3]))
})

unpooled_df$n <- ov$n[which(ov$full_id %in% c(6, 14, 16, 23, 31))]
unpooled_df$model <- "unpooled"
pooled_df$n <- ov$n[which(ov$full_id %in% c(6, 14, 16, 23, 31))]
pooled_df$model <- "partially pooled"
comb <- rbind(unpooled_df, pooled_df)
comb$set <- factor(comb$set)

library(ggplot2)
ggplot(comb, aes(x=set, y=median_B10, ymin=lower_B10, ymax=upper_B10, color=model))+
  geom_pointrange(position=position_dodge(.5)) + theme_bw()
