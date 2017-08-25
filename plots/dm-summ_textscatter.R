#dm summaries

library(ggplot2)

dat <- readRDS("BB_data/clean_unit_summaries.rds")
ov <- readRDS("BB_data/overview.rds")

id <- which(ov$f>=3)
ov$stan_id <- NA
ov$stan_id[id] <- 1:44

dat <- dplyr::filter(dat, model %in% ov$model[id])

# model now matches stan group ids
dat$model <- ov$stan_id[match(dat$model, ov$model)]

library(plyr)
my_summary <- ddply(dat, .(model), function(x){
  tt <- sum(x$endtime - x$starttime)
  fs <- sum(x$censored == 0)
  data.frame(tt=tt, fs=fs)
})

my_summary$tt100k <- my_summary$tt/100000

#in yrs
#xbrks = c(10, 100, 1000, 10000)
xbrks2 = c(1, 10, 100, 1000)


setwd("paper/fig/")
pdf("dm-summ-scatter.pdf", width=7, height=5)
ggplot(my_summary, aes(x=tt100k, y=fs)) + geom_text(aes(label=model), size=4) + 
  scale_x_continuous(name="Total time observed (hundred thousands of hours)", trans="log", breaks = xbrks2, labels = xbrks2) +
  scale_y_continuous(name="Observed failures", trans="log", breaks = c(1, 10, 100, 1000)) +
  theme_bw(base_size=14) #+geom_smooth(method = "lm", se = F, linetype="dashed")
dev.off()
