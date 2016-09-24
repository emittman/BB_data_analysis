dat <- readRDS("clean_unit_summaries.rds")
dat$model <- as.integer(dat$model)

library(plyr)
library(dplyr)
overview <- ddply(dat, .(model), summarise,
                  n=length(model),
                  f=sum(failed>0),
                  early_f = sum(failed>0 & end_time<365*24*1),
                  late_f = sum(failed>0 & end_time>365*24*2))

saveRDS(overview, file="overview.rds")
