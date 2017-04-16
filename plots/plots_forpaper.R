library(backblaze)
library(plyr)
library(dplyr)
library(rstan)

source("../workflow/functions.R")
source("../plotting_fns/KM_plot.R")
source("../plotting_fns/greenwood_errors.R")

dat <- backblaze


dat$model <- as.integer(dat$model)
dat$failed[dat$failed==2] <- 1

#Look at Distribution of Failures and Total Time on Test
#Total Filures By Model
tot<- dat %>% select(model, censored) %>% group_by(model) %>% summarise(sum(censored==FALSE))
names(tot)[2]<-"count"
theme_set(theme_bw(base_size=18))
ggplot(tot, aes(x = count)) + geom_histogram() + scale_x_log10()
ggplot(tot, aes(x=factor(model),y=count))+geom_bar(stat="identity")+xlab("Model")+ylab("Total Failures")+scale_x_discrete(breaks=seq(1, 62, by=2))+
  ggtitle("Total Failures by Model") + scale_y_log10()

#Total Time on Test
tt<-all %>% group_by(model) %>% mutate(tt=(endtime-starttime)/(24*365))
tt2<-tt %>% group_by(model) %>% summarize(sum(tt))
names(tt2)[2]<-"total"
theme_set(theme_bw(base_size=18))
ggplot(tt2, aes(x=factor(model),y=total))+geom_bar(stat="identity")+xlab("Model")+ylab("Total Time (Years)")+scale_x_discrete(breaks=seq(1, 62, by=2))+ggtitle("Total Time on Test by Model")+scale_y_log10()