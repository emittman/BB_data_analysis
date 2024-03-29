library(plyr)
library(dplyr)
library(rstan)
library(cowplot)

dat <- readRDS("../BB_data/clean_unit_summaries.rds")
dat$model <- as.integer(dat$model)
dat$failed[dat$failed==2] <- 1

overview <- ddply(dat, .(model), summarise,
                  n=length(model),
                  f=sum(failed>0),
                  early_f = sum(failed>0 & endtime<365*24*1),
                  late_f = sum(failed>0 & endtime>365*24*2))
# id <- unique(dat$model)
id <- with(overview, which(overview$early >= 0 & overview$late_f >= 0 & f >=3))
overview$stan_id <- NA
overview[id,]$stan_id <- 1:length(id)

#new id



#Look at Distribution of Failures and Total Time on Test
#Total Filures By Model

#Subset based on Model with >=3 failures
datsub = subset(dat, model %in% id)

tot <- datsub %>% select(model, censored) %>% group_by(model) %>% summarise(sum(censored==FALSE))
names(tot)[2]<-"count"

#Sort by Failures
tot <- as.data.frame(tot)
tot$model <- 1:length(id)  #Make New ID based on Subset
tot$model <- factor(tot$model,levels=tot$model[order(tot$count)])  #Sort by Total Failures


theme_set(theme_bw(base_size=18))
#ggplot(tot, aes(x = count)) + geom_histogram(bins = 15) + scale_x_log10()
fail <- ggplot(tot, aes(x=model, y=count)) + geom_bar(stat="identity") + xlab("Drive Brand Model") + ylab("Total Failures") + 
   ggtitle("Total Failures") + scale_y_log10(breaks=c(10,100,1000)) + coord_flip() 

#Total Time on Test
tt <- datsub %>% group_by(model) %>% mutate(tt=(endtime-starttime)/(24*365))
tt2 <- tt %>% group_by(model) %>% summarize(sum(tt))
names(tt2)[2]<-"total"

tt2 <- as.data.frame(tt2)
tt2$count <- tot$count
tt2$model <- 1:length(id)
tt2$model <- factor(tt2$model,levels=tt2$model[order(tt2$count)])  #Sort by Total Failures

theme_set(theme_bw(base_size=18))
testtime <- ggplot(tt2, aes(x=factor(model),y=total)) + geom_bar(stat="identity") + xlab("Drive Brand Model") + 
  ylab("Total Time (Years)") + ggtitle("Total Observed Running Time") + scale_y_log10(breaks=c(10,100,1000, 10000)) +
  coord_flip() + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

#Make Cow Plot
plot_grid(fail, testtime, ncol = 2, nrow = 1)



