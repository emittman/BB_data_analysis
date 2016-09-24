d <- readRDS("unit_summaries_ALL.rds")

d$model <- as.factor(d$model)

## subset by:
#### at least 1 failure..
num_failed <- ddply(d, .(model), summarise, fails = sum(failed > 0))

model_names <- num_failed$model[num_failed$fails > 0]

id_1failure <- d$model %in% model_names


#### logically possible
id_possible <- d$end_time != d$start_time & d$end_time != 0

id_inc <- id_1failure & id_possible

#subset and refactor d!
d <- d[id_inc,]

d <- droplevels(d)

saveRDS(d, "clean_unit_summaries.rds")