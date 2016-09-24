library(RSQLite)
library(plyr)
library(dplyr)
name13 <- "C:/Users/emittman/GitHub/533/docs_2013/docs_2013/drive_stats.db"
name14 <- "C:/Users/emittman/GitHub/533/docs_2014/docs_2014/drive_stats.db"
name15 <- "C:/Users/emittman/GitHub/533/docs_2015/docs_2015/drive_stats.db"
name16 <- "C:/Users/emittman/GitHub/533/docs_Q1_2016/drive_stats.db"
#need to check how unit was withdrawn
q <- "SELECT serial_number, model, MIN(Date), MAX(Date), MIN(smart_9_raw), MAX(smart_9_raw), MAX(failure)   
      FROM drive_stats   
      GROUP BY serial_number;  
      "
subset <- NA
names <- c(name13, name14, name15,name16)
for(name in names){
  conn <- dbConnect(drv = SQLite(),
                    dbname = name)
  
  
  temp <- dbGetQuery(conn, q)
  
  
  subset <- rbind(subset, temp)
  
  dbDisconnect(conn)
}

#format dates
subset[,3] <- as.Date(subset[,3])
subset[,4] <- as.Date(subset[,4])
names(subset)[3:7] <- c("start_date","end_date","start_time","end_time", "failed")

stats <- ddply(subset, .(serial_number), summarize, 
               model = model[1],
               start_date = min(start_date),
               end_date = max(end_date),
               start_time = min(start_time),
               end_time = max(end_time),
               failed = sum(failed)
)

saveRDS(stats, file="unit_summaries_ALL.rds")

