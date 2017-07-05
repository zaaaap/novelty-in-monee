library(readr)
library(data.table)
library(dplyr)

setwd("~/git/monee-novelty/Rstats")

sensor_vals <- fread("~/git/monee-novelty/Rstats/output2.txt")
sensor_vals[,V25:=NULL]
sensor_vals <- unique(sensor_vals)

samples <- sample_n(sensor_vals, 300)

write.table(samples, file = "sample_sensors-300.csv", col.names=FALSE, row.names=FALSE)