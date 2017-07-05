#!/usr/bin/env Rscript
library(readr)
library(data.table)
library(dplyr)

setwd("~/git/monee-novelty/stats/behaviour-sampling")

sensor_vals <- fread("~/git/monee-novelty/stats/behaviour-sampling/output.txt")
sensor_vals[,V25:=NULL]
sensor_vals <- unique(sensor_vals)

samples <- sample_n(sensor_vals, 20)

write.table(samples, file = "sample_sensors-20.csv", col.names=FALSE, row.names=FALSE)








