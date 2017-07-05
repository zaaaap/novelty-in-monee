#!/usr/bin/env Rscript
library(readr)
library(data.table)
library(dplyr)

setwd("~/git/monee-novelty/stats/behaviour_sampling")

sensor_vals <- fread("~/git/monee-novelty/stats/behaviour_sampling/output.txt")
sensor_vals[,V25:=NULL]
sensor_vals <- unique(sensor_vals)

samples <- sample_n(sensor_vals, 50)

write.table(samples, file = "sample_sensors-300.csv", col.names=FALSE, row.names=FALSE)








