#!/usr/bin/env Rscript
library(readr)
library(data.table)
library(dplyr)
require(gridExtra)
require(reshape2)
require(ggplot2)
library(ggthemes)

calc_mean = function(data) {
  # 1 = distance, 2 life, 3 puck1, 4 puck2, 5 pool, 6 archiveadd
  data_mean = data[[6]] #define distance
  data_means_blocked <- matrix(data_mean, nrow = 1000, ncol = 10)
  avg = rowMeans(data_means_blocked, na.rm=TRUE) / 100
  avg_overall = mean(avg)

  show(avg_overall)
  return (avg_overall)
}

setwd("~/git/monee-novelty/analysis/archivesize")


d0 <- fread("genstats_729_1.csv")
d1 <- fread("genstats_729_20.csv")
d2 <- fread("genstats_729_50.csv")
d3 <- fread("genstats_729_100.csv")

# Read in data and transform -------
d0meaned = calc_mean(d0)
d1meaned = calc_mean(d1)
d2meaned = calc_mean(d2)
d3meaned = calc_mean(d3)

