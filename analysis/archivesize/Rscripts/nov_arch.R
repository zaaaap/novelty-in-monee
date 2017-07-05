#!/usr/bin/env Rscript
library(readr)
library(data.table)
library(dplyr)
require(gridExtra)
require(reshape2)
require(ggplot2)
library(ggthemes)
library(scales)

sem = function(x) {
  return (sd(x)/sqrt(length(x)))
}

calc_mean_diff = function(data2) {
  diff2 = mean(data2[[4]]) - mean(data2[[3]])
  
  show(diff2)
  return (diff2)
}

calc_mean = function(data) {
  iterations = seq(1000, 1000000, 1000)
  
  data_mean = data[[1]]
  data_means_blocked <- matrix(data_mean, nrow = 1000, ncol = 10)
  avg = rowMeans(data_means_blocked, na.rm=TRUE)
  l = loess(avg~iterations, span = 0.25)
  smoothed = predict(l)
  mins = apply(data_means_blocked, 1, min)
  maxs = apply(data_means_blocked, 1, max)
  return (list(avg, smoothed, mins, maxs))
}

get_stats = function(data, name) {
  print(paste("--------- summary for", name, sep=" "))
  show(summary(data))
  show(shapiro.test(data)) #h0: samples from normal distribution
  # normality
}

multi.ttest <- function(data0, data1, data2, data3) {
  mat = matrix(, nrow=1000, ncol=0)
  mat <- cbind(mat, data0)
  mat <- cbind(mat, data1)
  mat <- cbind(mat, data2)
  mat <- cbind(mat, data3)
  mat <- as.matrix(mat)
  
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      test <- t.test(mat[, i], mat[, j])
      show(test)
      p.mat[i, j] <- p.mat[j, i] <- test$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  
  upper <- signif(p.mat,2)
  upper[upper.tri(upper)]<-""
  upper<-as.data.frame(upper)
  show(upper)
}

plotter = function(data0, data1, data2, data3, plotname, yname, varname) {
  iterations = seq(1000, 1000000, 1000)
  
  iterations = seq(1000, 1000000, 1000)
  
  # df <- data.table(it=iterations, d0=data0[[1]], d1=data1[[1]], d2=data2[[1]], d3=data3[[1]])
  df <- data.table(it=iterations, 
                   d0=data0[[1]], d0up = data0[[3]], d0low = data0[[4]], 
                   d1=data1[[1]], d1up = data1[[3]], d1low = data1[[4]], 
                   d2=data2[[1]], d2up = data2[[3]], d2low = data2[[4]], 
                   d3=data3[[1]], d3up = data3[[3]], d3low = data3[[4]]
  )
  
  myPalette <- c("#000000", "#FF0000", "#0000FF", "#32CD32")
  
  mean0 = mean(df$d0)
  mean1 = mean(df$d1)
  mean2 = mean(df$d2)
  mean3 = mean(df$d3)
  
  gtest <- ggplot(df, aes(x = it)) +
    # geom_hline(yintercept = mean1, linetype="dashed", color = 2)+
    geom_line(aes(y=d0, colour = "1"),  size = 0.3, alpha=0.5) +
    geom_line(aes(y=d1, colour = "20"),  size = 0.3, alpha=0.5) +
    # geom_hline(yintercept = mean2, linetype="dashed", color = 3)+
    geom_line(aes(y=d2, colour = "50"), size = 0.3, alpha=0.5) +
    # geom_hline(yintercept = mean3, linetype="dashed", color = 4)+
    geom_line(aes(y=d3, colour = "100"), size = 0.3, alpha=0.5) +
    geom_smooth(aes(y=d0, color = "1"), method=lm, formula = y ~ splines::bs(x), se = FALSE) +
    geom_smooth(aes(y=d0up, color = "1"), method=lm, formula = y ~ splines::bs(x), se = FALSE, linetype=2, size = 0.5) +
    geom_smooth(aes(y=d0low, color = "1"), method=lm, formula = y ~ splines::bs(x), se = FALSE, linetype=2, size = 0.5) +
    geom_smooth(aes(y=d1, color = "20"), method=lm, formula = y ~ splines::bs(x), se = FALSE) +
    geom_smooth(aes(y=d1up, color = "20"), method=lm, formula = y ~ splines::bs(x), se = FALSE, linetype=2, size = 0.5) +
    geom_smooth(aes(y=d1low, color = "20"), method=lm, formula = y ~ splines::bs(x), se = FALSE, linetype=2, size = 0.5) +
    geom_smooth(aes(y=d2, color = "50"), method=lm, formula = y ~ splines::bs(x), se = FALSE) +
    geom_smooth(aes(y=d2up, color = "50"), method=lm, formula = y ~ splines::bs(x), se = FALSE, linetype=2, size = 0.5) +
    geom_smooth(aes(y=d2low, color = "50"), method=lm, formula = y ~ splines::bs(x), se = FALSE, linetype=2, size = 0.5) +
    geom_smooth(aes(y=d3, color = "100"), method=lm, formula = y ~ splines::bs(x), se = FALSE)+
    geom_smooth(aes(y=d3up, color = "100"), method=lm, formula = y ~ splines::bs(x), se = FALSE, linetype=2, size = 0.5) +
    geom_smooth(aes(y=d3low, color = "100"), method=lm, formula = y ~ splines::bs(x), se = FALSE, linetype=2, size = 0.5) +
    scale_x_continuous(limits = c(0, 1000000), labels = comma) + 
    scale_y_continuous(limits = c(NA, NA)) +
    labs(x = "timesteps", y = yname, colour = varname) + 
    theme_linedraw() + 
    scale_colour_manual(values=myPalette) +
    theme(legend.justification = c(1.1, 1.1), legend.position = c(1, 1), legend.direction = "horizontal",
          plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))

  
  # Save plots ---
  ggsave(plotname, plot = gtest, device = pdf, path = ".", width = 120, height = 120, units = "mm")
  
}

setwd("~/git/monee-novelty/analysis/archivesize")

plotname = "arch_nov.pdf"
yname = "avg. novelty"
varname = "archive size"

d0 <- fread("div_729_1.csv")
d1 <- fread("div_729_20.csv")
d2 <- fread("div_729_50.csv")
d3 <- fread("div_729_100.csv")

# Read in data and transform -------
d0meaned = calc_mean(d0)
d1meaned = calc_mean(d1)
d2meaned = calc_mean(d2)
d3meaned = calc_mean(d3)

# Statistics ----
get_stats(d0meaned[[1]], "d1")
get_stats(d1meaned[[1]], "d20")
get_stats(d2meaned[[1]], "d50")
get_stats(d3meaned[[1]], "d100")

calc_mean_diff(d0meaned)
calc_mean_diff(d1meaned)
calc_mean_diff(d2meaned)
calc_mean_diff(d3meaned)



result = multi.ttest(d0meaned[[1]], d1meaned[[1]], d2meaned[[1]], d3meaned[[1]])

# Plots ----
plotter(d0meaned, d1meaned, d2meaned, d3meaned, plotname, yname, varname)



