#!/usr/bin/env Rscript
library(readr)
library(data.table)
library(dplyr)
require(gridExtra)
require(reshape2)
require(ggplot2)
library(ggthemes)

sem = function(x) {
  return (sd(x)/sqrt(length(x)))
}

calc_mean = function(data) {
  iterations = seq(1000, 1000000, 1000)
  
  # 1 = distance, 2 life, 3 puck1, 4 puck2, 5 pool, 6 archiveadd
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

multi.ttest <- function(data1, data2, data3) {
  mat = matrix(, nrow=1000, ncol=0)
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

plot3 = function(data1, data2, data3, plotname, yname, varname) {
  iterations = seq(1000, 1000000, 1000)
  
  df <- data.table(it=iterations, sm1=data1[[2]], d1=data1[[1]], sm2=data2[[2]], d2=data2[[1]], sm3=data3[[2]], d3=data3[[1]])
  
  mean1 = mean(df$d1)
  mean2 = mean(df$d2)
  mean3 = mean(df$d3)
  
  gtest <- ggplot(df, aes(x = it)) +
    # geom_hline(yintercept = mean(mean1), linetype="dashed", color = 2)+
    geom_line(aes(y = d1, colour = "20"),  size = 0.5, alpha=0.2) +
    geom_smooth(aes(y=d1, color = "20"), method=lm, formula = y ~ splines::bs(x), se = FALSE) +
    # geom_hline(yintercept = mean(mean2), linetype="dashed", color = 3)+
    geom_line(aes(y = d2, colour = "50"), alpha = 0.5, alpha=0.2) +
    geom_smooth(aes(y=d2, color = "50"), method=lm, formula = y ~ splines::bs(x), se = FALSE) +
    # geom_hline(yintercept = mean(mean3), linetype="dashed", color = 4)+
    geom_line(aes(y = d3, colour = "100"), alpha = 0.5, alpha=0.2) +
    geom_smooth(aes(y=d3, color = "100"), method=lm, formula = y ~ splines::bs(x), se = FALSE)+
    scale_x_continuous(limits = c(0, 1000000)) + 
    scale_y_continuous(limits = c(0, NA)) +
    labs(x = "timesteps", y = yname, colour = varname) + 
    theme_linedraw() +
    theme(legend.justification = c(1.1, 1.1), legend.position = c(1, 1), legend.direction = "horizontal")
  
  gtest2 = ggplot(df) +
    geom_vline(xintercept = mean(mean1), linetype="dashed", color = 3)+
    geom_vline(xintercept = mean(mean2), linetype="dashed", color = 4)+
    geom_vline(xintercept = mean(mean3), linetype="dashed", color = 2)+
    geom_density(aes(d1, fill = "20"), alpha = 0.5, size = 0.2, trim=F)+
    geom_density(aes(d2, fill = "50"), alpha = 0.5, size = 0.2, trim=F)+
    geom_density(aes(d3, fill = "100"), alpha = 0.5, size = 0.2, trim=F)+
    coord_flip() + xlim(0,NA) +
    theme_void() + theme(legend.position = "none", plot.margin = unit(c(0,0,2,0), "lines"))
  
  arr = grid.arrange(gtest, gtest2, ncol=2, nrow=1, widths=c(4, 1.5), heights=c(4))
  show(arr)
  
  # Save plots ---
  ggsave(plotname, plot = arr, device = pdf, path = ".", width = 160, height = 100, units = "mm")
  
}

setwd("~/git/monee-novelty/analysis/samplesize")

plotname = "test.pdf"
yname = "avg test"
varname = "testvar"

d1 <- fread("div_729_20.csv")
d2 <- fread("div_729_50.csv")
d3 <- fread("div_729_100.csv")

# Read in data and transform -------
d1meaned = calc_mean(d1)
d2meaned = calc_mean(d2)
d3meaned = calc_mean(d3)

# Statistics ----
get_stats(d1meaned[[1]], "d20")
get_stats(d2meaned[[1]], "d50")
get_stats(d3meaned[[1]], "d100")

result = multi.ttest(d1meaned[[1]], d2meaned[[1]], d3meaned[[1]])

# Plots ----
plot3(d1meaned, d2meaned, d3meaned, plotname, yname, varname)

