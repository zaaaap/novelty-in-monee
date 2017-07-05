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

calc_mean_dist = function(data) {
  iterations = seq(1000, 1000000, 1000)
  
  # 1 = distance, 2 life, 3 puck1, 4 puck2, 5 pool, 6 archiveadd
  data_mean = data[[1]] #define distance
  data_means_blocked <- matrix(data_mean, nrow = 1000, ncol = 10)
  avg = rowMeans(data_means_blocked, na.rm=TRUE)
  l = loess(avg~iterations, span = 0.25)
  smoothed = predict(l)
  mins = apply(data_means_blocked, 1, min)
  maxs = apply(data_means_blocked, 1, max)
  return (list(avg, smoothed, mins, maxs))
}

calc_mean_nov = function(data) {
  iterations = seq(1000, 1000000, 1000)
  
  data_mean = rowMeans(data, na.rm=TRUE)
  data_means_blocked <- matrix(data_mean, nrow = 1000, ncol = 10)
  avg = rowMeans(data_means_blocked, na.rm=TRUE)
  l = loess(avg~iterations, span = 0.25)
  smoothed = predict(l)
  mins = apply(data_means_blocked, 1, min)
  maxs = apply(data_means_blocked, 1, max)
  return (list(avg, smoothed, mins, maxs))
}

calc_mean_pucks = function(data) {
  iterations = seq(1000, 1000000, 1000)
  
  # 1 = distance, 2 life, 3 puck1, 4 puck2, 5 pool, 6 archiveadd
  data_mean = data[[3]] + data[[4]] #define distance
  data_means_blocked <- matrix(data_mean, nrow = 1000, ncol = 10)
  avg = rowMeans(data_means_blocked, na.rm=TRUE)
  l = loess(avg~iterations, span = 0.25)
  smoothed = predict(l)
  mins = apply(data_means_blocked, 1, min)
  maxs = apply(data_means_blocked, 1, max)
  return (list(avg, smoothed, mins, maxs))
}


testplot_cor = function(nov, distORpucks, yname, plotname, archsize, xlimit, ylimit) {
  # in steps of 100k
  iterations = sort(rep(seq(100000,1000000,100000), 10))
  dcor = data.table(it=iterations, nov=nov[[1]], var2=distORpucks[[1]])
  
  result = cor.test(dcor$nov, dcor$var2, method = "spearman", conf.level = 0.95)
  show(result)
  
  rlabel = paste("R =", round(result$estimate, 2))
  archlabel = paste("archive size = ", archsize)
  
  gcor = ggplot(dcor, aes(var2, nov)) + 
    geom_point(aes(colour=it), size=0.5) + 
    geom_smooth(method="lm", se=F, color = "#0000FF") +
    scale_colour_gradient(low = "lightgrey", high = "black") +
    theme_linedraw() +
    scale_x_continuous(limits = c(0, xlimit)) + 
    scale_y_continuous(limits = c(0, ylimit)) +
    labs(y = "avg. novelty", x = yname, colour = "time-steps") +
    theme(aspect.ratio=1) +
    geom_text(x =250, y =4, aes(label = rlabel), hjust = 0)
    
  
  # show(gcor)
  # ggsave(plotname, plot = gcor, device = pdf, path = ".", width = 100, height = 100, units = "mm")
  return(gcor)
}

setwd("~/git/monee-novelty/analysis/archivesize")

plotname = "cor_nd0.pdf"
yname = "avg. distance travelled"

d0 <- fread("genstats_729_1.csv")
d1 <- fread("genstats_729_20.csv")
d2 <- fread("genstats_729_50.csv")
d3 <- fread("genstats_729_100.csv")

n0 <- fread("div_729_1.csv")
n1 <- fread("div_729_20.csv")
n2 <- fread("div_729_50.csv")
n3 <- fread("div_729_100.csv")

# Read in data and transform -------
p0meaned = calc_mean_pucks(d0)
p1meaned = calc_mean_pucks(d1)
p2meaned = calc_mean_pucks(d2)
p3meaned = calc_mean_pucks(d3)

n0meaned = calc_mean_nov(n0)
n1meaned = calc_mean_nov(n1)
n2meaned = calc_mean_nov(n2)
n3meaned = calc_mean_nov(n3)

d0meaned = calc_mean_dist(d0)
d1meaned = calc_mean_dist(d1)
d2meaned = calc_mean_dist(d2)
d3meaned = calc_mean_dist(d3)

# DO the correlation thingy -----
# is there a correlation with novelty ?

nd0plot = testplot_cor(n0meaned, d0meaned, "avg. distance travelled", "cor_nd0.pdf", "1", 3500, 25)
nd1plot = testplot_cor(n1meaned, d1meaned, "avg. distance travelled", "cor_nd1.pdf", "20", 3500, 25)
nd2plot = testplot_cor(n2meaned, d2meaned, "avg. distance travelled", "cor_nd2.pdf", "50", 3500, 25)
nd3plot = testplot_cor(n3meaned, d3meaned, "avg. distance travelled", "cor_nd3.pdf", "100", 3500, 25)

arr_dist = grid.arrange(nd0plot, nd1plot, nd2plot,nd3plot, ncol=2, nrow=2, widths=c(4, 4), heights=c(4,4))
show(arr_dist)

np0plot = testplot_cor(n0meaned, p0meaned, "no. pucks collected", "cor_np0.pdf", "1", 1000, 25)
np1plot = testplot_cor(n1meaned, p1meaned, "no. pucks collected", "cor_np1.pdf", "20", 1000, 25)
np2plot = testplot_cor(n2meaned, p2meaned, "no. pucks collected", "cor_np2.pdf", "50", 1000, 25)
np3plot = testplot_cor(n3meaned, p3meaned, "no. pucks collected", "cor_np3.pdf", "100", 1000, 25)


arr_both = grid.arrange(np0plot, nd0plot, ncol=2, nrow=1, widths=c(4, 4))
arr_pucks = grid.arrange(np0plot, np1plot, np2plot,np3plot, ncol=2, nrow=2, widths=c(4, 4), heights=c(4, 4))
# ggsave("arch_both_cor.pdf", plot = arr_both, device = pdf, path = ".", width = 200, height = 100, units = "mm")
ggsave("arch_pucks_cor.pdf", plot = arr_pucks, device = pdf, path = ".", width = 200, height = 100, units = "mm")
ggsave("arch_dist_cor.pdf", plot = arr_dist, device = pdf, path = ".", width = 200, height = 100, units = "mm")
show(arr_pucks)

