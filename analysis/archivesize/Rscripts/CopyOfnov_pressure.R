#!/usr/bin/env Rscript
library(readr)
library(data.table)
library(dplyr)
library(ggplot2)
require(gridExtra)
require(reshape2)
library(ggthemes)
library(car)
library(gridExtra)

# define functions ----

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
  
  upper <- p.mat
  upper[upper.tri(upper)]<-""
  upper<-as.data.frame(upper)
  show(upper)
}

plotter = function(data0, data1, data2, data3, plotname, yname, varname) {
  iterations = seq(1000, 1000000, 1000)
  
  df <- data.table(it=iterations, d0=data0, d1=data1, d2=data2, d3=data3)
  
  myPalette <- c("#000000", "#FF0000", "#0000FF", "#32CD32")
  mean0 = mean(df$d0)
  mean1 = mean(df$d1)
  mean2 = mean(df$d2)
  mean3 = mean(df$d3)
  
  gtest <- ggplot(df, aes(x = it)) +
    geom_line(aes(y = d0, colour = "1"),  size = 0.3, alpha=0.5) +
    # geom_hline(yintercept = mean1, linetype="dashed", color = 2)+
    geom_line(aes(y=d1, colour = "20"),  size = 0.3, alpha=0.5) +
    # geom_hline(yintercept = mean2, linetype="dashed", color = 3)+
    geom_line(aes(y=d2, colour = "50"), size = 0.3, alpha=0.5) +
    # geom_hline(yintercept = mean3, linetype="dashed", color = 4)+
    geom_line(aes(y=d3, colour = "100"), size = 0.3, alpha=0.5) +
    geom_smooth(aes(y=d0, color = "1"), method=lm, formula = y ~ splines::bs(x), se = FALSE) +
    geom_smooth(aes(y=d1, color = "20"), method=lm, formula = y ~ splines::bs(x), se = FALSE) +
    geom_smooth(aes(y=d2, color = "50"), method=lm, formula = y ~ splines::bs(x), se = FALSE) +
    geom_smooth(aes(y=d3, color = "100"), method=lm, formula = y ~ splines::bs(x), se = FALSE)+
    scale_x_continuous(limits = c(0, 1000000)) + 
    scale_y_continuous(limits = c(NA, NA)) +
    labs(x = "timesteps", y = yname, colour = varname) + 
    theme_linedraw() + 
    # scale_colour_brewer(palette = "Set1") +
    scale_colour_manual(values=myPalette) +
    theme(legend.justification = c(1.1, 1.1), legend.position = c(1, 1), legend.direction = "horizontal")
  
  gtest2 = ggplot(df) +
    geom_vline(xintercept = mean(d0), linetype="dashed", colour = "black")+
    geom_vline(xintercept = mean(d1), linetype="dashed", colour = "#0000FF")+
    geom_vline(xintercept = mean(d2), linetype="dashed", colour = "#32CD32")+
    geom_vline(xintercept = mean(d3), linetype="dashed", colour = "#FF0000")+
    geom_density(aes(d0, fill = "1"), alpha = 0.5, size = 0.3, trim=F)+
    geom_density(aes(d1, fill = "20"), alpha = 0.5, size = 0.3, trim=F)+
    geom_density(aes(d2, fill = "50"), alpha = 0.5, size = 0.3, trim=F)+
    geom_density(aes(d3, fill = "100"), alpha = 0.5, size = 0.3, trim=F)+
    coord_flip() +
    scale_fill_manual(values=myPalette) + 
    scale_colour_manual(values=myPalette) +
    theme_void() + theme(legend.position = "none", plot.margin = unit(c(0,0,2,0), "lines"))
  
  arr = grid.arrange(gtest, gtest2, ncol=2, nrow=1, widths=c(4, 1.5), heights=c(4))
  show(gtest)
  
  # Save plots ---
  ggsave(plotname, plot = arr, device = pdf, path = ".", width = 160, height = 100, units = "mm")
}

calc_pressure <- function(windata, inddata){
  ps <- rep(0,10000)
  # 1000 times 100 ind stats for 10 runs = 1mio lines
  # -> still have to run that 10 times 
  # add offspring coutn to individual stats
  for(j in 0:9) {
    print(j)
    
    for (i in 0:999) {
      a <- j*100000 + i*100 + 1
      z <- j*100000 + i*100 + 100
      
      currentdd <- inddata[a:z,]
      current_count <- as.data.table(table(as.numeric(windata[i+1,]), useNA = "no"))
      current_count$V1 <- as.numeric(current_count$V1)
      currentdd <- merge(currentdd, current_count, by = "V1", all.x = TRUE)
      currentdd[is.na(currentdd)] <- 0
      # choose V:    1 = id, 2 = dist, 3 = pucks, 4 = nov, 5 = offspring
      mean_pucks <- mean(currentdd$V4)
      mean_offspring <- mean(currentdd$N)
      #create contingency table and calc pressure
      cont <- matrix(c(nrow(currentdd[which(currentdd$N <= mean_offspring & currentdd$V4 <= mean_pucks)]),
                       nrow(currentdd[which(currentdd$N > mean_offspring & currentdd$V4 <= mean_pucks)]),
                       nrow(currentdd[which(currentdd$N <= mean_offspring & currentdd$V4 > mean_pucks)]),
                       nrow(currentdd[which(currentdd$N > mean_offspring & currentdd$V4 > mean_pucks)])),
                     nrow = 2, ncol = 2)
      
      if(cont[2,1] > 0) {
        cont[2,1] = cont[2,1] - 1
        cont[1,1] = cont[1,1] + 1
      } 
      if(cont[1,2] > 0) {
        cont[1,2] = cont[1,2] - 1
        cont[2,2] = cont[2,2] + 1
      }
      
      ps[j*1000 + i + 1] <- fisher.test(cont)$p.value
    }  
  }
  ps2 <- -1 * log(ps*2)
  means_blocked <- matrix(ps2, nrow = 1000, ncol = 10)
  avg_pressure <- rowMeans(means_blocked, na.rm=TRUE)
  return(avg_pressure)
}

calc_pressure2 <- function(pooldata){
  time = 1000
  
  for(i in 1:nrow(pooldata)) {
    temp_m = matrix(, nrow=0, ncol=6)
    temp = pooldata[i,]
    for(j in 1:length(temp)) {
      if(!is.na(temp[j])) {
        ind_data = read_csv(paste(temp[j], "\n"), col_names = )
        temp_m = rbind(temp_m, ind_data)
      }
    }
    temp_m = as.data.frame(temp_m)   
    mean_nov <- mean(temp_m$X5)
    mean_offspring <- mean(temp_m$X6)
    show(mean_nov)
    # cont <- matrix(c(nrow(temp_m[which(temp_m$N <= mean_offspring & currentdd$V4 <= mean_pucks)]),
    #                  nrow(currentdd[which(currentdd$N > mean_offspring & currentdd$V4 <= mean_pucks)]),
    #                  nrow(currentdd[which(currentdd$N <= mean_offspring & currentdd$V4 > mean_pucks)]),
    #                  nrow(currentdd[which(currentdd$N > mean_offspring & currentdd$V4 > mean_pucks)])),
    #                nrow = 2, ncol = 2)
    # 
    # if(cont[2,1] > 0) {
    #   cont[2,1] = cont[2,1] - 1
    #   cont[1,1] = cont[1,1] + 1
    # }
    # if(cont[1,2] > 0) {
    #   cont[1,2] = cont[1,2] - 1
    #   cont[2,2] = cont[2,2] + 1

  }
}


# Read in data -----
setwd("~/git/monee-novelty/analysis/archivesize")

wd0 <- read_delim("poolstats_729_20.csv", delim = ";", col_names = F)
wd0 = as.matrix(wd0)


calc_pressure2(wd0[1:100,])



