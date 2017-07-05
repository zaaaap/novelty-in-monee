plotter = function(data0, data1, data2, data3, plotname, yname, varname) {
  iterations = seq(1000, 1000000, 1000)
  
  # df <- data.table(it=iterations, d0=data0[[1]], d1=data1[[1]], d2=data2[[1]], d3=data3[[1]])
  df <- data.table(it=iterations, 
                   d0=data0[[1]], d0up = data0[[3]], d0low = data0[[4]], d0sm = data0[[2]], 
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
    geom_line(aes(y = d0, colour = "1"),  size = 0.3, alpha=0.5) +
    # geom_hline(yintercept = mean1, linetype="dashed", color = 2)+
    geom_line(aes(y=d1, colour = "20"),  size = 0.3, alpha=0.5) +
    # geom_hline(yintercept = mean2, linetype="dashed", color = 3)+
    geom_line(aes(y=d2, colour = "50"), size = 0.3, alpha=0.5) +
    # geom_hline(yintercept = mean3, linetype="dashed", color = 4)+
    geom_line(aes(y=d3, colour = "100"), size = 0.3, alpha=0.5) +
    geom_smooth(aes(y=d0, color = "1"), method=lm, formula = y ~ splines::bs(x), se = FALSE) +
    geom_smooth(aes(y=d0up, color = "1"), method=lm, formula = y ~ splines::bs(x), se = FALSE, linetype="dotted", size = 0.5) +
    geom_smooth(aes(y=d0low, color = "1"), method=lm, formula = y ~ splines::bs(x), se = FALSE, linetype="dotted", size = 0.5) +
    geom_smooth(aes(y=d1, color = "20"), method=lm, formula = y ~ splines::bs(x), se = FALSE) +
    geom_smooth(aes(y=d1up, color = "20"), method=lm, formula = y ~ splines::bs(x), se = FALSE, linetype="dotted", size = 0.5) +
    geom_smooth(aes(y=d1low, color = "20"), method=lm, formula = y ~ splines::bs(x), se = FALSE, linetype="dotted", size = 0.5) +
    geom_smooth(aes(y=d2, color = "50"), method=lm, formula = y ~ splines::bs(x), se = FALSE) +
    geom_smooth(aes(y=d2up, color = "50"), method=lm, formula = y ~ splines::bs(x), se = FALSE, linetype="dotted", size = 0.5) +
    geom_smooth(aes(y=d2low, color = "50"), method=lm, formula = y ~ splines::bs(x), se = FALSE, linetype="dotted", size = 0.5) +
    geom_smooth(aes(y=d3, color = "100"), method=lm, formula = y ~ splines::bs(x), se = FALSE)+
    geom_smooth(aes(y=d3up, color = "100"), method=lm, formula = y ~ splines::bs(x), se = FALSE, linetype="dotted", size = 0.5) +
    geom_smooth(aes(y=d3low, color = "100"), method=lm, formula = y ~ splines::bs(x), se = FALSE, linetype="dotted", size = 0.5) +
    scale_x_continuous(limits = c(0, 1000000)) + 
    scale_y_continuous(limits = c(NA, NA)) +
    labs(x = "timesteps", y = yname, colour = varname) + 
    theme_linedraw() + 
    scale_colour_manual(values=myPalette) +
    theme(legend.justification = c(1.1, 1.1), legend.position = c(1, 1), legend.direction = "horizontal")
  
  gtest2 = ggplot(df) +
    geom_density(aes(d0, fill = "1"), alpha = 0.5, size = 0.2, trim=F)+
    geom_density(aes(d1, fill = "20"), alpha = 0.5, size = 0.2, trim=F)+
    geom_density(aes(d2, fill = "50"), alpha = 0.5, size = 0.2, trim=F)+
    geom_density(aes(d3, fill = "100"), alpha = 0.5, size = 0.2, trim=F)+
    coord_flip() + xlim(0,NA) +
    scale_fill_manual(values=myPalette) + 
    scale_colour_manual(values=myPalette) +
    theme_void() + theme(legend.position = "none", plot.margin = unit(c(8.3,0,3,0), "lines"))
  
  arr = grid.arrange(gtest, gtest2, ncol=2, nrow=1, widths=c(4, 1.5), heights=c(4))
  
  # Save plots ---
  ggsave(plotname, plot = arr, device = pdf, path = ".", width = 160, height = 100, units = "mm")
  
}