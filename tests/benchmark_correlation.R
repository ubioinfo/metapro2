n = 100; d = 40;

rho = 0.2; 
c10_40 = simul_corr(N = n, d = d, r = 0, e = 0, rho = rho, cortype = "const")
c11_40 = simul_corr(N = n, d = d, r = 3, e = 2.7, rho = rho, cortype = "const")
c12_40 = simul_corr(N = n, d = d, r = 15, e = 1.7, rho = rho, cortype = "const")
c13_40 = simul_corr(N = n, d = d, r = 24, e = 1.5, rho = rho, cortype = "const")

rho = 0.5; 
c20_40 = simul_corr(N = n, d = d, r = 0, e = 0, rho = rho, cortype = "const")
c21_40 = simul_corr(N = n, d = d, r = 3, e = 2.8, rho = rho, cortype = "const")
c22_40 = simul_corr(N = n, d = d, r = 15, e = 1.8, rho = rho, cortype = "const")
c23_40 = simul_corr(N = n, d = d, r = 24, e = 1.6, rho = rho, cortype = "const")

rho = 0.8; 
c30_40 = simul_corr(N = n, d = d, r = 0, e = 0, rho = rho, cortype = "const")
c31_40 = simul_corr(N = n, d = d, r = 3, e = 2.8, rho = rho, cortype = "const")
c32_40 = simul_corr(N = n, d = d, r = 15, e = 2, rho = rho, cortype = "const")
c33_40 = simul_corr(N = n, d = d, r = 24, e = 1.7, rho = rho, cortype = "const")


n = 100; d = 100;

rho = 0.2;
c10_100 = simul_corr(N = n, d = d, r = 0, e = 0, rho = rho, cortype = "const")
c11_100 = simul_corr(N = n, d = d, r = 5, e = 2.4, rho = rho, cortype = "const")
c12_100 = simul_corr(N = n, d = d, r = 15, e = 1.7, rho = rho, cortype = "const")
c13_100 = simul_corr(N = n, d = d, r = 60, e = 1.2, rho = rho, cortype = "const")

rho = 0.5; 
c20_100 = simul_corr(N = n, d = d, r = 0, e = 0, rho = rho, cortype = "const")
c21_100 = simul_corr(N = n, d = d, r = 5, e = 2.8, rho = rho, cortype = "const")
c22_100 = simul_corr(N = n, d = d, r = 15, e = 1.8, rho = rho, cortype = "const")
c23_100 = simul_corr(N = n, d = d, r = 60, e = 1.3, rho = rho, cortype = "const")

rho = 0.8; 
c30_100 = simul_corr(N = n, d = d, r = 0, e = 0, rho = rho, cortype = "const")
c31_100 = simul_corr(N = n, d = d, r = 5, e = 2.8, rho = rho, cortype = "const")
c32_100 = simul_corr(N = n, d = d, r = 15, e = 1.9, rho = rho, cortype = "const")
c33_100 = simul_corr(N = n, d = d, r = 60, e = 1.4, rho = rho, cortype = "const")


# Combine 4 simulation results into 2x2 boxplots
library(ggplot2)
library(reshape2)
library(gridExtra)  # for grid.arrange
library(ggplot2)
library(reshape2)
library(gridExtra)  # for grid.arrange

plot_combined_results <- function(a1, a2, a3, a4) {
  to_df_long <- function(pv_matrix, label) {
    df <- as.data.frame(-log(pv_matrix))
    df$sim <- 1:nrow(df)
    df_long <- melt(df, id.vars = "sim", variable.name = "method", value.name = "neglogp")
    df_long$setting <- label
    return(df_long)
  }
  
  df1 <- to_df_long(a1, "%True Effects = 0%")    # (1,1)
  df2 <- to_df_long(a2, "%True Effects = 5%")    # (1,2)
  df3 <- to_df_long(a3, "%True Effects = 15%")   # (2,1)
  df4 <- to_df_long(a4, "%True Effects = 60%")   # (2,2)
  
  plots <- list(df1, df2, df3, df4); #h = 14
  plots <- lapply(plots, function(df_sub) {
    ggplot(df_sub, aes(x = method, y = neglogp, fill = method)) +
      geom_boxplot(width = 0.4, outlier.shape = 8, outlier.size = 2) + 
      #geom_jitter(width = 0.2, size = 0.6, alpha = 0.4) +
      geom_hline(yintercept = -log(0.05), color = "red", linetype = "dashed") +
      labs(title = unique(df_sub$setting), y = "-log(p-value)") + #ylim(0, h) +
      theme_minimal() + theme(legend.position = "none")
  })
  grid.arrange(grobs = plots, ncol = 2)
}

## Three figures
library(ggplot2)
library(reshape2)
library(gridExtra)

plot_combined_results <- function(a1, a2, a3, compare_set = c("full", "reduced")) {
  compare_set <- match.arg(compare_set)
  
  if (compare_set == "reduced") {
    methods_to_plot <- c("wARTP", "wARTPgam", "wARTPlog", "ARTP")
  } else {
    methods_to_plot <- c("wARTP", "ARTP", "wFisher", "Fisher", "ordmeta", "GBJ", "GHC", "minP")
  }
  
  to_df_long <- function(pv_matrix, label) {
    df <- as.data.frame(-log(pv_matrix))
    df$sim <- 1:nrow(df)
    df_long <- melt(df, id.vars = "sim", variable.name = "method", value.name = "neglogp")
    df_long <- df_long[df_long$method %in% methods_to_plot, ]
    df_long$setting <- label
    return(df_long)
  }
  
  df1 <- to_df_long(a1, "%True Effects = 5%")
  df2 <- to_df_long(a2, "%True Effects = 15%")
  df3 <- to_df_long(a3, "%True Effects = 60%")
  
  plots <- list(df1, df2, df3)
  
  plots <- lapply(plots, function(df_sub) {
    ggplot(df_sub, aes(x = method, y = neglogp, fill = method)) +
      geom_boxplot(width = 0.4, outlier.shape = 8, outlier.size = 2) +
      geom_hline(yintercept = -log(0.05), color = "red", linetype = "dashed") +
      labs(title = unique(df_sub$setting), y = "-log(p-value)", x = "") +
      theme_minimal() + theme(legend.position = "none")
  })
  
  grid.arrange(grobs = plots, ncol = 1)  # 3 x 1
}


#plot_combined_results(c00_40, c10_40, c20_40, c30_40)
plot_combined_results(c01_40, c02_40, c03_40, c04_40)
plot_combined_results(c11_40, c12_40, c13_40, c14_40)
plot_combined_results(c21_40, c22_40, c23_40, c24_40)
plot_combined_results(c31_40, c32_40, c33_40, c34_40)

#plot_combined_results(c00_100, c10_100, c20_100, c30_100)
plot_combined_results(c01_100, c02_100, c03_100, c04_100)
plot_combined_results(c11_100, c12_100, c13_100, c14_100)
plot_combined_results(c21_100, c22_100, c23_100, c24_100)
plot_combined_results(c31_100, c32_100, c33_100, c34_100)

