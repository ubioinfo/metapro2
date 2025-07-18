n = 1000; 
d = 40; w0_40 = simul_weight(N = n, d = d, r = 0, e = 0)
d = 100; w0_100 = simul_weight(N = n, d = d, r = 0, e = 0)
d = 500; w0_500 = simul_weight(N = n, d = d, r = 0, e = 0)


n =200; d = 40;
w1_40 = simul_weight(N = n, d = d, r = 2, e = 0.40)
w2_40 = simul_weight(N = n, d = d, r = 6, e = 0.27)
w3_40 = simul_weight(N = n, d = d, r = 24, e = 0.14)

n = 200; d = 100;
w1_100 = simul_weight(N = n, d = d, r = 5, e = 0.32)
w2_100 = simul_weight(N = n, d = d, r = 15, e = 0.20)
w3_100 = simul_weight(N = n, d = d, r = 60, e = 0.10)

n = 200; d = 500;
w1_500 = simul_weight(N = n, d = d, r = 25, e = 0.20)
w2_500 = simul_weight(N = n, d = d, r = 75, e = 0.12)
w3_500 = simul_weight(N = n, d = d, r = 300, e = 0.06)


# Plot results
library(ggplot2)
library(reshape2)
library(gridExtra)
library(ggpubr)

plot_combined_results <- function(a1, a2, a3, compare_set = c("full", "reduced")) {
  compare_set <- match.arg(compare_set)
  
  if (compare_set == "reduced") {
    methods_to_plot <- c("wARTP", "GBJ", "GHC")
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
      labs(title = unique(df_sub$setting), y = "-log(p-value)", x = "") + coord_cartesian(ylim = c(0, 11)) +
      theme_minimal() + theme(legend.position = "none")
  })
  
  grid.arrange(grobs = plots, nrow = 1)  # 3 x 1
}


plot_combined_results(w1_40, w2_40, w3_40, compare_set = "full")
plot_combined_results(w1_40, w2_40, w3_40, compare_set = "reduced")

plot_combined_results(w1_100, w2_100, w3_100, compare_set = "full")
plot_combined_results(w1_100, w2_100, w3_100, compare_set = "reduced")

plot_combined_results(w1_500, w2_500, w3_500, compare_set = "full")
plot_combined_results(w1_500, w2_500, w3_500, compare_set = "reduced")