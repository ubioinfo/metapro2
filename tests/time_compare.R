library(GBJ)
library(ggplot2)
library(dplyr)

time_compare_mixed <- function(dims = c(40, 100, 500)) {
  times <- data.frame()
  
  for (d in dims) {
    for (i in 1:10) {
      cat("wARTP repeat:", i, "dim:", d, "\n")
      z <- rnorm(d)
      p <- 2 * (1 - pnorm(abs(z))) / 5
      weight <- c(rep(1, d - 1), 1.1)
      
      # Unweighted B = 50K
      t1 <- system.time({
        res <- tryCatch(wARTP(p, cor_mat = diag(d), B2 = 50000)$p_value, error = function(e) NA)
      })[3]
      if (!is.na(res)) {
        times <- rbind(times, data.frame(method = "ARTP_50K", dimension = d, time = t1))
      }
      
      # Unweighted B = 500K
      t2 <- system.time({
        res <- tryCatch(wARTP(p, cor_mat = diag(d), B2 = 500000)$p_value, error = function(e) NA)
      })[3]
      if (!is.na(res)) {
        times <- rbind(times, data.frame(method = "ARTP_500K", dimension = d, time = t2))
      }
      
      # Weighted B = 50K
      t3 <- system.time({
        res <- tryCatch(wARTP(p, cor_mat = diag(d), weight = weight, B2 = 50000)$p_value, error = function(e) NA)
      })[3]
      if (!is.na(res)) {
        times <- rbind(times, data.frame(method = "wARTP_50K", dimension = d, time = t3))
      }
      
      # Weighted B = 500K
      t4 <- system.time({
        res <- tryCatch(wARTP(p, cor_mat = diag(d), weight = weight, B2 = 500000)$p_value, error = function(e) NA)
      })[3]
      if (!is.na(res)) {
        times <- rbind(times, data.frame(method = "wARTP_500K", dimension = d, time = t4))
      }
    }
    
    for (i in 1:10) {
      cat("Other methods repeat:", i, "dim:", d, "\n")
      z <- rnorm(d)
      p <- 2 * (1 - pnorm(abs(z))) / 5
      
      corr <- 0.1
      cor_mat <- matrix(corr, nrow = d, ncol = d)
      diag(cor_mat) <- 1
      
      t5 <- system.time({
        res <- tryCatch(GBJ(z, cor_mat = cor_mat)$GBJ_pvalue, error = function(e) NA)
      })[3]
      times <- rbind(times, data.frame(method = "GBJ", dimension = d, time = t5))
      
      t6 <- system.time({
        res <- tryCatch(wFisher(p), error = function(e) NA)
      })[3]
      times <- rbind(times, data.frame(method = "wFisher", dimension = d, time = t6))
      
      t7 <- system.time({
        res <- tryCatch(ordmeta(p, B = 50000)$p.value, error = function(e) NA)
      })[3]
      times <- rbind(times, data.frame(method = "ordmeta", dimension = d, time = t7))
    }
  }
  
  return(times)
}

# Run
set.seed(1)
result <- time_compare_mixed(dims = c(40, 100, 500))

# Summarize & Plot
summary_df <- result %>%
  mutate(log_time = log10(1 + time)) %>%
  group_by(method, dimension) %>%
  summarise(
    mean_time = mean(log_time),
    sd_time = sd(log_time),
    .groups = "drop"
  )

# Plotting setup
summary_df$method <- factor(summary_df$method,
                            levels = c("ARTP_50K", "ARTP_500K",
                                       "wARTP_50K", "wARTP_500K",
                                       "GBJ", "ordmeta", "wFisher"))

method_colors <- c(
  "ARTP_50K" = "orange",
  "ARTP_500K" = "orange",
  "wARTP_50K" = "darkred",
  "wARTP_500K" = "darkred",
  "GBJ" = "blue",
  "ordmeta" = "purple",
  "wFisher" = "forestgreen"
)

method_linetypes <- c(
  "ARTP_50K" = "dashed",
  "ARTP_500K" = "solid",
  "wARTP_50K" = "dashed",
  "wARTP_500K" = "solid",
  "GBJ" = "solid",
  "ordmeta" = "solid",
  "wFisher" = "solid"
)

method_shapes <- c(
  "ARTP_50K" = 16,
  "ARTP_500K" = 17,
  "wARTP_50K" = 15,
  "wARTP_500K" = 18,
  "GBJ" = 16,
  "ordmeta" = 16,
  "wFisher" = 16
)

library(scales)
ggplot(summary_df, aes(x = dimension, y = mean_time,
                       color = method, linetype = method, shape = method, group = method)) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(aes(ymin = mean_time - sd_time,
                    ymax = mean_time + sd_time),
                width = 5, alpha = 0.6) +
  scale_color_manual(values = method_colors) +
  scale_linetype_manual(values = method_linetypes) +
  scale_shape_manual(values = method_shapes) +
  scale_x_continuous(breaks = c(40, 100, 500)) +
  scale_y_continuous(
    name = "Computation Time (sec)",
    breaks = log10(1 + c(0.1, 1, 10, 100, 300)),
    labels = c(0.1, 1, 10, 100, 300)
  ) +
  labs(x = "Number of p-values combined",
       color = "Method", linetype = "Method", shape = "Method") +
  theme_minimal(base_size = 18) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 15),
    legend.key.width = unit(1.5, "cm")
  )