simul_corr <- function(N = 15, d = 50, r = 10, e = 2, rho = 0.85, cortype = c("pd", "const")) {
  library(GBJ)
  library(mvtnorm)
  library(ggplot2)
  library(reshape2)
  library(parallel)
  
  cortype <- match.arg(cortype)
  
  pd_matrix <- function(n, rho = 0.5) {
    A <- matrix(rnorm(n * n), n, n)
    A <- crossprod(A) + rho * diag(n)
    cov2cor(A)
  }
  
  const_matrix <- function(n, rho = 0.5) {
    mat <- matrix(rho, n, n)
    diag(mat) <- 1
    return(mat)
  }
  
  base_cor_mat <- switch(cortype,
                         pd = NULL,
                         const = const_matrix(d, rho))
  
  sim_results <- mclapply(1:N, function(i) {
    tryCatch({
      true_effect <- sample(c(rep(e, r), rep(0, d - r)))
      
      cor_mat <- if (cortype == "pd") pd_matrix(d, rho) else base_cor_mat
      z_scores <- as.numeric(rmvnorm(1, mean = true_effect, sigma = cor_mat))
      p_values <- 2 * (1 - pnorm(abs(z_scores)))
      
      gbj_result <- tryCatch(GBJ::GBJ(test_stats = z_scores, cor_mat = cor_mat)$GBJ_pvalue, error = function(e) NA)
      ghc_result <- tryCatch(GBJ::GHC(test_stats = z_scores, cor_mat = cor_mat)$GHC_pvalue, error = function(e) NA)
      wartp <- tryCatch(wARTP(p_values, cor_mat = cor_mat, weight = rep(1, d))$p_value, error = function(e) NA)
      
      return(c(wARTP = wartp, GBJ = gbj_result, GHC = ghc_result))
    }, error = function(e) {
      message(sprintf("Simulation %d failed: %s", i, e$message))
      return(rep(NA, 3))
    })
  }, mc.cores = detectCores() - 1)
  
  pv <- do.call(rbind, sim_results)
  colnames(pv) <- c("wARTP", "GBJ", "GHC")
  
  df <- as.data.frame(-log(pv))
  df$sim <- 1:nrow(df)
  df_long <- melt(df, id.vars = "sim", variable.name = "method", value.name = "neglogp")
  df_long$method <- factor(df_long$method, levels = c("wARTP", "GBJ", "GHC"))
  
  print(
    ggplot(df_long, aes(x = method, y = neglogp, fill = method)) +
      geom_boxplot(width = 0.4, outlier.shape = 8, outlier.size = 2) +
      geom_hline(yintercept = -log(0.05), color = "red", linetype = "dashed") +
      labs(y = "-log(p-value)", x = "", title = paste("Multivariate Normal Simulation (cortype =", cortype, ")")) +
      theme_minimal() + theme(legend.position = "none")
  )
  
  print(apply(pv < 0.05, 2, sum))
  print(apply(pv, 2, median))
  return(pv)
}
