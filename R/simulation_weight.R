# Simulation benchmark comparing weighted vs unweighted p-value combination methods under heteroscedastic sample sizes (weights)

simul_weight <- function(N = 15, d = 50, r = 10, e = 0.2) {
  library(GBJ)
  library(ggplot2)
  library(reshape2)
  library(parallel)
  library(mvtnorm)
  
  if (r > d) stop("r (number of signals) must be <= d")
  
  sim_results <- function(i) {
    n1 <- sample(20:200, d, replace = TRUE)
    mu <- c(rep(e, r), rep(0, d - r))
    
    p_values <- sapply(1:d, function(j) {
      g1 <- rnorm(n1[j], mean = 0, sd = 1)
      g2 <- rnorm(n1[j], mean = mu[j], sd = 1)
      t.test(g1, g2, var.equal = TRUE)$p.value
    })
    
    weights <- n1
    
    wartp <- tryCatch(wARTP(p_values, cor_mat = diag(d), weight = weights)$p_value, error = function(e) NA)
    artp <- tryCatch(wARTP(p_values, cor_mat = diag(d))$p_value, error = function(e) NA)
    wf_w <- tryCatch(wFisher(p_values, weight = weights)$p, error = function(e) NA)
    wf_u <- tryCatch(wFisher(p_values)$p, error = function(e) NA)
    minp <- tryCatch(minP(test_stats = qnorm(1 - p_values / 2), cor_mat = diag(d))$minP_pvalue, error = function(e) NA)
    gbj <- tryCatch(GBJ(test_stats = qnorm(1 - p_values / 2), cor_mat = diag(d))$GBJ_pvalue, error = function(e) NA)
    ghc <- tryCatch(GHC(test_stats = qnorm(1 - p_values / 2), cor_mat = diag(d))$GHC_pvalue, error = function(e) NA)
    ord <- tryCatch(ordmeta(p_values)$p_value, error = function(e) NA)
    
    return(c(wARTP = wartp, ARTP = artp,
             wFisher = wf_w, Fisher = wf_u, ordmeta = ord, GBJ = gbj, GHC = ghc, minP = minp))
  }
  
  cat("Running simulations in parallel...\n")
  results <- mclapply(1:N, function(i) {
    tryCatch(sim_results(i), error = function(e) {
      message(sprintf("Error in sim %d: %s", i, e$message))
      return(rep(NA, 8))
    })
  }, mc.cores = detectCores() - 1)
  
  pv <- do.call(rbind, results)
  colnames(pv) <- c("wARTP", "ARTP",
                    "wFisher", "Fisher", "ordmeta", "GBJ", "GHC", "minP")
  
  df <- as.data.frame(-log(pv))
  df$sim <- 1:nrow(df)
  df_long <- reshape2::melt(df, id.vars = "sim", variable.name = "method", value.name = "neglogp")
  
  df_long$method <- factor(df_long$method,
                           levels = c("wARTP", "ARTP",
                                      "wFisher", "Fisher", "ordmeta", "GBJ", "GHC", "minP"))
  
  print(
    ggplot(df_long, aes(x = method, y = neglogp, fill = method)) +
      geom_boxplot(width = 0.4, outlier.shape = 8, outlier.size = 2) +
      geom_hline(yintercept = -log(0.05), color = "red", linetype = "dashed") +
      labs(y = "-log(p-value)", x = "", title = "Weighted vs Unweighted P-value Combination Comparison") +
      theme_minimal() + theme(legend.position = "none")
  )
  
  print(apply(pv < 0.05, 2, sum, na.rm = TRUE))
  print(apply(pv, 2, median, na.rm = TRUE))
  return(pv)
}
