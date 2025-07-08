P-value combination is commonly used in a variety of genomic studies such as meta-analysis, multiomics integration, and set-based association analysis, with or without weights and correlation structures. 
In most cases, the signal sparsity among the genomic entities or datasets is not known, where powerful adaptive p-value combination is required.

To this end, metapro2 package implements a weighted version of the empirical test, adaptive rank truncated product (wARTP) in Rcpp, as well as weighted Fisher and Monte Carlo-based ordmeta methods. 
wARTP and ordmeta also provide estimates of signal sparsity among the combined p-values.

Original ARTP method is a powerful empirical p-value combination method for both indenepdent and dependent p-values, but its long computation time is a major drawback. 
Thus, we implement simulation-based ARTP in Rcpp which is an order of magnitude faster than ARTP R code. For independent p-values, sample-size weighting is tested for both wARTP and wFisher. 
wARTP exhibits the highest power for sparse signals, while wFisher performs the best for dense signals. wARTP is a versatile algorithm that exhibit high power for various sparsity and correlation levels.


# metapro2

The `metapro2` R package provides efficient methods for combining p-values in genomic studies, particularly when signal sparsity or correlation structures are present. It includes a fast Rcpp implementation of adaptive rank truncated product methods (`wARTP`), weighted Fisher's method (`wFisher`), and signal sparsity estimation (`ordmeta`). These methods are useful in meta-analysis, gene-set testing, and multi-omics integration.

---

## Installation

#install.packages("devtools")  # if not already installed
devtools::install_github("ubioinfo/metapro2")

library(metapro2)


# Generate 10 random p-values
set.seed(123)
p <- runif(10)

# No correlation assumed (identity matrix)
cor_mat <- diag(10)

# Run wARTP with 10,000 permutations
res <- wARTP(p, cor_mat = cor_mat, B2 = 10000)

# View p-value
res$p_value
res$signal_sparsity     # The optimal number of p-values divided by the total number of p-values


# Simulated p-values and weights
p <- c(0.01, 0.20, 0.05, 0.02)
weights <- c(1, 2, 1, 1)

# wFisher and ordmeta
res <- wFisher(p, weights)
res$p_value

res <- ordmeta(p)
res$p_value
res$opt_ord     # Optimal order indicating the number of core p-values, similar as signal_sparsity in wARTP
