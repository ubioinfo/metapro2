P-value combination is commonly used in a variety of genomic studies such as meta-analysis, multiomics integration, and set-based association analysis, with or without weights and correlation structures. 
In most cases, the signal sparsity among the genomic entities or datasets is not known, where powerful adaptive p-value combination is required.

To this end, metapro2 package implements a weighted version of the empirical test, adaptive rank truncated product (wARTP) in Rcpp, as well as weighted Fisher and Monte Carlo-based ordmeta methods. 
wARTP and ordmeta also provide estimates of signal sparsity among the combined p-values.

Original ARTP method is a powerful empirical p-value combination method for both indenepdent and dependent p-values, but its long computation time is a major drawback. 
Thus, we implement simulation-based ARTP in Rcpp which is an order of magnitude faster than ARTP R code. For independent p-values, sample-size weighting is tested for both wARTP and wFisher. 
wARTP exhibits the highest power for sparse signals, while wFisher performs the best for dense signals. wARTP is a versatile algorithm that exhibit high power for various sparsity and correlation levels.


# metapro2

The `metapro2` R package provides efficient methods for combining p-values in omics studies, particularly when signal sparsity or correlation structures are present. It includes a fast Rcpp implementation of:

- `wARTP`: weighted adaptive rank truncated product method  
- `wFisher`: weighted Fisher’s method  
- `ordmeta`: signal sparsity estimation  

These are useful in meta-analysis, gene-set testing, and multi-omics integration.

---

## Installation
You can install **metapro2** from GitHub:

```r
# install.packages("devtools")  # if not installed yet
devtools::install_github("ubioinfo/metapro2")

```

## Running wARTP
```r
library(metapro2)

# Generate 10 random p-values
set.seed(123)
p <- runif(10)

# Identity correlation matrix
cor_mat <- diag(10)

# Run wARTP with 50,000 permutations (default)
res <- wARTP(p, cor_mat = cor_mat, B2 = 50000) #the same as wARTP(p) in this case.

# View result
res$p_value
res$signal_sparsity  # Estimated proportion of informative p-values

```
 

## Running wFisher and ordmeta
```r
library(metapro2)

# Simulated p-values and weights
p <- c(0.01, 0.20, 0.05, 0.02)
weights <- c(1, 2, 1, 1)

# Weighted Fisher method
res <- wFisher(p, weights)
res$p_value

# Estimate sparsity using ordmeta
res <- ordmeta(p)
res$p_value
res$opt_ord  # Estimated number of informative p-values
```

## ✅ Platform-specific Installation Notes

### macOS & Linux

For most macOS and Linux users, installing from GitHub should work smoothly if the following are met:

- **Xcode Command Line Tools** installed (macOS only):

```sh
xcode-select --install
```

- **A working C++ compiler** (e.g., `g++`, `clang++`) with support for **C++17**

- **Internet access** for downloading dependencies

If you encounter errors related to `RcppArmadillo.h` or `LAPACK/BLAS`, try reinstalling the packages from source:

```r
install.packages("RcppArmadillo", type = "source")
install.packages("Rcpp", type = "source")
```

You may also define a `~/.R/Makevars` file with the following content:

```make
CXXFLAGS = -O2 -Wall -pedantic
CXX11FLAGS = -std=c++11
```

---

### ✅ Windows

Installing from GitHub on Windows may require some additional configuration:

- Install **Rtools**, and verify it’s properly configured:

```r
Sys.which("make")
```

- Install dependencies from source:

```r
install.packages("RcppArmadillo", type = "source")
install.packages("Rcpp", type = "source")
```

- Create or edit the file below:

```
C:/Users/YourName/Documents/.R/Makevars.win
```

With the following content (adjust paths as needed):

```make
CXXFLAGS = -I"C:/Users/YourName/AppData/Local/R/win-library/4.5/Rcpp/include" -I"C:/Users/YourName/AppData/Local/R/win-library/4.5/RcppArmadillo/include"
CXX14FLAGS = $(CXXFLAGS)
PKG_LIBS = -llapack -lblas -lgfortran -lquadmath
```

If you encounter errors like `undefined reference to 'dpotrf_'`, `quadmath_snprintf`, or `RcppArmadillo.h not found`, the above settings usually resolve them.

---

## 🚀 Installation (All Platforms)

### Option 1: Install from GitHub (latest version)

```r
devtools::install_github("ubioinfo/metapro2")
```

### Option 2: Install from release source tarball

```r
install.packages("https://github.com/ubioinfo/metapro2/releases/download/v0.1.0/metapro2_0.1.0.tar.gz",
                 repos = NULL, type = "source")
```
