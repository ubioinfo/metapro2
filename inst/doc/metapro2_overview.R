## -----------------------------------------------------------------------------
library(metapro2)
set.seed(123)
pvals <- runif(50) / 5
res <- wARTP(pvals)
print(res)

