
\name{wFisher}
\alias{wFisher}
\title{Weighted Fisher's Method for Combining P-values}
\description{Combines p-values using a weighted version of Fisher's method with Gamma transformation.}
\usage{
wFisher(p, weight = NULL, is.onetail = TRUE)
}
\arguments{
  \item{p}{A numeric vector of p-values}
  \item{weight}{A numeric vector of weights (optional). Default = equal weights}
}
\value{
A list with:
\itemize{
  \item \code{p} Combined p-value
}
}
\examples{
set.seed(42)
p <- runif(10) / 2
wFisher(p)
}
