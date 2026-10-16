\name{sp_environmental_folds}
\alias{sp_environmental_folds}
\alias{print.spmodels_environmental_folds}
\title{Construct Environmental-Space Validation Folds}
\usage{
sp_environmental_folds(
  object,
  variables,
  folds = 5L,
  scale = TRUE,
  seed = NULL,
  nstart = 20L
)

\method{print}{spmodels_environmental_folds}(x, ...)
}
\arguments{
  \item{object}{A fitted spatial \code{spmodels_fit}.}
  \item{variables}{Unique character column names or a one-sided additive
  formula selecting numeric environmental variables from fitted data.}
  \item{folds}{Number of environmental clusters.}
  \item{scale}{Standardize variables before clustering?}
  \item{seed}{Optional integer seed local to k-means.}
  \item{nstart}{Positive number of random k-means starts.}
  \item{x}{An object returned by \code{sp_environmental_folds()}.}
  \item{...}{Reserved for print-method compatibility.}
}
\value{An integer fold vector in fitted-row order with class
\code{spmodels_environmental_folds}. Environmental centres, variables,
scaling values and seed are retained as attributes.}
\description{Clusters unique fitted locations in a declared numeric
environmental space. All rows at one coordinate remain together and must have
identical environmental values, preventing exact-location leakage.}
\details{Variables are standardized by default so measurement units do not
arbitrarily determine Euclidean k-means distance. The resulting folds test
transfer across environmental regimes; they are not spatial distance buffers
and need not be geographically compact.}
