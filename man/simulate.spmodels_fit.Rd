\name{simulate.spmodels_fit}
\alias{simulate.spmodels_fit}
\title{Simulate from a Fitted Spatial Model}
\usage{
\method{simulate}{spmodels_fit}(object, nsim = 1L, seed = NULL, ...)
}
\arguments{
  \item{object}{A fitted \code{spmodels_fit} object.}
  \item{nsim}{Positive number of response simulations.}
  \item{seed}{Optional integer seed. The caller's random-number state is
  restored after seeded simulation.}
  \item{...}{Reserved for method compatibility.}
}
\value{A data frame with one simulated response per column.}
\description{Generates unconditional Gaussian, binomial or Poisson response
replicates at the fitted observation sites using the fitted mean, family and
covariance parameters. Binomial output contains simulated success counts with
the fitted trial totals.}
