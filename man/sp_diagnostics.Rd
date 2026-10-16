\name{sp_diagnostics}
\alias{sp_diagnostics}
\alias{plot.spmodels_diagnostics}
\title{Fitted Spatial-Model Diagnostics}
\usage{
sp_diagnostics(object, seed = NULL)

\method{plot}{spmodels_diagnostics}(x, ...)
}
\arguments{
  \item{object}{A fitted \code{spmodels_fit} object.}
  \item{seed}{Optional integer seed for non-Gaussian randomized quantile
  residuals. The caller's random-number state is restored.}
  \item{x}{An object returned by \code{sp_diagnostics()}.}
  \item{...}{Additional graphical arguments for the residual-index panel.}
}
\value{For Gaussian fits, \code{sp_diagnostics()} returns mean-model,
conditional and decorrelated residuals. For binomial and Poisson fits it
returns response, Pearson, signed deviance and randomized quantile residuals at
the fitted latent mode. Coordinates, distance units and family are retained. The plot method
invisibly returns its input.}
\description{Provides family-appropriate residual quantities for checking a
fitted spatial model. Gaussian decorrelated residuals and non-Gaussian Pearson,
deviance or randomized quantile residuals have different interpretations and are labelled
separately.}
