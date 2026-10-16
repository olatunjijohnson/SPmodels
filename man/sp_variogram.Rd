\name{sp_variogram}
\alias{sp_variogram}
\alias{plot.spmodels_variogram}
\title{Empirical Residual Variogram from a Fitted Spatial Model}
\usage{
sp_variogram(
  object,
  bins = 10L,
  max_distance = NULL,
  residual_type = NULL,
  seed = NULL
)

\method{plot}{spmodels_variogram}(x, ...)
}
\arguments{
  \item{object}{A fitted \code{spmodels_fit} object.}
  \item{bins}{Positive number of equal-width distance bins.}
  \item{max_distance}{Optional positive maximum pair distance.}
  \item{residual_type}{For Gaussian fits, use \code{"mean"} or
  \code{"conditional"} residuals. For binomial/Poisson fits, use
  \code{"randomized_quantile"}, \code{"pearson"} or \code{"deviance"}.
  The family-specific first choice is the default.}
  \item{seed}{Optional integer seed local to randomized-quantile residuals.}
  \item{x}{An object returned by \code{sp_variogram()}.}
  \item{...}{Additional graphical arguments passed to \code{plot}.}
}
\value{\code{sp_variogram()} returns an object of class
\code{spmodels_variogram} containing binned empirical semivariances and, for
Gaussian fits, the fitted model semivariance. The plot method invisibly
returns its input.}
\description{Calculates binned pairwise residual semivariances with pair
counts. Gaussian results include the fitted covariance-model variogram.
Discrete-family residual variograms deliberately omit a fitted response-scale
curve rather than treating link covariance as residual semivariance.}
