\name{sp_variogram}
\alias{sp_variogram}
\alias{plot.spmodels_variogram}
\title{Empirical Residual Variograms Before and After Spatial Fitting}
\usage{
sp_variogram(
  object = NULL,
  bins = 10L,
  max_distance = NULL,
  residual_type = NULL,
  seed = NULL,
  coords = NULL,
  values = NULL,
  data = NULL,
  distance_units = NULL
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
  \item{coords}{For an ordinary \code{lm}/\code{glm} or direct residual
  vector, either two numeric coordinate columns or a one-sided coordinate
  formula.}
  \item{values}{Numeric residual values when \code{object = NULL}.}
  \item{data}{Data used to evaluate a coordinate formula for an ordinary
  model. Supply the complete source data when coordinates are not in the
  model formula.}
  \item{distance_units}{Optional target unit for distances. Spatial and
  non-spatial SPmodels fits are converted from their stored units; projected
  \code{sf} coordinates are converted when their CRS units are recognised.
  Numeric coordinate columns are assumed already to use this unit.}
  \item{x}{An object returned by \code{sp_variogram()}.}
  \item{...}{Additional graphical arguments passed to \code{plot}.}
}
\value{\code{sp_variogram()} returns an object of class
\code{spmodels_variogram} containing binned empirical semivariances and, for
Gaussian fits, the fitted model semivariance. The plot method invisibly
returns its input.}
\description{Calculates binned pairwise residual semivariances with pair
counts. Before fitting, ordinary residuals from \code{lm}, \code{glm}, or a
non-spatial SPmodels fit can be used to assess whether spatial structure
motivates fitting a spatial process. After fitting, conditional residuals can
be used to assess remaining spatial structure. These are distinct diagnostic
questions. Gaussian spatial-fit results can include the fitted
covariance-model variogram. Discrete-family residual variograms deliberately
omit a fitted response-scale curve rather than treating link covariance as
residual semivariance.}
