\name{sp_control}
\alias{sp_control}
\title{Control Gaussian Spatial Model Fitting}
\usage{
sp_control(
  smoothness = 0.5,
  spatial_variance = NULL,
  scale = NULL,
  nugget_variance = NULL,
  random_variance = NULL,
  estimate = TRUE,
  estimate_nugget = NULL,
  distance_units = NULL,
  maxit = 500L,
  reltol = 1e-8
)
}
\arguments{
  \item{smoothness}{Fixed positive Matern smoothness.}
  \item{spatial_variance, scale}{Optional positive starting values when
  \code{estimate = TRUE}. Fixed spatial fitting requires both.}
  \item{nugget_variance}{Optional non-negative nugget variance. It is a
  starting value when estimated and a fixed value when
  \code{estimate_nugget = FALSE}. Supplying zero fixes the nugget at zero
  unless \code{estimate_nugget = TRUE}, which requires a positive start.}
  \item{random_variance}{Optional positive scalar or vector of starting or
  fixed Gaussian random-intercept variances. Named values must match the
  grouping variables in \code{random}.}
  \item{estimate}{Whether to estimate the spatial variance, scale, and
  supported additional variance components by maximum likelihood.}
  \item{estimate_nugget}{Optional logical value controlling nugget estimation
  when \code{estimate = TRUE}. The default estimates a missing or positive
  nugget and fixes a supplied zero nugget. Set it to \code{FALSE} to fix any
  supplied non-negative nugget value.}
  \item{distance_units}{A non-empty description of the units shared by both
  coordinate axes and used for model distances and covariance scales, such as
  \code{"km"}. Numeric coordinates are assumed to use these units already.
  Projected \code{sf} coordinates with a recognised linear CRS unit are
  converted before fitting.}
  \item{maxit}{Positive maximum number of optimisation iterations.}
  \item{reltol}{Positive optimisation tolerance.}
}
\value{An object of class \code{spmodels_control}.}
\description{Creates controls for Gaussian geostatistical fitting and shared
non-Gaussian covariance settings. The nugget may be fixed independently of
the other covariance parameters, including at exactly zero.}
