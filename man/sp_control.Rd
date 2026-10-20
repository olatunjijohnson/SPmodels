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
  distance_units = NULL,
  maxit = 500L,
  reltol = 1e-8
)
}
\arguments{
  \item{smoothness}{Fixed positive Matern smoothness.}
  \item{spatial_variance, scale, nugget_variance}{Optional positive starting
  values during estimation. Fixed spatial fitting requires all three when
  \code{estimate = FALSE}; fixed non-spatial fitting requires only
  \code{nugget_variance}.}
  \item{random_variance}{Optional positive scalar or vector of starting or
  fixed Gaussian random-intercept variances. Named values must match the
  grouping variables in \code{random}.}
  \item{estimate}{Whether to estimate covariance and declared random-effect
  variance parameters by maximum likelihood.}
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
non-Gaussian covariance settings.}
