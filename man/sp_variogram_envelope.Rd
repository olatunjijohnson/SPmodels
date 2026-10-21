\name{sp_variogram_envelope}
\alias{sp_variogram_envelope}
\alias{print.spmodels_variogram_envelope}
\alias{plot.spmodels_variogram_envelope}
\title{Monte Carlo Envelope for a Residual Variogram}
\usage{
sp_variogram_envelope(
  object = NULL,
  nsim = 999L,
  method = c("permutation", "parametric"),
  level = 0.95,
  bins = 10L,
  max_distance = NULL,
  residual_type = NULL,
  seed = NULL,
  coords = NULL,
  values = NULL,
  data = NULL,
  distance_units = NULL
)

\method{print}{spmodels_variogram_envelope}(x, ...)

\method{plot}{spmodels_variogram_envelope}(x, ...)
}
\arguments{
\item{object}{A spatial or non-spatial \code{spmodels_fit}, an
\code{lm}/\code{glm} object, or \code{NULL} when supplying \code{coords}
and \code{values}. The permutation method supports pre-fit inputs; the
parametric method requires a fitted spatial model.}
\item{nsim}{Positive number of Monte Carlo curves.}
\item{method}{Permute residuals across fixed locations, or simulate from the
fitted model while keeping fitted parameters fixed.}
\item{level}{Pointwise envelope probability strictly between zero and one.}
\item{bins}{Positive number of equal-width distance bins.}
\item{max_distance}{Optional positive largest pair distance.}
\item{residual_type}{Passed to \code{sp_variogram()}: Gaussian mean or
conditional residuals, or binomial/Poisson randomized-quantile, Pearson or
deviance residuals.}
\item{seed}{Optional integer seed local to this call.}
\item{coords}{Two numeric coordinate columns, a one-sided coordinate
formula, or projected \code{sf} data, passed to \code{sp_variogram()}.}
\item{values}{Numeric residual values when \code{object = NULL}.}
\item{data}{Data used to evaluate a coordinate formula for an ordinary
model.}
\item{distance_units}{Optional target unit for distances. Stored model or
projected \code{sf} coordinates are converted when their source units are
known; numeric coordinate values are assumed already in this unit.}
\item{x}{A variogram-envelope object.}
\item{...}{Additional graphical arguments.}
}
\value{An object containing observed and simulated variogram curves, pointwise
limits, standardized global statistics and a Monte Carlo p-value. Discrete
parametric envelopes additionally retain every refit's convergence code and
captured warning messages.}
\description{Builds a pointwise residual-variogram envelope and a global
curve-discrepancy test from residual permutations or fitted-model parametric
simulations. Pre-fit permutation envelopes are supported for ordinary model
residuals, supplied residuals and non-spatial SPmodels fits.}
\details{For a pre-fit diagnostic, the permutation null is exchangeability of
the selected ordinary-model residual values across fixed locations. This
provides the week 7 teaching comparison before fitting a spatial process.
Gaussian parametric curves keep fitted
parameters fixed. Binomial and Poisson parametric curves simulate a response,
repeat the declared Laplace or Monte Carlo ML fitting route, and calculate the
selected residuals from that refit. This includes fitting variability but is
more computationally expensive.}
