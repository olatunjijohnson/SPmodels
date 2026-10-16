\name{sp_variogram_envelope}
\alias{sp_variogram_envelope}
\alias{print.spmodels_variogram_envelope}
\alias{plot.spmodels_variogram_envelope}
\title{Monte Carlo Envelope for a Residual Variogram}
\usage{
sp_variogram_envelope(
  object,
  nsim = 999L,
  method = c("permutation", "parametric"),
  level = 0.95,
  bins = 10L,
  max_distance = NULL,
  residual_type = NULL,
  seed = NULL
)

\method{print}{spmodels_variogram_envelope}(x, ...)

\method{plot}{spmodels_variogram_envelope}(x, ...)
}
\arguments{
\item{object}{A fitted spatial \code{spmodels_fit}.}
\item{nsim}{Positive number of Monte Carlo curves.}
\item{method}{Permutation null or parametric simulations from the fitted model.}
\item{level}{Pointwise envelope probability strictly between zero and one.}
\item{bins}{Positive number of equal-width distance bins.}
\item{max_distance}{Optional positive largest pair distance.}
\item{residual_type}{Passed to \code{sp_variogram()}: Gaussian mean or
conditional residuals, or binomial/Poisson randomized-quantile, Pearson or
deviance residuals.}
\item{seed}{Optional integer seed local to this call.}
\item{x}{A variogram-envelope object.}
\item{...}{Additional graphical arguments.}
}
\value{An object containing observed and simulated variogram curves, pointwise
limits, standardized global statistics and a Monte Carlo p-value. Discrete
parametric envelopes additionally retain every refit's convergence code and
captured warning messages.}
\description{Builds a pointwise residual-variogram envelope and a global
curve-discrepancy test from residual permutations or fitted-model parametric
simulations for any supported spatial family.}
\details{The permutation null is exchangeability of the selected residual
values across fixed locations. Gaussian parametric curves keep fitted
parameters fixed. Binomial and Poisson parametric curves simulate a response,
repeat the declared Laplace or Monte Carlo ML fitting route, and calculate the
selected residuals from that refit. This includes fitting variability but is
more computationally expensive.}
