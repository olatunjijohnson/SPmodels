\name{sp_calibration}
\alias{sp_calibration}
\alias{print.spmodels_calibration}
\alias{plot.spmodels_calibration}
\title{Assess and display predictive interval calibration}
\usage{
sp_calibration(
  object,
  truth,
  levels = c(0.5, 0.8, 0.9, 0.95),
  confidence = 0.95
)

\method{print}{spmodels_calibration}(x, ...)

\method{plot}{spmodels_calibration}(
  x,
  xlab = "Nominal coverage",
  ylab = "Observed coverage",
  main = "Predictive interval calibration",
  ...
)
}
\arguments{
\item{object}{An object returned by \code{sp_predict_draws()}.}
\item{truth}{One finite value for every prediction location.}
\item{levels}{Central predictive-interval probabilities strictly between zero and one.}
\item{confidence}{Confidence level for pointwise Wilson bounds around empirical coverage.}
\item{x}{An object returned by \code{sp_calibration()}.}
\item{xlab, ylab, main}{Plot labels.}
\item{...}{Further arguments passed to the principal base-graphics call.}
}
\value{
\code{sp_calibration()} returns a \code{spmodels_calibration} object containing
the unrounded nominal and observed coverage values, Wilson bounds and the
number of locations. The print and plot methods return \code{x} invisibly.
}
\description{
Forms central type-8 empirical predictive intervals at each requested nominal
level, evaluates their coverage of the supplied truths and calculates
pointwise Wilson uncertainty bounds for the observed proportions.
}
\details{
The bounds quantify finite-location uncertainty in the empirical coverage
proportion. They are not simultaneous bands and do not incorporate dependence
between spatial locations.
}
