\name{sp_predict_target}
\alias{sp_predict_target}
\alias{print.spmodels_prediction_target}
\title{Summarise Targets from Joint Predictive Draws}
\usage{
sp_predict_target(
  object,
  weights = NULL,
  groups = NULL,
  threshold = NULL,
  probs = c(0.025, 0.5, 0.975)
)

\method{print}{spmodels_prediction_target}(x, ...)
}
\arguments{
  \item{object}{An object returned by \code{sp_predict_draws()}.}
  \item{weights}{Optional finite non-negative weight for every prediction
  location, or an object returned by \code{sp_polygon_weights()}. Weights are
  normalised within each group.}
  \item{groups}{Optional non-missing group identifier for every location. It
  requires \code{weights}; the default is one aggregate target.}
  \item{threshold}{Optional finite scalar used to calculate exceedance
  probabilities.}
  \item{probs}{Finite quantile probabilities between zero and one.}
  \item{x}{An object returned by \code{sp_predict_target()}.}
  \item{...}{Reserved for print-method compatibility.}
}
\value{An object of class \code{spmodels_prediction_target} containing
target-level draws, a summary table, normalised weights, groups, threshold and
scale metadata.}
\description{Calculates means, standard deviations, quantiles and optional
exceedance probabilities from joint predictive draws. Optional grouped
weighted means preserve uncertainty because aggregation is performed within
every joint draw before summarisation.}
