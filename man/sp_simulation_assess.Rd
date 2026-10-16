\name{sp_simulation_assess}
\alias{sp_simulation_assess}
\alias{print.spmodels_simulation_assessment}
\title{Assess Predictive Draws Against Known Values}
\usage{
sp_simulation_assess(object, truth, interval = 0.95, threshold = NULL)

\method{print}{spmodels_simulation_assessment}(x, ...)
}
\arguments{
  \item{object}{An object returned by \code{sp_predict_draws()}.}
  \item{truth}{One finite value for every prediction location.}
  \item{interval}{Central interval probability strictly between zero and one.}
  \item{threshold}{Optional finite threshold for predicted exceedance
  probabilities and observed threshold indicators.}
  \item{x}{An object returned by \code{sp_simulation_assess()}.}
  \item{...}{Reserved for print-method compatibility.}
}
\value{An object of class \code{spmodels_simulation_assessment} containing
location-level errors, intervals, coverage and empirical CRPS, plus overall
bias, RMSE, MAE, coverage, interval width and mean CRPS.}
\description{Assesses a joint predictive draw ensemble against known simulated
truth or held-out observations. All interval and score definitions are
calculated directly from the retained draws.}
