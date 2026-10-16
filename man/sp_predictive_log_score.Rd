\name{sp_predictive_log_score}
\alias{sp_predictive_log_score}
\title{Integrated Pointwise Predictive Log Scores}
\usage{
sp_predictive_log_score(
  object,
  newdata,
  observed,
  trials = NULL,
  unseen = c("error", "marginal")
)
}
\arguments{
  \item{object}{A fitted \code{spmodels_fit} model.}
  \item{newdata}{Prediction data containing all required model fields.}
  \item{observed}{One finite observed response per row. Binomial values are
  success counts.}
  \item{trials}{Optional non-negative integer binomial trial totals. Defaults
  to one for binomial models and must be omitted otherwise.}
  \item{unseen}{Policy for grouping levels absent from the fitted data.}
}
\value{A numeric vector of negative log predictive densities.}
\description{Evaluates Gaussian predictive normal scores or integrates a
binomial/Poisson response probability over the approximate Gaussian predictive
distribution on the link scale.}
\details{The non-Gaussian integral is one-dimensional for each prediction row
and conditional on fitted parameters. Binomial intervention multipliers are
included. A failed or numerically underflowed integral produces an explicit
error rather than a clipped score.}
