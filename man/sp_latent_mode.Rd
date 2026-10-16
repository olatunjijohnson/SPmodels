\name{sp_latent_mode}
\alias{sp_latent_mode}
\title{Calculate a Latent Gaussian Mode}
\usage{
sp_latent_mode(
  y,
  linear_predictor,
  covariance,
  family = stats::binomial(),
  trials = NULL,
  probability_multiplier = NULL,
  initial = NULL,
  maxit = 100L,
  tolerance = 1e-8
)
}
\arguments{
  \item{y}{Numeric outcomes: binomial successes or Poisson counts.}
  \item{linear_predictor}{Finite fixed linear predictor before the latent vector.}
  \item{covariance}{Finite symmetric positive-definite latent covariance matrix.}
  \item{family}{A binomial-logit or Poisson-log family.}
  \item{trials}{Binomial integer trial counts; defaults to one per outcome.}
  \item{probability_multiplier}{Optional known binomial prevalence multiplier
  in \code{(0, 1]}; not used for Poisson outcomes.}
  \item{initial}{Optional finite initial latent vector.}
  \item{maxit}{Positive maximum Newton iterations.}
  \item{tolerance}{Positive convergence tolerance for the maximum step.}
}
\value{A \code{spmodels_latent_mode} list containing the mode, inverse local
curvature, gradient, iteration information, conditional means and log posterior.}
\description{Finds the posterior mode of a latent Gaussian vector for
independent binomial-logit or Poisson-log observations. This is a deterministic
building block for later Monte Carlo maximum-likelihood fitting and does not
itself estimate model parameters.}
