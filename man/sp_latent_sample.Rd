\name{sp_latent_sample}
\alias{sp_latent_sample}
\title{Sample a Non-Gaussian Latent Spatial Vector}
\usage{
sp_latent_sample(
  y,
  linear_predictor,
  covariance,
  family = stats::binomial(),
  trials = NULL,
  initial = NULL,
  control = sp_mcmc_control(),
  seed = NULL
)
}
\arguments{
  \item{y}{Numeric binomial successes or Poisson counts.}
  \item{linear_predictor}{Finite fixed linear predictor before the latent vector.}
  \item{covariance}{Finite symmetric positive-definite latent covariance matrix.}
  \item{family}{A binomial-logit or Poisson-log family.}
  \item{trials}{Binomial integer trial counts; defaults to one per outcome.}
  \item{initial}{Optional finite initial latent vector, or a matrix with one
  row per chain and one column per latent component.}
  \item{control}{Settings from \code{sp_mcmc_control()}.}
  \item{seed}{Optional integer seed local to this call.}
}
\value{A \code{spmodels_latent_sample} object containing pooled draws, a
draw-by-component-by-chain array, chain identities, chain-specific acceptance
information, initial states, the latent mode and resolved controls.}
\description{Runs a preconditioned Metropolis-adjusted Langevin algorithm for
the supplied latent posterior. This samples latent effects conditional on
supplied parameters; it is not by itself a fitted parameter model.}
