\name{sp_mcml_control}
\alias{sp_mcml_control}
\alias{print.spmodels_mcml_control}
\title{Control Monte Carlo Maximum-Likelihood Correction}
\usage{
sp_mcml_control(
  nsim = 2000L,
  maxit = 200L,
  reltol = 1e-06,
  proposal_scale = 1,
  ess_warning = 0.1,
  adapt = FALSE,
  pilot_nsim = 500L,
  adapt_shrinkage = 0.5
)
}
\arguments{
  \item{nsim}{Number of independent importance samples; at least 100.}
  \item{maxit}{Positive maximum number of optimization iterations.}
  \item{reltol}{Positive relative optimization tolerance.}
  \item{proposal_scale}{One or more distinct positive standard-deviation
  multipliers for the local Gaussian latent-posterior proposal. Multiple
  values define a normalized deterministic mixture with balanced allocation.}
  \item{ess_warning}{Importance-sampling effective-sample fraction strictly
  between zero and one below which \code{sp_mcml()} warns. When adaptation is
  requested, a pilot below this threshold retains the original proposal
  moments.}
  \item{adapt}{Whether to use an independent discarded pilot sample to
  estimate regularized proposal moments.}
  \item{pilot_nsim}{Number of pilot importance draws; at least 100.}
  \item{adapt_shrinkage}{Number greater than zero and no greater than one.
  It is the weight assigned to the original Laplace moments; the remaining
  weight is assigned to the pilot importance-weighted moments.}
}
\value{An object of class \code{spmodels_mcml_control}.}
\description{Configures the normalized importance distribution, optional
discarded-pilot adaptation and deterministic common-random-number optimization
used by \code{sp_mcml()}.}
