\name{sp_mcmc_control}
\alias{sp_mcmc_control}
\title{Control Latent MCMC Sampling}
\usage{
sp_mcmc_control(
  iterations = 4000L,
  burnin = 1000L,
  thin = 3L,
  step_size = 0.7,
  chains = 1L,
  initial_scale = 2
)
}
\arguments{
  \item{iterations}{Total positive number of MALA iterations.}
  \item{burnin}{Non-negative initial iterations discarded.}
  \item{thin}{Positive retention interval after burn-in.}
  \item{step_size}{Positive MALA step size.}
  \item{chains}{Positive number of independently evolved chains.}
  \item{initial_scale}{Non-negative multiplier for dispersed initial states
  around the latent mode when multiple chains are requested.}
}
\value{An object of class \code{spmodels_mcmc_control}.}
\description{Defines explicit iteration and proposal settings for the
preconditioned Metropolis-adjusted Langevin sampler.}
