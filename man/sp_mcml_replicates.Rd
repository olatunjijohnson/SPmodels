\name{sp_mcml_replicates}
\alias{sp_mcml_replicates}
\alias{print.spmodels_mcml_replicates}
\title{Repeat Monte Carlo maximum-likelihood fits independently}
\usage{
sp_mcml_replicates(
  object,
  nrep = 3L,
  control = sp_mcml_control(),
  seed = NULL,
  verbose = FALSE
)

\method{print}{spmodels_mcml_replicates}(x, ...)
}
\arguments{
\item{object}{A fitted binomial-logit or Poisson-log \code{spmodels_fit}.}
\item{nrep}{Number of independent fits, at least two.}
\item{control}{Settings returned by \code{sp_mcml_control()}.}
\item{seed}{Optional integer defining consecutive, call-local replicate seeds.}
\item{verbose}{Logical scalar passed to \code{sp_mcml()}.}
\item{x}{An object returned by \code{sp_mcml_replicates()}.}
\item{...}{Unused.}
}
\value{
An object of class \code{spmodels_mcml_replicates} containing every fit,
replicate estimates, between-replicate summaries and likelihood diagnostics.
It also includes empirical Monte Carlo, mean fitted-Hessian sampling, and
combined parameter covariance matrices on both optimizer and natural model
scales. When adaptation is requested, diagnostics also contain pilot ESS,
maximum weight, relative MCSE and the used/fallback decision.
}
\description{
Runs \code{sp_mcml()} with independent importance samples. Replication
diagnoses Monte Carlo sensitivity; it does not replace scientific sampling
uncertainty or combine the fits into a new estimator.
The empirical between-replicate covariance measures parameter-level Monte
Carlo sensitivity. When every fitted Hessian covariance is usable, the object
also reports their mean and its sum with the Monte Carlo covariance; otherwise
those sampling and combined entries remain missing.
}
