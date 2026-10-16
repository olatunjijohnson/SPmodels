\name{sp_mcml_precision}
\alias{sp_mcml_precision}
\title{Escalate MCML sample size to declared precision targets}
\usage{
sp_mcml_precision(
  object,
  control = sp_mcml_control(),
  target_relative_mcse = 0.02,
  min_ess_fraction = 0.1,
  max_nsim = 32000L,
  growth = 2,
  seed = NULL,
  refresh = TRUE,
  verbose = FALSE
)
}
\arguments{
\item{object}{A fitted binomial-logit or Poisson-log \code{spmodels_fit}.}
\item{control}{Initial settings returned by \code{sp_mcml_control()}.}
\item{target_relative_mcse}{Positive target for relative likelihood Monte Carlo standard error.}
\item{min_ess_fraction}{Required effective-sample fraction strictly between zero and one.}
\item{max_nsim}{Largest permitted importance-sample count.}
\item{growth}{Multiplicative sample-size growth greater than one.}
\item{seed}{Optional integer defining consecutive, call-local attempt seeds.}
\item{refresh}{Whether an unsuccessful fit supplies the next proposal centre and curvature.}
\item{verbose}{Logical scalar passed to \code{sp_mcml()}.}
}
\value{
The final \code{spmodels_fit}, with an \code{mcml_precision} component that
records every attempted sample size, diagnostic and stopping reason.
Adaptive attempts include pilot ESS, maximum weight, relative MCSE and the
used/fallback decision.
}
\description{
Refits using geometrically increasing independent importance-sample sizes
until both declared diagnostic targets are met or \code{max_nsim} is reached.
}
