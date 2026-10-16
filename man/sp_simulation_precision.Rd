\name{sp_simulation_precision}
\alias{sp_simulation_precision}
\alias{print.spmodels_simulation_precision}
\title{Precision-Directed Simulation Recovery}
\usage{
sp_simulation_precision(
  object,
  initial_nsim = 20L,
  max_nsim = 1000L,
  growth = 2,
  target_mcse,
  min_success = 10L,
  seed = NULL,
  level = 0.95,
  on_error = c("record", "stop"),
  workers = 1L
)

\method{print}{spmodels_simulation_precision}(x, ...)
}
\arguments{
  \item{object}{A fitted \code{spmodels_fit}.}
  \item{initial_nsim}{Initial positive replicate count.}
  \item{max_nsim}{Maximum replicate count, not below the initial count.}
  \item{growth}{Finite multiplicative increase strictly greater than one.}
  \item{target_mcse}{Finite positive target for every recovered parameter's
  bias Monte Carlo standard error.}
  \item{min_success}{Minimum successful refits required per parameter.}
  \item{seed}{Optional base seed. Attempt seeds advance deterministically.}
  \item{level}{Wald interval probability passed to recovery.}
  \item{on_error}{Record failed refits or stop immediately.}
  \item{workers}{Positive number of replicate workers passed to every complete
  recovery attempt.}
  \item{x}{An object returned by \code{sp_simulation_precision()}.}
  \item{...}{Reserved for print-method compatibility.}
}
\value{An object containing the final
\code{spmodels_simulation_recovery}, every attempted replicate count and seed,
parameter-wise achieved MCSE and success count, and the target status.}
\description{Repeats a complete fitted-model recovery experiment at
geometrically increasing replicate counts until the Monte Carlo standard error
of mean bias passes for every recovered parameter or the exact cap is
evaluated.}
\details{For parameter \eqn{j}, the monitored quantity is the empirical
standard deviation of successful replicate estimates divided by the square
root of their count. Each attempt is a complete independently seeded
experiment. Failure at the cap remains explicitly false; criteria are never
relaxed to manufacture a pass.}
