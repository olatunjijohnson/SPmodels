\name{sp_simulation_recovery}
\alias{sp_simulation_recovery}
\alias{print.spmodels_simulation_recovery}
\title{Repeated Fitted-Model Simulation and Parameter Recovery}
\usage{
sp_simulation_recovery(
  object,
  nsim = 100L,
  seed = NULL,
  level = 0.95,
  on_error = c("record", "stop"),
  workers = 1L
)

\method{print}{spmodels_simulation_recovery}(x, ...)
}
\arguments{
  \item{object}{A fitted \code{spmodels_fit}.}
  \item{nsim}{Positive number of response simulation/refit replicates.}
  \item{seed}{Optional integer seed local to the operation.}
  \item{level}{Wald interval probability strictly between zero and one.}
  \item{on_error}{Record failed refits and continue, or stop on failure.
  Sequential execution stops immediately; already-dispatched parallel work
  may finish before the error is returned.}
  \item{workers}{Positive number of replicate workers. Values above one use
  forked processes on non-Windows systems. MCML recovery then requires an
  explicit seed.}
  \item{x}{An object returned by \code{sp_simulation_recovery()}.}
  \item{...}{Reserved for print-method compatibility.}
}
\value{An object containing simulated responses, replicate estimates and
standard errors, convergence codes, recorded failures, fitted truth and a
parameter-level recovery summary.}
\description{Simulates response datasets from the supplied fitted model and
refits the same declared estimation route to every replicate. The summary
reports bias, RMSE, empirical standard deviation, mean reported standard error
and Wald interval coverage. MCML is repeated when the supplied fit used MCML.}
\details{The fitted parameter values are treated as simulation truth. This is
a conditional fitted-model assessment, not evidence of performance across a
scientific parameter space. Named responses and two-column binomial responses
are supported. External offset arguments must instead be written in the model
formula so every refit is self-contained. The recovery vector includes every
estimated mean, spatial/nugget covariance, grouped random-intercept variance,
temporal-scale and intervention parameter; fixed quantities are not reported
as recovered estimates.}
