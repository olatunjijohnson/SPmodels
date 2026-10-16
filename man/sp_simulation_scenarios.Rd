\name{sp_simulation_scenarios}
\alias{sp_simulation_scenarios}
\alias{print.spmodels_simulation_scenarios}
\title{Run Parameter-Recovery Experiments Across Fitted Scenarios}
\usage{
sp_simulation_scenarios(
  objects,
  nsim = 100L,
  seed = NULL,
  level = 0.95,
  on_error = c("record", "stop"),
  workers = 1L
)

\method{print}{spmodels_simulation_scenarios}(x, ...)
}
\arguments{
  \item{objects}{A named, non-empty list of fitted \code{spmodels_fit}
  objects. Names are retained as scenario identifiers.}
  \item{nsim}{One positive integer replicate count for every scenario, or one
  count per scenario. Named values are matched exactly to scenario names.}
  \item{seed}{Optional base integer seed. Successive scenarios use successive
  deterministic seeds.}
  \item{level}{Wald interval probability strictly between zero and one.}
  \item{on_error}{Record failed refits and continue, or stop immediately.}
  \item{workers}{Positive number of scenario workers. Values above one use
  forked processes on non-Windows systems; Windows currently requires one.}
  \item{x}{An object returned by \code{sp_simulation_scenarios()}.}
  \item{...}{Reserved for print-method compatibility.}
}
\value{An object containing every named
\code{spmodels_simulation_recovery} result, their row-bound parameter summary,
and scenario-level family, replicate, seed and successful-refit metadata.}
\description{Runs the existing fitted-model simulation and refitting contract
over multiple explicitly named fitted scenarios. Each result remains available
for audit; the combined summary adds a scenario column without pooling unlike
parameters or replacing individual failures.}
\details{The scenario models define their own fitted simulation truths and
designs. This function orchestrates separate conditional experiments; it does
not claim that the supplied scenarios form a scientifically sufficient design.
When a base seed is supplied, each scenario is reproducible and the caller's
random-number state remains unchanged. Scenario seeds are fixed before
dispatch, so sequential and parallel execution are identical.}
