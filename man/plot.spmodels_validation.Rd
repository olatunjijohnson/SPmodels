\name{plot.spmodels_cross_validation}
\alias{plot.spmodels_cross_validation}
\alias{plot.spmodels_simulation_assessment}
\alias{plot.spmodels_simulation_recovery}
\title{Display validation and simulation results}
\usage{
\method{plot}{spmodels_cross_validation}(
  x,
  type = c("observed", "error", "fold"),
  xlab = NULL,
  ylab = NULL,
  main = NULL,
  ...
)

\method{plot}{spmodels_simulation_assessment}(
  x,
  type = c("prediction", "error"),
  xlab = NULL,
  ylab = NULL,
  main = NULL,
  ...
)

\method{plot}{spmodels_simulation_recovery}(
  x,
  parameter = 1L,
  intervals = TRUE,
  xlab = "Simulation replicate",
  ylab = NULL,
  main = NULL,
  ...
)
}
\arguments{
\item{x}{A cross-validation, predictive-assessment or simulation-recovery object.}
\item{type}{The requested view for cross-validation or predictive assessment.}
\item{parameter}{Name or index of one recorded recovery parameter.}
\item{intervals}{Whether to draw reported Wald intervals when available.}
\item{xlab, ylab, main}{Plot labels.}
\item{...}{Further arguments passed to the principal base-graphics call.}
}
\value{The input object, invisibly.}
\description{
Base-graphics views over the unrounded values already stored in SPmodels
validation and simulation results. No estimates or summaries are recomputed
for display, apart from the displayed Wald endpoints in a recovery plot.
}
