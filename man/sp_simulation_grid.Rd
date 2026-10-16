\name{sp_simulation_grid}
\alias{sp_simulation_grid}
\title{Construct a Cartesian Grid of Simulation-Recovery Scenarios}
\usage{
sp_simulation_grid(object, values, prefix = "scenario")
}
\arguments{
  \item{object}{A fitted \code{spmodels_fit} object defining the design,
  response family and refitting specification.}
  \item{values}{A named list of non-empty finite numeric vectors. Names must
  identify parameters that the fitted specification estimates. Random-effect
  variances use names such as \code{random::clinic}; estimated grouped slope
  covariance entries use their qualified \code{random_slope::} model-table
  names.}
  \item{prefix}{A non-empty prefix for deterministic scenario names.}
}
\value{A named list of generative \code{spmodels_fit} scenario objects. Its
\code{grid} attribute records the Cartesian product in scenario order.}
\description{Builds a reproducible Cartesian parameter design for
\code{sp_simulation_scenarios()}. Each object retains the fitted design and
refitting specification while its generative mean and covariance state are
rebuilt at the declared truth.}
\details{These are simulation scenario objects, not new analyses of the
observed response. Parameters fixed by the model specification cannot be
varied because they would not be recoverable on refitting. Variance, scale and
decay values must be positive; an estimated intervention alpha must be strictly
between zero and one. Every grouped slope covariance combination must remain
positive definite.}
