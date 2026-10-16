\name{sp_temporal}
\alias{sp_temporal}
\alias{print.spmodels_temporal}
\title{Specify Survey Time for a Separable Space-Time Model}
\usage{
sp_temporal(time, scale, time_units, estimate = FALSE)

\method{print}{spmodels_temporal}(x, ...)
}
\arguments{
\item{time}{One-sided formula selecting one finite numeric time column.}
\item{scale}{Positive exponential temporal scale or its starting value.}
\item{time_units}{Non-empty label for the time and scale units.}
\item{estimate}{Estimate temporal scale jointly with other parameters?}
\item{x}{A temporal specification.}
\item{...}{Unused.}
}
\value{An object of class \code{spmodels_temporal}.}
\description{Declares survey times and the exponential temporal scale for a separable Matern-exponential space-time covariance.}
