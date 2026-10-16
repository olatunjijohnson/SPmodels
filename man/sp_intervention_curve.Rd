\name{sp_intervention_curve}
\alias{sp_intervention_curve}
\alias{print.spmodels_intervention_curve}
\alias{plot.spmodels_intervention_curve}
\title{Summarise a Decay-Adjusted Intervention Effect over Time}
\usage{
sp_intervention_curve(object, lag, level = 0.95)

\method{print}{spmodels_intervention_curve}(x, ...)

\method{plot}{spmodels_intervention_curve}(x, ...)
}
\arguments{
\item{object}{An intervention specification or fitted model containing one.}
\item{lag}{Finite non-negative lags in the intervention time units.}
\item{level}{Pointwise confidence level strictly between zero and one.}
\item{x}{An intervention curve.}
\item{...}{Additional graphical arguments for the plot method; unused by print.}
}
\value{An object of class \code{spmodels_intervention_curve} containing an unrounded values table and metadata.}
\description{Calculates the remaining proportional reduction and one-intervention prevalence multiplier. Fitted models with usable intervention-parameter covariance receive delta-method pointwise intervals.}
