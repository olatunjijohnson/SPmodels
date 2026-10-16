\name{spacetime_cor}
\alias{spacetime_cor}
\title{Separable Matern-Exponential Space-Time Correlation}
\usage{spacetime_cor(distance, time_lag, spatial_scale, smoothness, temporal_scale)}
\arguments{
\item{distance}{Non-negative finite spatial distances.}
\item{time_lag}{Non-negative finite time lags with dimensions matching \code{distance}.}
\item{spatial_scale}{Positive Matern spatial scale.}
\item{smoothness}{Positive Matern smoothness.}
\item{temporal_scale}{Positive exponential temporal scale in the same units as \code{time_lag}.}
}
\value{Correlations with dimensions and names inherited from \code{distance}.}
\description{Multiplies stationary Matern spatial correlation by stationary exponential temporal correlation.}
