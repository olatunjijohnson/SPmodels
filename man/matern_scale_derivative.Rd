\name{matern_scale_derivative}
\alias{matern_scale_derivative}
\title{Matérn Correlation Derivatives with Respect to Scale}
\usage{
matern_scale_derivative(distance, scale, smoothness, order = 1L)
}
\arguments{
  \item{distance}{Numeric finite non-negative distances.}
  \item{scale}{One finite positive scale in the distance unit.}
  \item{smoothness}{One finite positive Matérn smoothness.}
  \item{order}{Derivative order, either 1 or 2.}
}
\value{A numeric object preserving the dimensions and names of
\code{distance}.}
\description{Calculates the first or second derivative of the SPmodels Matérn
correlation with respect to scale, holding distance and smoothness fixed.}
