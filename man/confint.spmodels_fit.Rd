\name{confint.spmodels_fit}
\alias{confint.spmodels_fit}
\title{Confidence Intervals for Fitted SPmodels Parameters}
\usage{
\method{confint}{spmodels_fit}(object, parm, level = 0.95, ...)
}
\arguments{
  \item{object}{A fitted \code{spmodels_fit} object.}
  \item{parm}{Optional numeric positions or qualified names such as
  \code{"mean::(Intercept)"} or \code{"covariance::scale"}.}
  \item{level}{Confidence probability strictly between zero and one.}
  \item{...}{Reserved for method compatibility.}
}
\value{A matrix containing lower and upper limits. Its row names identify the
model component and term. Attributes record the level and transformation.}
\description{Constructs local Wald intervals from the fitted standard errors.
Mean coefficients and random intercept-slope covariance entries use the
identity scale, positive parameters use a logarithmic scale, and an estimated
intervention impact uses the logit scale.}
\details{These intervals rely on the fitted local observed-Hessian
approximation. A missing standard error produces missing limits rather than a
fabricated interval. Fixed parameters are omitted. The intervals do not include
Monte Carlo error or repair weak identifiability.}
