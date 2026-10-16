\name{sp_fit_neighbor}
\alias{sp_fit_neighbor}
\alias{coef.spmodels_neighbor_fit}
\alias{confint.spmodels_neighbor_fit}
\alias{vcov.spmodels_neighbor_fit}
\alias{logLik.spmodels_neighbor_fit}
\alias{fitted.spmodels_neighbor_fit}
\alias{residuals.spmodels_neighbor_fit}
\alias{print.spmodels_neighbor_fit}
\alias{summary.spmodels_neighbor_fit}
\alias{print.summary.spmodels_neighbor_fit}
\alias{predict.spmodels_neighbor_fit}
\title{Fit an approximate Gaussian nearest-predecessor model}
\usage{
sp_fit_neighbor(
  formula,
  data,
  coords,
  neighbors = 30L,
  ordering = c("input", "maximin"),
  observation_order = NULL,
  subset = NULL,
  na.action = stats::na.fail,
  offset = NULL,
  control = sp_control(),
  verbose = FALSE
)

\method{coef}{spmodels_neighbor_fit}(object, ...)
\method{confint}{spmodels_neighbor_fit}(object, parm, level = 0.95, ...)
\method{vcov}{spmodels_neighbor_fit}(object, ...)
\method{logLik}{spmodels_neighbor_fit}(object, ...)
\method{fitted}{spmodels_neighbor_fit}(object, ...)
\method{residuals}{spmodels_neighbor_fit}(object, ...)
\method{print}{spmodels_neighbor_fit}(x, ...)
\method{summary}{spmodels_neighbor_fit}(object, ...)
\method{print}{summary.spmodels_neighbor_fit}(x, ...)
\method{predict}{spmodels_neighbor_fit}(
  object,
  newdata = NULL,
  se.fit = FALSE,
  target = c("latent", "observation"),
  ...
)
}
\arguments{
  \item{formula, data}{A model formula and data frame.}
  \item{coords}{A one-sided formula selecting exactly two numeric coordinate
  columns. Coordinates must already be projected.}
  \item{neighbors}{Positive maximum number of earlier observations retained
  in each conditional density.}
  \item{ordering}{Use input-row order or deterministic maximin order.}
  \item{observation_order}{Optional permutation of observation-row positions.
  Supplying it overrides \code{ordering}.}
  \item{subset, na.action, offset}{Standard model-frame arguments.}
  \item{control}{Settings from \code{sp_control()}. The initial approximation
  supports the Gaussian spatial covariance only, without grouped or temporal
  covariance components.}
  \item{verbose}{Logical scalar controlling optimization tracing.}
  \item{object, x}{A fitted nearest-predecessor model or its summary.}
  \item{parm}{Optional fitted-parameter names or positions for confidence
  intervals.}
  \item{level}{Confidence probability strictly between zero and one.}
  \item{newdata}{Optional prediction data containing predictors and point
  coordinates. With \code{NULL}, return sequential training fitted values.}
  \item{se.fit}{Return approximate standard errors?}
  \item{target}{Predict the smooth latent target or a future observation.}
  \item{...}{Reserved for generic compatibility.}
}
\value{
\code{sp_fit_neighbor()} returns a \code{spmodels_neighbor_fit} object. It
retains the processing order, predecessor graph, local coefficients and
conditional variances. Generic methods return the corresponding coefficients,
covariance, likelihood, fitted values, residuals, summaries or local
conditional predictions.
}
\description{
Fits a Gaussian identity-link geostatistical regression by multiplying ordered
conditional densities that retain at most \code{neighbors} earlier
observations. Only covariance blocks of at most that size are factorized.
Retaining every predecessor recovers the dense Gaussian likelihood, regression
coefficient estimates and coefficient covariance up to numerical precision.
The approximation is explicit and never silently replaces \code{sp_fit()}.
}
\details{
The ordering and neighbor budget are part of the statistical approximation and
are retained in the fitted object. Maximin ordering starts at row one and then
selects the point farthest from the selected set. Prediction conditions each
new point on its nearest training observations and does not claim a joint
prediction covariance. Fitted-parameter standard errors use local observed
curvature when it is positive definite.

This initial route deliberately rejects non-Gaussian, grouped, random-slope and
space-time covariance models. Use the ordinary fitters for those models.
}
\examples{
d <- data.frame(
  y = c(1.2, 0.8, 1.7, 1.1, 1.5, 0.9),
  z = c(0, 1, 0, 1, 0, 1),
  east = c(0, 1, 2, 0, 1, 2),
  north = c(0, 0, 0, 1, 1, 1)
)
fit <- sp_fit_neighbor(
  y ~ z, d, coords = ~ east + north, neighbors = 2,
  control = sp_control(
    estimate = FALSE, smoothness = 0.5, spatial_variance = 1,
    scale = 1, nugget_variance = 0.1, distance_units = "km"
  )
)
summary(fit)
predict(fit, d[1:2, ], se.fit = TRUE)
}
