\name{sp_predict_draws}
\alias{sp_predict_draws}
\alias{print.spmodels_prediction_draws}
\title{Joint Predictive Draws from a Fitted Spatial Model}
\usage{
sp_predict_draws(
  object,
  newdata = NULL,
  nsim = 1000L,
  target = c("latent", "observation"),
  type = c("response", "link"),
  seed = NULL,
  parameter_uncertainty = FALSE,
  unseen = c("error", "marginal")
)

\method{print}{spmodels_prediction_draws}(x, ...)
}
\arguments{
  \item{object}{A fitted \code{spmodels_fit} object.}
  \item{newdata}{Optional prediction data containing all predictors and point
  coordinates. The fitted locations are used when this is \code{NULL}.}
  \item{nsim}{Positive number of joint draws.}
  \item{target}{Draw the smooth latent target or a future observation-site
  target. The latter includes the independent variance component.}
  \item{type}{Draw on the response or linear-predictor scale. For binomial and
  Poisson fits, latent response draws are probabilities or intensities;
  observation response draws are Bernoulli outcomes or Poisson counts.}
  \item{seed}{Optional integer seed local to the call.}
  \item{parameter_uncertainty}{Include a local Gaussian approximation to
  fitted regression and covariance-parameter uncertainty? This requires a
  usable fitted Hessian covariance. Gaussian covariance and regression blocks
  are treated as asymptotically independent. Fitted intervention parameters
  rebuild the intervention multiplier and latent state for every draw.}
  \item{unseen}{Policy for grouping levels absent during fitting.
  \code{"error"} rejects them; \code{"marginal"} draws new zero-mean random
  intercepts with the fitted variance, shared by rows carrying the same new
  level.}
  \item{x}{An object returned by \code{sp_predict_draws()}.}
  \item{...}{Reserved for print-method compatibility.}
}
\value{An object of class \code{spmodels_prediction_draws}. Its \code{draws}
component is an \code{nsim} by number-of-locations matrix. The object also
contains the predictive link mean and covariance and target metadata.}
\description{Draws jointly across prediction locations from the exact
conditional Gaussian distribution or the local Gaussian predictive
approximation for a binomial-logit or Poisson-log fit. Joint drawing preserves
the spatial dependence needed for grid and regional targets.}
\details{For non-Gaussian models, response-scale latent draws transform the
joint linear-predictor draw into probabilities or intensities. Observation
draws additionally apply Bernoulli or Poisson response sampling. Regression-
and covariance-parameter uncertainty is omitted by default. When requested,
the returned \code{parameter_draws} record the regression coefficients,
log covariance parameters and any logit-impact/log-decay values used for each
predictive draw. Draw-specific intervention multipliers are also retained.}
\examples{
d <- data.frame(y = c(1, 2, 1.5, 3), x = c(0, 1, 0, 1),
                east = c(0, 1, 0, 1), north = c(0, 0, 1, 1))
fit <- sp_fit(y ~ x, d, coords = ~ east + north,
              control = sp_control(spatial_variance = 1, scale = 1,
                                   nugget_variance = 0.2, estimate = FALSE,
                                   distance_units = "km"))
draws <- sp_predict_draws(fit, d[1:2, ], nsim = 10, seed = 1)
dim(draws$draws)
}
