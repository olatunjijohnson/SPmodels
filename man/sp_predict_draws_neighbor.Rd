\name{sp_predict_draws_neighbor}
\alias{sp_predict_draws_neighbor}
\title{Bounded-Memory Approximate Joint Predictive Draws}
\usage{
sp_predict_draws_neighbor(
  object,
  newdata = NULL,
  nsim = 1000L,
  neighbors = 30L,
  ordering = c("input", "maximin"),
  prediction_order = NULL,
  target = c("latent", "observation"),
  type = c("response", "link"),
  seed = NULL,
  parameter_uncertainty = FALSE,
  unseen = c("error", "marginal")
)
}
\arguments{
  \item{object}{A fitted spatial \code{spmodels_fit} object.}
  \item{newdata}{Optional prediction data containing all predictors and point
  coordinates. The fitted locations are used when this is \code{NULL}.}
  \item{nsim}{Positive number of joint draws.}
  \item{neighbors}{Positive maximum number of earlier prediction locations
  retained in each conditional distribution.}
  \item{ordering}{Use input-row order or a deterministic maximin order that
  starts at row one and repeatedly selects the row farthest from the selected
  set. Ignored when \code{prediction_order} is supplied.}
  \item{prediction_order}{Optional permutation of prediction-row positions.
  Input order is used by default; returned columns retain input order.}
  \item{target}{Draw the smooth latent target or a future observation-site
  target.}
  \item{type}{Draw on the response or linear-predictor scale.}
  \item{seed}{Optional integer seed local to the call.}
  \item{parameter_uncertainty}{Include the fitted local Gaussian parameter
  approximation? Every parameter draw rebuilds the fitted state and complete
  nearest-predecessor factorization.}
  \item{unseen}{Policy for grouping levels absent during fitting.}
}
\value{An object of class \code{spmodels_prediction_draws}. Its full
\code{predictive_covariance} is \code{NULL}; the conditional variances,
neighbor indices and processing order define the approximation.}
\description{Generates joint draws from a nearest-predecessor factorization of
the fitted Gaussian predictive distribution without allocating the complete
prediction-site covariance matrix.}
\details{For each prediction row, the function conditions on at most
\code{neighbors} earlier rows selected by Euclidean coordinate distance. The
model's exact Gaussian or local non-Gaussian predictive mean and covariance
formula is used for that small subset. Retaining every predecessor reproduces
the ordinary joint Gaussian factorization up to numerical precision. Smaller
neighborhoods are an order-dependent approximation. Maximin ordering uses
linear working storage but quadratic distance work. Optional parameter
uncertainty uses the same fitted-Hessian approximation as
\code{sp_predict_draws()} and can be computationally expensive because every
parameter draw rebuilds the fitted state and all local conditionals. The
training fit remains dense. Singular positive-semidefinite predecessor blocks,
such as coincident latent prediction rows, use a checked eigen-based
Moore--Penrose solve without artificial jitter.}
\examples{
d <- data.frame(y = c(1, 2, 1.5, 3), x = c(0, 1, 0, 1),
                east = c(0, 1, 0, 1), north = c(0, 0, 1, 1))
fit <- sp_fit(y ~ x, d, coords = ~ east + north,
  control = sp_control(spatial_variance = 1, scale = 1,
    nugget_variance = 0.2, estimate = FALSE, distance_units = "km"))
draws <- sp_predict_draws_neighbor(fit, d, nsim = 10, neighbors = 2,
                                   seed = 1)
dim(draws$draws)
}
