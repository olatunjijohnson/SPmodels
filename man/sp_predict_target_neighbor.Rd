\name{sp_predict_target_neighbor}
\alias{sp_predict_target_neighbor}
\title{Stream Bounded-Memory Approximate Predictive Targets}
\usage{
sp_predict_target_neighbor(
  object,
  newdata = NULL,
  weights,
  groups = NULL,
  nsim = 1000L,
  neighbors = 30L,
  ordering = c("input", "maximin"),
  prediction_order = NULL,
  target = c("latent", "observation"),
  type = c("response", "link"),
  threshold = NULL,
  probs = c(0.025, 0.5, 0.975),
  seed = NULL,
  unseen = c("error", "marginal")
)
}
\arguments{
  \item{object}{A fitted spatial \code{spmodels_fit} object.}
  \item{newdata}{Optional point prediction data.}
  \item{weights}{Finite non-negative weights matching prediction rows, or a
  polygon-weight object.}
  \item{groups}{Optional target identifier per prediction row. It requires
  numeric \code{weights}; the default is one target.}
  \item{nsim}{Positive number of draws.}
  \item{neighbors}{Positive predecessor budget.}
  \item{ordering}{Input or deterministic maximin processing order.}
  \item{prediction_order}{Optional explicit permutation of prediction rows.}
  \item{target}{Latent or future observation-site target.}
  \item{type}{Response or link scale.}
  \item{threshold}{Optional finite target exceedance threshold.}
  \item{probs}{Finite target quantile probabilities in zero--one.}
  \item{seed}{Optional integer seed local to the call.}
  \item{unseen}{Policy for grouping levels absent during fitting.}
}
\value{An \code{spmodels_prediction_target} containing target-level draws and
summaries, sparse normalized weights, the neighbor graph, conditional
variances, processing order and peak number of point-draw vectors retained.}
\description{Generates nearest-predecessor joint prediction draws and
immediately accumulates grouped weighted targets, avoiding storage of a
complete location-level draw matrix when only aggregate targets are needed.}
\details{A point-draw vector is discarded after its final use in the declared
neighbor graph. Working draw storage is therefore proportional to the graph's
active width rather than automatically to every prediction location, in
addition to the returned number-of-draws by number-of-targets matrix. Numeric
and polygon allocation weights are normalized within target. Parameter
uncertainty is not included. Observation-scale binomial/Poisson sampling is
generated during streaming, so its seeded draws need not equal a post-hoc call
that consumes random numbers in a different order.}
\examples{
d <- data.frame(y = c(1, 2, 1.5, 3), x = c(0, 1, 0, 1),
                east = c(0, 1, 0, 1), north = c(0, 0, 1, 1))
fit <- sp_fit(y ~ x, d, coords = ~ east + north,
  control = sp_control(spatial_variance = 1, scale = 1,
    nugget_variance = 0.2, estimate = FALSE, distance_units = "km"))
target <- sp_predict_target_neighbor(
  fit, d, weights = rep(1, nrow(d)), nsim = 10,
  neighbors = 2, seed = 1
)
target$summary
}
