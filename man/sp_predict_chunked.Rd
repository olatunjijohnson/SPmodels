\name{sp_predict_chunked}
\alias{sp_predict_chunked}
\title{Predict a Large Point Set in Bounded-Memory Chunks}
\usage{
sp_predict_chunked(
  object,
  newdata,
  chunk_size = 1000L,
  se.fit = FALSE,
  target = c("latent", "observation"),
  type = c("response", "link"),
  unseen = c("error", "marginal")
)
}
\arguments{
  \item{object}{A fitted \code{spmodels_fit} object.}
  \item{newdata}{Non-empty prediction data containing every required
  predictor, coordinate, time and grouping column.}
  \item{chunk_size}{Positive maximum number of prediction rows per chunk.}
  \item{se.fit}{Return pointwise standard errors as well as fitted values?}
  \item{target}{Predict the latent process or a future observation.}
  \item{type}{Return response- or link-scale predictions.}
  \item{unseen}{Policy for grouping levels absent during fitting.}
}
\value{The same numeric vector or pointwise list schema as \code{predict()},
with attributes \code{chunk_size} and \code{chunks}.}
\description{Applies ordinary point prediction to consecutive row chunks and
combines results in original order. This is an exact computational
reorganization, not a covariance approximation.}
\details{The complete training-by-prediction cross-covariance matrix is never
held in memory. Peak cross-covariance storage is proportional to the number of
training rows times \code{chunk_size}. The function returns pointwise means
and optional marginal standard errors; use \code{sp_predict_draws()} when a
joint covariance or dependent draws across all prediction rows are required.}
\examples{
d <- data.frame(y = c(1, 2, 1.5, 3), x = c(0, 1, 0, 1),
                east = c(0, 1, 0, 1), north = c(0, 0, 1, 1))
fit <- sp_fit(y ~ x, d, coords = ~ east + north,
              control = sp_control(spatial_variance = 1, scale = 1,
                nugget_variance = .2, estimate = FALSE,
                distance_units = "km"))
sp_predict_chunked(fit, d, chunk_size = 2L)
}
