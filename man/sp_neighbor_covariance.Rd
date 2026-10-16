\name{sp_neighbor_covariance}
\alias{sp_neighbor_covariance}
\alias{print.spmodels_neighbor_covariance}
\title{Inspect Nearest-Predecessor Predictive Covariance}
\usage{
sp_neighbor_covariance(
  object,
  newdata = NULL,
  neighbors = 30L,
  ordering = c("input", "maximin"),
  prediction_order = NULL,
  target = c("latent", "observation"),
  unseen = c("error", "marginal"),
  reference = FALSE
)

\method{print}{spmodels_neighbor_covariance}(x, ...)
}
\arguments{
  \item{object}{A fitted spatial \code{spmodels_fit} object.}
  \item{newdata}{Optional point prediction data.}
  \item{neighbors}{Positive predecessor budget.}
  \item{ordering}{Input or deterministic maximin order.}
  \item{prediction_order}{Optional explicit prediction-row permutation.}
  \item{target}{Latent or future observation-site target.}
  \item{unseen}{Policy for grouping levels absent during fitting.}
  \item{reference}{Also construct the ordinary exact fitted predictive
  covariance and approximation-error metrics?}
  \item{x}{An object returned by \code{sp_neighbor_covariance()}.}
  \item{...}{Reserved for print-method compatibility.}
}
\value{An object containing the implied covariance, predictive mean, local
conditional coefficients/variances, graph and order. With
\code{reference = TRUE}, it also contains the exact fitted predictive
covariance and relative Frobenius, maximum absolute covariance and maximum
marginal-variance errors.}
\description{Reconstructs the complete covariance implied by the
nearest-predecessor sequential conditionals for diagnostic comparison.}
\details{This diagnostic deliberately allocates a complete prediction-site
matrix and is intended for manageable validation subsets. It is not the
bounded-memory production route. The covariance recursion uses each local
conditional coefficient and variance; retaining every predecessor reproduces
the ordinary fitted predictive covariance up to numerical precision.}
\examples{
d <- data.frame(y = c(1, 2, 1.5, 3), x = c(0, 1, 0, 1),
                east = c(0, 1, 0, 1), north = c(0, 0, 1, 1))
fit <- sp_fit(y ~ x, d, coords = ~ east + north,
  control = sp_control(spatial_variance = 1, scale = 1,
    nugget_variance = 0.2, estimate = FALSE, distance_units = "km"))
diagnostic <- sp_neighbor_covariance(fit, d, neighbors = 2,
                                     reference = TRUE)
diagnostic$errors
}
