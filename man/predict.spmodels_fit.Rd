\name{predict.spmodels_fit}
\alias{predict.spmodels_fit}
\title{Predict from a Fitted SPmodels Model}
\usage{
\method{predict}{spmodels_fit}(
  object,
  newdata = NULL,
  se.fit = FALSE,
  target = c("latent", "observation"),
  type = c("response", "link"),
  unseen = c("error", "marginal"),
  ...
)
}
\arguments{
  \item{object}{A fitted \code{spmodels_fit} object.}
  \item{newdata}{Optional data frame containing predictors and coordinates.}
  \item{se.fit}{Whether to return prediction standard errors.}
  \item{target}{Predict the latent spatial process or a future noisy
  observation-site target. Observation uncertainty includes nugget variance
  and, on the non-Gaussian response scale, response sampling variation.}
  \item{type}{Return predictions on the response or linear-predictor scale.
  These coincide for the Gaussian identity link.}
  \item{unseen}{Policy for grouping levels absent during fitting.
  \code{"error"} rejects them. \code{"marginal"} integrates a new zero-mean
  random intercept with its fitted variance and preserves dependence between
  prediction rows sharing the same new level.}
  \item{...}{Reserved for method compatibility.}
}
\value{A numeric vector, or a list containing \code{fit}, \code{se.fit},
\code{target}, and, for non-Gaussian fits, \code{type}.}
\description{Computes conditional Gaussian predictions using universal
kriging equations. For binomial and Poisson Laplace fits, computes an
approximate Gaussian predictive linear predictor and integrates it through
the response link. Prediction uncertainty is conditional on fitted regression
and covariance parameters.}
