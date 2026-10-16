\name{sp_prediction_update}
\alias{sp_prediction_update}
\alias{print.spmodels_prediction_update}
\title{Compare predictions under an updated predictor scenario}
\usage{
sp_prediction_update(
  object,
  newdata,
  changes,
  rows = NULL,
  se.fit = FALSE,
  target = c("latent", "observation"),
  type = c("response", "link")
)

\method{print}{spmodels_prediction_update}(x, ...)
}
\arguments{
\item{object}{A fitted \code{spmodels_fit} object.}
\item{newdata}{Non-empty prediction data accepted by \code{predict()}.}
\item{changes}{A non-empty named list replacing fitted predictor or
in-formula offset variables. Each value has length one or the number of
selected rows. Named values must exactly match the selected row names.}
\item{rows}{Rows to update: \code{NULL} for all rows, a logical vector, unique
integer indices, or unique row names.}
\item{se.fit}{Include baseline and updated pointwise standard errors.}
\item{target}{Predict the latent process or a future observation.}
\item{type}{Return response- or link-scale predictions.}
\item{x}{An object returned by \code{sp_prediction_update()}.}
\item{...}{Reserved for print-method compatibility.}
}
\value{
An object of class \code{spmodels_prediction_update}. It contains the updated
data; baseline and updated design matrices and offsets; row-level baseline,
updated and difference values; optional standard errors; changed rows and
columns; and preservation checks.
}
\description{
Constructs an explicit prediction scenario without refitting or mutating the
fitted model. Coordinates, geometry, row order and unmentioned columns remain
unchanged. Formula transformations, factors and \code{offset()} expressions
are rebuilt through the fitted model's terms and contrast definitions.
}
\details{
Coordinate variables cannot be changed. Offsets supplied separately to the
original fitting call are rejected because they cannot be reconstructed safely
from prediction data; write an offset in the model formula for scenario use.
The returned differences compare point predictions conditional on the same
fitted parameter values and fitted latent state.
}
