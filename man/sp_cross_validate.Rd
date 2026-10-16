\name{sp_cross_validate}
\alias{sp_cross_validate}
\alias{print.spmodels_cross_validation}
\title{Fixed-Fold Cross-Validation for an SPmodels Fit}
\usage{
sp_cross_validate(
  object,
  folds = 5L,
  seed = NULL,
  unseen = c("error", "marginal")
)

\method{print}{spmodels_cross_validation}(x, ...)
}
\arguments{
  \item{object}{A fitted \code{spmodels_fit} object retaining aligned source
  data and controls.}
  \item{folds}{An integer number of balanced randomly assigned folds, or one
  fixed non-missing fold identifier per fitted observation, or an object from
  \code{sp_buffered_folds()}, \code{sp_block_folds()} or
  \code{sp_environmental_folds()}.}
  \item{seed}{Optional local seed used only when \code{folds} is a count.}
  \item{unseen}{Policy for held-out grouping levels absent from a training
  fold. The safe default errors; \code{"marginal"} integrates new zero-mean
  grouped effects.}
  \item{x}{An object returned by \code{sp_cross_validate()}.}
  \item{...}{Reserved for print-method compatibility.}
}
\value{An object of class \code{spmodels_cross_validation} containing held-out
row predictions, overall metrics, fold metrics, resolved fold identifiers,
per-fold training sizes and optional buffer distance.}
\description{Refits the original model after withholding every declared fold,
then predicts held-out observations. Reports RMSE, mean absolute error and mean
negative log score. Fixed user-supplied fold identifiers are preserved.}
\details{Cross-validation currently requires any offset to be written inside
the model formula. Random folds are a convenience for reproducible checks;
spatially separated folds should normally be supplied for spatial validation.
Buffered folds remove near-test rows from training independently for each
fold. Group-constrained folds that withhold complete fitted random-effect
levels generally require the explicit marginal unseen-level policy.}
