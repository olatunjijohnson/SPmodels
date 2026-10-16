\name{sp_buffered_folds}
\alias{sp_buffered_folds}
\alias{print.spmodels_buffered_folds}
\title{Construct Buffered Spatial Validation Folds}
\usage{
sp_buffered_folds(object, folds = 5L, buffer, seed = NULL, nstart = 20L)

\method{print}{spmodels_buffered_folds}(x, ...)
}
\arguments{
  \item{object}{A fitted spatial \code{spmodels_fit} model.}
  \item{folds}{Number of geometry-driven test folds.}
  \item{buffer}{Finite positive exclusion distance in fitted coordinate units.}
  \item{seed}{Optional local seed for fold construction.}
  \item{nstart}{Positive number of k-means starts.}
  \item{x}{An object returned by \code{sp_buffered_folds()}.}
  \item{...}{Reserved for print-method compatibility.}
}
\value{An integer fold vector carrying the exclusion distance, coordinate
units, cluster centres, seed and k-means metadata.}
\description{Constructs geometry-driven held-out folds and declares a spatial
buffer. During cross-validation, rows outside the held-out fold but no farther
than the buffer from any test location are excluded from that fold's training
fit.}
\details{The buffer changes training membership independently for every fold;
buffered rows are neither trained on nor scored in that fold. Test rows remain
exactly those identified by the fold vector. An over-large buffer that leaves
no usable training data fails explicitly.}
