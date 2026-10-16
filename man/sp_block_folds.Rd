\name{sp_block_folds}
\alias{sp_block_folds}
\alias{print.spmodels_block_folds}
\title{Construct Regular Spatial-Block Folds}
\usage{
sp_block_folds(
  object,
  block_size,
  folds = 5L,
  origin = NULL,
  assignment = c("systematic", "random"),
  seed = NULL
)

\method{print}{spmodels_block_folds}(x, ...)
}
\arguments{
  \item{object}{A fitted spatial \code{spmodels_fit}.}
  \item{block_size}{One finite positive block width, or separate widths for
  the two coordinate axes, in fitted distance units.}
  \item{folds}{Number of validation folds, not exceeding the number of
  occupied blocks.}
  \item{origin}{Optional two-coordinate block-grid origin. Coordinate-wise
  minima are used by default.}
  \item{assignment}{Allocate canonically ordered blocks round-robin, or
  randomly permute whole blocks before balanced allocation.}
  \item{seed}{Optional integer seed local to random allocation.}
  \item{x}{An object returned by \code{sp_block_folds()}.}
  \item{...}{Reserved for print-method compatibility.}
}
\value{An integer fold vector in fitted-row order with class
\code{spmodels_block_folds}. Occupied-block indices, allocation, block size,
origin, units and seed are retained as attributes.}
\description{Partitions projected fitted coordinates into axis-aligned
rectangular blocks and assigns whole occupied blocks to validation folds.
Repeated coordinates and all rows in one block remain together.}
\details{Systematic assignment follows blocks from top to bottom and left to
right, cycling through fold labels. Random assignment permutes whole blocks
under a caller-local seed before balanced allocation. This controls the
support withheld together; it does not by itself impose a distance buffer.}
