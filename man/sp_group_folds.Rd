\name{sp_group_folds}
\alias{sp_group_folds}
\alias{print.spmodels_group_folds}
\title{Construct Group-Constrained Validation Folds}
\usage{
sp_group_folds(
  object,
  group,
  folds = 5L,
  assignment = c("balanced", "random"),
  seed = NULL
)

\method{print}{spmodels_group_folds}(x, ...)
}
\arguments{
\item{object}{A fitted \code{spmodels_fit} object.}
\item{group}{Name of one grouping column in the fitted data.}
\item{folds}{Number of folds, not exceeding the observed group count.}
\item{assignment}{Deterministic size-first or locally seeded randomized group
order, followed by least-loaded allocation.}
\item{seed}{Optional local seed for randomized assignment.}
\item{x}{A group-fold object.}
\item{...}{Reserved for print-method compatibility.}
}
\value{An integer fold vector in fitted-row order. Group names, row counts and
allocations are retained in the \code{groups} attribute.}
\description{Allocates whole observed groups to validation folds, preventing
the same group from appearing in both training and held-out rows. Allocation
greedily assigns each next group to the fold with the fewest rows.}
