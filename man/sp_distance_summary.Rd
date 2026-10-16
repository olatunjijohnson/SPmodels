\name{sp_distance_summary}
\alias{sp_distance_summary}
\alias{print.spmodels_distance_summary}
\title{Summarise Pairwise Spatial Distances}
\usage{
sp_distance_summary(
  data,
  coords = NULL,
  distance_units = NULL,
  duplicates = c("include", "exclude", "error"),
  probs = c(0, 0.25, 0.5, 0.75, 1)
)

\method{print}{spmodels_distance_summary}(x, ...)
}
\arguments{
\item{data}{A spatial SPmodels fit, projected point \code{sf}, or data frame.}
\item{coords}{A two-column coordinate formula for an ordinary data frame.}
\item{distance_units}{Explicit coordinate units for an ordinary data frame.}
\item{duplicates}{Include, exclude, or reject repeated-location pairs.}
\item{probs}{Finite quantile probabilities between zero and one.}
\item{x}{A distance-summary object.}
\item{...}{Reserved for method compatibility.}
}
\value{An object containing retained pair distances, a one-row summary, units
and the duplicate policy.}
\description{Calculates all unordered Euclidean coordinate-pair distances and
reports counts, duplicate pairs, mean, standard deviation and requested
type-8 quantiles.}
