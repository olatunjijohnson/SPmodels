\name{sp_target_geometry}
\alias{sp_target_geometry}
\alias{plot.spmodels_target_sf}
\title{Attach Prediction-Target Summaries to Spatial Geometry}
\usage{
sp_target_geometry(object, geometry, id = NULL)

\method{plot}{spmodels_target_sf}(x, value = "mean", ...)
}
\arguments{
\item{object}{An object returned by \code{sp_predict_target()}.}
\item{geometry}{A non-empty \code{sf} object without empty geometries.}
\item{id}{Optional geometry column matched exactly to target identifiers.}
\item{x}{An \code{spmodels_target_sf} object.}
\item{value}{Numeric target-summary column to display.}
\item{...}{Arguments passed to the ordinary \code{sf} plot method.}
}
\value{An \code{sf} object containing unrounded prediction-target summaries
and audit metadata. The plot method returns it invisibly.}
\description{Joins calculated point/grid or grouped polygon summaries to
geometry using row order or explicit unique identifiers.}
