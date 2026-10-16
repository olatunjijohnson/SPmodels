\name{sp_polygon_weights}
\alias{sp_polygon_weights}
\title{Assign Point Prediction Weights to Polygon Regions}
\usage{
sp_polygon_weights(
  points,
  polygons,
  region,
  weight = NULL,
  unmatched = c("error", "drop")
)
}
\arguments{
  \item{points}{A non-empty projected point \code{sf} object in the same row
  order as the prediction draws.}
  \item{polygons}{A non-empty projected polygon \code{sf} object using the
  same CRS.}
  \item{region}{One column name in \code{polygons} containing non-missing
  region identifiers.}
  \item{weight}{Optional numeric point weights or one numeric column name in
  \code{points}. Equal raw weights are used by default.}
  \item{unmatched}{Fail for points outside every polygon, or drop them from
  subsequent aggregation.}
}
\value{An object of class \code{spmodels_polygon_weights} containing retained
point indices, raw weights, region identifiers, CRS and unmatched indices. It
can be passed directly as \code{weights} to \code{sp_predict_target()}.}
\description{Uses projected point-in-polygon membership to define regional
weights for joint predictive draws. Points assigned to multiple distinct
region identifiers fail explicitly. This is centre-point allocation, not an
exact cell-area overlay.}
