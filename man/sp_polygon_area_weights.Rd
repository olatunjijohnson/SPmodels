\name{sp_polygon_area_weights}
\alias{sp_polygon_area_weights}
\title{Construct Exact Fractional Cell-Area Polygon Weights}
\usage{
sp_polygon_area_weights(
  cells,
  polygons,
  region,
  density = NULL,
  unmatched = c("error", "drop"),
  empty_regions = c("error", "drop")
)
}
\arguments{
\item{cells}{Projected polygon \code{sf} prediction cells in draw-column order.}
\item{polygons}{Projected non-overlapping target polygons in the same CRS.}
\item{region}{Column name holding non-missing region identifiers.}
\item{density}{Optional finite non-negative cell values or a column name;
weights are intersection area times density.}
\item{unmatched}{Error for or drop cells with no positive-area intersection.}
\item{empty_regions}{Error for or record/drop regions with no positive-area
cell intersection.}
}
\value{An \code{spmodels_polygon_weights} object accepted by
\code{sp_predict_target()}, retaining cell indices, region groups, raw weights,
intersection areas, cell fractions and excluded cells/regions.}
\description{Uses exact projected polygon intersections so a prediction cell
can contribute proportionally to more than one target region.}
