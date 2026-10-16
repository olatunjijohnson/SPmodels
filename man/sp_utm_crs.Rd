\name{sp_utm_crs}
\alias{sp_utm_crs}
\alias{print.spmodels_crs_choice}
\title{Choose a Projected UTM Coordinate Reference System}
\usage{
sp_utm_crs(
  data,
  zone_policy = c("error", "centre"),
  hemisphere_policy = c("error", "centre")
)

\method{print}{spmodels_crs_choice}(x, ...)
}
\arguments{
\item{data}{A non-empty \code{sf} or \code{sfc} object with known CRS.}
\item{zone_policy}{Error on multiple zones or use the bounding-box centre.}
\item{hemisphere_policy}{Error on an equator crossing or use the bounding-box
centre.}
\item{x}{A CRS-choice object.}
\item{...}{Reserved for method compatibility.}
}
\value{An object containing an \code{sf::crs}, EPSG code, UTM zone,
hemisphere and audit metadata.}
\description{Selects a WGS84 UTM CRS after checking latitude support,
antimeridian crossing, standard Norway/Svalbard exceptions, zone boundaries
and hemisphere boundaries.}
