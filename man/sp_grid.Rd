\name{sp_grid}
\alias{sp_grid}
\title{Construct a Regular Prediction Grid}
\usage{
sp_grid(data, cellsize, boundary = NULL)
}
\arguments{
  \item{data}{A non-empty projected \code{sf} object containing points or
  polygonal geometry.}
  \item{cellsize}{One positive grid spacing or two positive x/y spacings, in
  the coordinate units of \code{data}.}
  \item{boundary}{Optional projected polygon \code{sf} or \code{sfc} object
  using the same CRS as \code{data}.}
}
\value{A point \code{sf} object with stable top-to-bottom, left-to-right
\code{grid_id} values.}
\description{Creates regular point centres within a projected polygon. If
point data are supplied without a boundary, their convex hull defines the
grid extent. Polygon data use their union. Geographic coordinates and
zero-area boundaries fail explicitly.}
\examples{
if (requireNamespace("sf", quietly = TRUE)) {
  points <- sf::st_as_sf(data.frame(x = c(0, 2, 0), y = c(0, 0, 2)),
                         coords = c("x", "y"), crs = 3857)
  grid <- sp_grid(points, cellsize = 0.5)
  nrow(grid)
}
}
