\name{sp_spatial_folds}
\alias{sp_spatial_folds}
\alias{print.spmodels_spatial_folds}
\title{Construct Geometry-Driven Spatial Folds}
\usage{
sp_spatial_folds(object, folds = 5L, seed = NULL, nstart = 20L)

\method{print}{spmodels_spatial_folds}(x, ...)
}
\arguments{
  \item{object}{A fitted spatial \code{spmodels_fit}.}
  \item{folds}{Integer number of spatial folds.}
  \item{seed}{Optional integer seed local to fold construction.}
  \item{nstart}{Positive number of random k-means starts.}
  \item{x}{An object returned by \code{sp_spatial_folds()}.}
  \item{...}{Reserved for print-method compatibility.}
}
\value{An integer vector in fitted-row order, with class
\code{spmodels_spatial_folds} and canonical cluster centres attached.}
\description{Clusters unique fitted coordinates into spatially compact folds.
Repeated coordinate locations always receive the same fold, preventing a
physical location from appearing in both training and held-out data. Cluster
labels are canonicalized from top to bottom and then left to right. The result
can be supplied directly to \code{sp_cross_validate()}.}
\details{This is k-means blocking in the fitted projected or declared
Euclidean coordinate system. It is not a guarantee of equal fold sizes or a
universal design for every scientific dependence structure.}
