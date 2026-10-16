\name{matern_cor}
\alias{matern_cor}
\title{Matérn Correlation}
\usage{
matern_cor(distance, scale, smoothness)
}
\arguments{
  \item{distance}{A numeric vector, matrix, or array of finite,
  non-negative distances.}
  \item{scale}{A finite positive numeric scalar in the same units as
  \code{distance}.}
  \item{smoothness}{A finite positive numeric scalar.}
}
\value{
A numeric object with the dimensions and names of \code{distance}, containing
correlations between zero and one.
}
\description{
Evaluates the unit-variance Matérn correlation with dimensionless argument
\code{distance / scale}. This scale convention is explicit; the argument is
not described as a practical or effective range.
}
\details{
At zero distance the function returns one by continuity. Positive arguments
are evaluated on the log scale using R's exponentially scaled modified Bessel
function. See the Matérn chapter in the package's technical derivations for
the assumptions, algorithm, and validation links.
}
\examples{
d <- as.matrix(dist(c(0, 0.5, 1)))
matern_cor(d, scale = 1, smoothness = 0.5)
}

