#' Matern correlation
#'
#' Evaluate the unit-variance Matern correlation using the dimensionless
#' distance `distance / scale`. The scale is positive and has the same units as
#' `distance`; `smoothness` is the positive, dimensionless Matern smoothness.
#'
#' @param distance A numeric vector, matrix, or array of finite, non-negative
#'   distances.
#' @param scale A finite positive numeric scalar in the same units as
#'   `distance`.
#' @param smoothness A finite positive numeric scalar.
#'
#' @return A numeric object with the dimensions and names of `distance`, whose
#'   entries are correlations in `[0, 1]`.
#' @export
matern_cor <- function(distance, scale, smoothness) {
  if (!is.numeric(distance)) {
    stop("`distance` must be numeric.", call. = FALSE)
  }
  if (anyNA(distance) || any(!is.finite(distance))) {
    stop("`distance` must contain only finite, non-missing values.", call. = FALSE)
  }
  if (any(distance < 0)) {
    stop("`distance` must be non-negative.", call. = FALSE)
  }
  if (!is.numeric(scale) || length(scale) != 1L || is.na(scale) ||
      !is.finite(scale) || scale <= 0) {
    stop("`scale` must be one finite positive number.", call. = FALSE)
  }
  if (!is.numeric(smoothness) || length(smoothness) != 1L ||
      is.na(smoothness) || !is.finite(smoothness) || smoothness <= 0) {
    stop("`smoothness` must be one finite positive number.", call. = FALSE)
  }

  x <- distance / scale
  answer <- x
  is_zero <- x == 0
  answer[is_zero] <- 1

  positive <- !is_zero
  if (any(positive)) {
    xp <- x[positive]
    value <- if (smoothness == 0.5) {
      exp(-xp)
    } else if (smoothness == 1.5) {
      (1 + xp) * exp(-xp)
    } else if (smoothness == 2.5) {
      (1 + xp + xp^2 / 3) * exp(-xp)
    } else {
      scaled_bessel <- besselK(xp, nu = smoothness, expon.scaled = TRUE)
      log_correlation <- (1 - smoothness) * log(2) - lgamma(smoothness) +
        smoothness * log(xp) + log(scaled_bessel) - xp
      exp(log_correlation)
    }

    # At extremely small positive arguments besselK() can overflow although
    # the limiting correlation is numerically indistinguishable from one.
    recoverable <- !is.finite(value) & xp < .Machine$double.eps
    value[recoverable] <- 1
    if (any(!is.finite(value))) {
      stop("Matern correlation could not be evaluated at the supplied values.",
           call. = FALSE)
    }

    # Round-off can move a theoretical correlation a few ulps outside its
    # mathematical bounds.
    answer[positive] <- pmin(1, pmax(0, value))
  }

  answer
}

# Internal covariance paths have already validated non-negative finite
# distances and positive scalar parameters. Avoid repeating full matrix scans
# inside every optimizer evaluation for the common closed-form smoothnesses.
.sp_matern_cor_unchecked <- function(distance, scale, smoothness) {
  x <- distance / scale
  if (smoothness == 0.5) return(exp(-x))
  if (smoothness == 1.5) return((1 + x) * exp(-x))
  if (smoothness == 2.5) return((1 + x + x^2 / 3) * exp(-x))
  matern_cor(distance, scale, smoothness)
}

