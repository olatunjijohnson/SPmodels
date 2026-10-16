#' Derivatives of Matern correlation with respect to scale
#'
#' Evaluate the first or second derivative of the SPmodels Matern correlation
#' with respect to its positive scale parameter. Smoothness and distance are
#' held fixed.
#'
#' @param distance A numeric vector, matrix, or array of finite non-negative
#'   distances.
#' @param scale A finite positive scale in the same units as `distance`.
#' @param smoothness A finite positive Matern smoothness.
#' @param order Derivative order, either 1 or 2.
#'
#' @return A numeric object with the dimensions and names of `distance`.
#' @export
matern_scale_derivative <- function(distance, scale, smoothness,
                                    order = 1L) {
  # Reuse the public correlation function for complete input validation.
  correlation <- matern_cor(distance, scale, smoothness)
  if (!is.numeric(order) || length(order) != 1L || is.na(order) ||
      !order %in% c(1, 2)) {
    stop("`order` must be 1 or 2.", call. = FALSE)
  }
  x <- distance / scale
  answer <- correlation
  answer[] <- 0
  positive <- x > sqrt(.Machine$double.eps)
  if (any(positive)) {
    xp <- x[positive]
    k_nu <- besselK(xp, nu = smoothness, expon.scaled = TRUE)
    k_previous <- besselK(xp, nu = smoothness - 1,
                          expon.scaled = TRUE)
    ratio <- k_previous / k_nu
    if (order == 1L) {
      value <- correlation[positive] * xp * ratio / scale
    } else {
      value <- correlation[positive] / scale^2 *
        (xp^2 - (2 * smoothness + 1) * xp * ratio)
    }
    if (any(!is.finite(value))) {
      stop("Matern scale derivative could not be evaluated at the supplied values.",
           call. = FALSE)
    }
    answer[positive] <- value
  }
  answer
}
