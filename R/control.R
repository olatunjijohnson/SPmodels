#' Control Gaussian spatial model fitting
#'
#' Set covariance values, optimisation settings, and the units attached to
#' Euclidean coordinate distances. Matérn smoothness is fixed in this first
#' fitting implementation. The remaining covariance parameters are either all
#' estimated by maximum likelihood or all supplied by the user.
#'
#' @param smoothness Fixed positive Matérn smoothness.
#' @param spatial_variance,scale,nugget_variance Optional positive starting
#'   values when `estimate = TRUE`. Fixed spatial fitting requires all three;
#'   fixed non-spatial fitting requires only `nugget_variance`.
#' @param random_variance Optional positive scalar or vector of starting/fixed
#'   random-intercept variances. Named values must match grouping variables.
#' @param estimate Logical scalar. Estimate covariance parameters by maximum
#'   likelihood?
#' @param distance_units Non-empty character scalar for model distances and
#'   covariance scales, for example `"km"`. Numeric coordinate columns are
#'   assumed to already use these units. For projected `sf` input with a
#'   recognised linear CRS unit, coordinates are converted from the CRS unit
#'   to this unit before fitting.
#' @param maxit Positive integer maximum number of optimisation iterations.
#' @param reltol Positive optimisation tolerance.
#'
#' @return An object of class `spmodels_control`.
#' @export
sp_control <- function(smoothness = 0.5, spatial_variance = NULL,
                       scale = NULL, nugget_variance = NULL,
                       random_variance = NULL,
                       estimate = TRUE, distance_units = NULL,
                       maxit = 500L, reltol = 1e-8) {
  positive_scalar <- function(x, name, allow_null = TRUE) {
    if (is.null(x) && allow_null) return(NULL)
    if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
        !is.finite(x) || x <= 0) {
      stop(sprintf("`%s` must be one finite positive number.", name),
           call. = FALSE)
    }
    as.numeric(x)
  }
  smoothness <- positive_scalar(smoothness, "smoothness", FALSE)
  spatial_variance <- positive_scalar(spatial_variance, "spatial_variance")
  scale <- positive_scalar(scale, "scale")
  nugget_variance <- positive_scalar(nugget_variance, "nugget_variance")
  if (!is.null(random_variance) &&
      (!is.numeric(random_variance) || !length(random_variance) ||
       anyNA(random_variance) || any(!is.finite(random_variance)) ||
       any(random_variance <= 0))) {
    stop("`random_variance` must be NULL or finite positive values.",
         call. = FALSE)
  }
  if (!is.null(random_variance) && !is.null(names(random_variance)) &&
      (anyNA(names(random_variance)) || any(names(random_variance) == "") ||
       anyDuplicated(names(random_variance)))) {
    stop("Names on `random_variance` must be non-empty and unique.",
         call. = FALSE)
  }
  if (!is.logical(estimate) || length(estimate) != 1L || is.na(estimate)) {
    stop("`estimate` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(distance_units) &&
      (!is.character(distance_units) || length(distance_units) != 1L ||
       is.na(distance_units) || !nzchar(trimws(distance_units)))) {
    stop("`distance_units` must be NULL or one non-empty character value.",
         call. = FALSE)
  }
  if (!is.numeric(maxit) || length(maxit) != 1L || is.na(maxit) ||
      !is.finite(maxit) || maxit < 1 || maxit != as.integer(maxit)) {
    stop("`maxit` must be one positive integer.", call. = FALSE)
  }
  reltol <- positive_scalar(reltol, "reltol", FALSE)
  structure(list(
    smoothness = smoothness,
    spatial_variance = spatial_variance,
    scale = scale,
    nugget_variance = nugget_variance,
    random_variance = random_variance,
    estimate = estimate,
    distance_units = if (is.null(distance_units)) NULL else trimws(distance_units),
    maxit = as.integer(maxit), reltol = reltol
  ), class = "spmodels_control")
}
