#' Control Gaussian spatial model fitting
#'
#' Set covariance values, optimisation settings, and the units attached to
#' Euclidean coordinate distances. Matérn smoothness is fixed in this first
#' fitting implementation. Spatial variance, scale, and supported additional
#' variance components are estimated together when `estimate = TRUE`; the
#' nugget can instead be held fixed with `estimate_nugget = FALSE`.
#'
#' @param smoothness Fixed positive Matérn smoothness.
#' @param spatial_variance,scale Optional positive starting values when
#'   `estimate = TRUE`. Fixed spatial fitting requires both.
#' @param nugget_variance Optional non-negative nugget variance. It is a
#'   starting value when estimated and a fixed value when
#'   `estimate_nugget = FALSE`. Supplying zero fixes the nugget at zero unless
#'   `estimate_nugget = TRUE`, which is invalid because zero cannot be a
#'   log-scale starting value.
#' @param random_variance Optional positive scalar or vector of starting/fixed
#'   random-intercept variances. Named values must match grouping variables.
#' @param estimate Logical scalar. Estimate covariance parameters by maximum
#'   likelihood? This controls the spatial variance, scale, and supported
#'   additional variance components.
#' @param estimate_nugget Optional logical scalar controlling nugget
#'   estimation when `estimate = TRUE`. The default estimates a missing or
#'   positive nugget and fixes a supplied zero nugget. Set it to `FALSE` to
#'   fix any supplied non-negative nugget value.
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
                       estimate = TRUE, estimate_nugget = NULL,
                       distance_units = NULL,
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
  if (!is.null(nugget_variance) &&
      (!is.numeric(nugget_variance) || length(nugget_variance) != 1L ||
       is.na(nugget_variance) || !is.finite(nugget_variance) ||
       nugget_variance < 0)) {
    stop("`nugget_variance` must be NULL or one finite non-negative number.",
         call. = FALSE)
  }
  if (!is.null(nugget_variance)) nugget_variance <- as.numeric(nugget_variance)
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
  if (!is.null(estimate_nugget) &&
      (!is.logical(estimate_nugget) || length(estimate_nugget) != 1L ||
       is.na(estimate_nugget))) {
    stop("`estimate_nugget` must be NULL, TRUE, or FALSE.", call. = FALSE)
  }
  if (!estimate && isTRUE(estimate_nugget)) {
    stop("`estimate_nugget = TRUE` requires `estimate = TRUE`.", call. = FALSE)
  }
  if (estimate && is.null(estimate_nugget)) {
    estimate_nugget <- is.null(nugget_variance) || nugget_variance > 0
  }
  if (!estimate) estimate_nugget <- FALSE
  if (identical(estimate_nugget, FALSE) && estimate &&
      is.null(nugget_variance)) {
    stop("Fixing the nugget requires a supplied `nugget_variance`.",
         call. = FALSE)
  }
  if (isTRUE(estimate_nugget) && !is.null(nugget_variance) &&
      nugget_variance == 0) {
    stop("An estimated nugget requires a positive starting value; use `estimate_nugget = FALSE` to fix it at zero.",
         call. = FALSE)
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
    estimate_nugget = isTRUE(estimate_nugget),
    distance_units = if (is.null(distance_units)) NULL else trimws(distance_units),
    maxit = as.integer(maxit), reltol = reltol
  ), class = "spmodels_control")
}

.sp_estimated_covariance_names <- function(control, spatial = TRUE) {
  if (!isTRUE(control$estimate)) return(character())
  if (!isTRUE(spatial)) {
    return(if (isTRUE(control$estimate_nugget)) "nugget_variance" else
      character())
  }
  c("spatial_variance", "scale",
    if (isTRUE(control$estimate_nugget)) "nugget_variance" else character())
}

.sp_fitted_covariance_names <- function(object) {
  if (!is.null(object$covariance_estimated_names)) {
    return(object$covariance_estimated_names)
  }
  if (!isTRUE(object$covariance_estimated)) return(character())
  if (isTRUE(object$spatial)) names(object$covariance) else "nugget_variance"
}
