#' Separable Matérn-exponential space-time correlation
#'
#' Calculate a stationary separable correlation from spatial distances and
#' absolute time lags. Spatial correlation is Matérn and temporal correlation
#' is exponential.
#'
#' @param distance Non-negative finite spatial distances.
#' @param time_lag Non-negative finite time lags with the same dimensions as
#'   `distance`.
#' @param spatial_scale Positive Matérn spatial scale.
#' @param smoothness Positive Matérn smoothness.
#' @param temporal_scale Positive exponential temporal scale in the same units
#'   as `time_lag`.
#'
#' @return Correlations with dimensions and names inherited from `distance`.
#' @export
spacetime_cor <- function(distance, time_lag, spatial_scale, smoothness,
                          temporal_scale) {
  if (!is.numeric(distance) || !is.numeric(time_lag) ||
      !identical(dim(distance), dim(time_lag)) ||
      length(distance) != length(time_lag)) {
    stop("`distance` and `time_lag` must be numeric with identical dimensions.",
         call. = FALSE)
  }
  if (anyNA(time_lag) || any(!is.finite(time_lag)) || any(time_lag < 0)) {
    stop("`time_lag` must contain finite non-negative values.", call. = FALSE)
  }
  if (!is.numeric(temporal_scale) || length(temporal_scale) != 1L ||
      is.na(temporal_scale) || !is.finite(temporal_scale) ||
      temporal_scale <= 0) {
    stop("`temporal_scale` must be one finite positive number.",
         call. = FALSE)
  }
  spatial <- matern_cor(distance, spatial_scale, smoothness)
  answer <- spatial * exp(-time_lag / temporal_scale)
  attributes(answer) <- attributes(spatial)
  answer
}

#' Specify survey time for a separable space-time model
#'
#' @param time One-sided formula selecting one finite numeric time column.
#' @param scale Positive exponential temporal scale. When estimated, this is
#'   the starting value.
#' @param time_units Non-empty label for the time and scale units.
#' @param estimate Estimate temporal scale jointly with the other fitted
#'   parameters?
#'
#' @return An object of class `spmodels_temporal`.
#' @export
sp_temporal <- function(time, scale, time_units, estimate = FALSE) {
  if (!inherits(time, "formula") || length(time) != 2L) {
    stop("`time` must be a one-sided formula selecting one column.",
         call. = FALSE)
  }
  if (length(all.vars(time)) != 1L) {
    stop("`time` must select exactly one time variable.", call. = FALSE)
  }
  if (!is.numeric(scale) || length(scale) != 1L || is.na(scale) ||
      !is.finite(scale) || scale <= 0) {
    stop("`scale` must be one finite positive number.", call. = FALSE)
  }
  if (!is.character(time_units) || length(time_units) != 1L ||
      is.na(time_units) || !nzchar(trimws(time_units))) {
    stop("`time_units` must be one non-empty character label.", call. = FALSE)
  }
  if (!is.logical(estimate) || length(estimate) != 1L || is.na(estimate)) {
    stop("`estimate` must be TRUE or FALSE.", call. = FALSE)
  }
  structure(list(time = time, variable = all.vars(time), scale = as.numeric(scale),
                 time_units = trimws(time_units), estimate = estimate),
            class = "spmodels_temporal")
}

.sp_temporal_values <- function(specification, data, rows = NULL) {
  if (is.null(specification)) return(NULL)
  if (!inherits(specification, "spmodels_temporal")) {
    stop("`time` must be created by `sp_temporal()`.", call. = FALSE)
  }
  model_data <- if (inherits(data, "sf")) sf::st_drop_geometry(data) else data
  frame <- stats::model.frame(specification$time, model_data,
                              na.action = stats::na.pass)
  values <- frame[[1L]]
  if (!is.null(rows)) values <- values[rows]
  if (!is.numeric(values) || anyNA(values) || any(!is.finite(values))) {
    stop("Survey times must be finite numeric values.", call. = FALSE)
  }
  as.numeric(values)
}

.sp_temporal_lag <- function(left, right = left) {
  if (is.null(left)) return(NULL)
  abs(outer(left, right, "-"))
}

#' @export
print.spmodels_temporal <- function(x, ...) {
  cat("SPmodels separable temporal specification\n")
  cat("Time:"); print(x$time)
  cat("Exponential scale:", format(x$scale), x$time_units,
      if (x$estimate) "(estimated)" else "(fixed)", "\n")
  invisible(x)
}
