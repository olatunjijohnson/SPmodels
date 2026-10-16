#' Specify a decay-adjusted intervention history
#'
#' Define an independently represented intervention history for a binomial
#' spatial model. At survey time `t`, every earlier intervention at time `u`
#' in the same unit multiplies baseline prevalence by
#' `1 - alpha * exp(-(t - u) / decay)`. An intervention at exactly the survey
#' time is excluded, corresponding to a survey immediately before rollout.
#'
#' @param history Data frame with one row per unit/intervention time.
#' @param unit Name of the unit-identifying column in both `history` and the
#'   survey/prediction data.
#' @param time Name of the intervention-time column in `history` and the
#'   survey-time column in survey/prediction data.
#' @param alpha Immediate proportional reduction, between zero and one.
#' @param decay Positive decay time in the same units as `time`.
#' @param time_units Character label for the units shared by `time` and
#'   `decay`, such as `"years"`.
#' @param estimate Estimate `alpha` and `decay` during Laplace fitting? Their
#'   supplied values are then starting values.
#' @param alpha_penalty Optional two positive beta-shape parameters used as a
#'   normalized penalty for `alpha` during estimation.
#' @param treated Optional name of a zero-one column in `history`. Rows with
#'   zero are retained for audit but contribute no intervention effect.
#'
#' @return An object of class `spmodels_intervention`.
#' @export
sp_intervention <- function(history, unit, time, alpha, decay, time_units,
                            estimate = FALSE, alpha_penalty = NULL,
                            treated = NULL) {
  if (!is.data.frame(history) || !nrow(history)) {
    stop("`history` must be a non-empty data frame.", call. = FALSE)
  }
  scalar_name <- function(value, name) {
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
        !nzchar(value) || !value %in% names(history)) {
      stop(sprintf("`%s` must name one column in `history`.", name),
           call. = FALSE)
    }
    value
  }
  unit <- scalar_name(unit, "unit")
  time <- scalar_name(time, "time")
  if (!is.null(treated)) treated <- scalar_name(treated, "treated")
  history_unit <- history[[unit]]
  if (anyNA(history_unit) || any(!nzchar(as.character(history_unit)))) {
    stop("Intervention unit identifiers must be non-missing and non-empty.",
         call. = FALSE)
  }
  history_time <- history[[time]]
  if (!is.numeric(history_time) || anyNA(history_time) ||
      any(!is.finite(history_time))) {
    stop("Intervention times must be finite numeric values.", call. = FALSE)
  }
  treatment <- if (is.null(treated)) rep(1, nrow(history)) else
    history[[treated]]
  if (!is.numeric(treatment) && !is.logical(treatment)) {
    stop("The intervention indicator must contain only zero and one.",
         call. = FALSE)
  }
  treatment <- as.numeric(treatment)
  if (anyNA(treatment) || any(!treatment %in% c(0, 1))) {
    stop("The intervention indicator must contain only zero and one.",
         call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      !is.finite(alpha) || alpha < 0 || alpha > 1) {
    stop("`alpha` must be one finite number between zero and one.",
         call. = FALSE)
  }
  if (!is.numeric(decay) || length(decay) != 1L || is.na(decay) ||
      !is.finite(decay) || decay <= 0) {
    stop("`decay` must be one finite positive number.", call. = FALSE)
  }
  if (!is.character(time_units) || length(time_units) != 1L ||
      is.na(time_units) || !nzchar(trimws(time_units))) {
    stop("`time_units` must be one non-empty character label.",
         call. = FALSE)
  }
  if (!is.logical(estimate) || length(estimate) != 1L || is.na(estimate)) {
    stop("`estimate` must be TRUE or FALSE.", call. = FALSE)
  }
  if (estimate && (alpha <= 0 || alpha >= 1)) {
    stop("Estimated `alpha` requires a starting value strictly between zero and one.",
         call. = FALSE)
  }
  if (!is.null(alpha_penalty)) {
    if (!is.numeric(alpha_penalty) || length(alpha_penalty) != 2L ||
        anyNA(alpha_penalty) || any(!is.finite(alpha_penalty)) ||
        any(alpha_penalty <= 0)) {
      stop("`alpha_penalty` must be NULL or two finite positive beta-shape parameters.",
           call. = FALSE)
    }
    alpha_penalty <- stats::setNames(as.numeric(alpha_penalty),
                                     c("shape1", "shape2"))
  }
  normalized <- data.frame(
    unit = as.character(history_unit), time = as.numeric(history_time),
    treated = treatment, source_row = seq_len(nrow(history)),
    stringsAsFactors = FALSE
  )
  if (anyDuplicated(normalized[c("unit", "time")])) {
    stop("`history` must have at most one row per unit and intervention time.",
         call. = FALSE)
  }
  normalized <- normalized[order(normalized$unit, normalized$time,
                                 normalized$source_row), , drop = FALSE]
  row.names(normalized) <- NULL
  structure(list(
    history = normalized, unit = unit, time = time,
    treated = treated, alpha = as.numeric(alpha), decay = as.numeric(decay),
    time_units = trimws(time_units),
    estimate = estimate, alpha_penalty = alpha_penalty,
    time_rule = "strictly before survey"
  ), class = "spmodels_intervention")
}

.sp_intervention_multiplier <- function(specification, data, rows = NULL) {
  if (is.null(specification)) return(NULL)
  if (!inherits(specification, "spmodels_intervention")) {
    stop("`intervention` must be created by `sp_intervention()`.",
         call. = FALSE)
  }
  model_data <- if (inherits(data, "sf")) sf::st_drop_geometry(data) else data
  required <- c(specification$unit, specification$time)
  if (!all(required %in% names(model_data))) {
    stop(sprintf("Intervention prediction data require columns: %s.",
                 paste(required, collapse = ", ")), call. = FALSE)
  }
  units <- model_data[[specification$unit]]
  times <- model_data[[specification$time]]
  if (!is.null(rows)) {
    units <- units[rows]
    times <- times[rows]
  }
  if (anyNA(units) || any(!nzchar(as.character(units)))) {
    stop("Survey unit identifiers must be non-missing and non-empty.",
         call. = FALSE)
  }
  if (!is.numeric(times) || anyNA(times) || any(!is.finite(times))) {
    stop("Survey times must be finite numeric values.", call. = FALSE)
  }
  history <- specification$history
  vapply(seq_along(times), function(index) {
    keep <- history$unit == as.character(units[[index]]) &
      history$time < times[[index]] & history$treated == 1
    if (!any(keep)) return(1)
    lag <- times[[index]] - history$time[keep]
    prod(1 - specification$alpha * exp(-lag / specification$decay))
  }, numeric(1L))
}

#' Calculate decay-adjusted intervention multipliers
#'
#' Evaluate the intervention-history multiplier for survey or prediction rows.
#'
#' @param intervention An object returned by [sp_intervention()].
#' @param data A data frame containing the declared unit and survey-time
#'   columns.
#'
#' @return One multiplier in `[0, 1]` per row of `data`.
#' @export
sp_intervention_multiplier <- function(intervention, data) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  .sp_intervention_multiplier(intervention, data)
}

#' @export
print.spmodels_intervention <- function(x, ...) {
  cat("SPmodels decay-adjusted intervention history\n")
  cat("Events:", nrow(x$history), " Units:",
      length(unique(x$history$unit)), "\n")
  cat("Immediate reduction alpha:", format(x$alpha),
      " Decay time:", format(x$decay), x$time_units, "\n")
  cat("Parameters:", if (x$estimate) "estimated" else "fixed", "\n")
  if (!is.null(x$alpha_penalty)) {
    cat("Alpha beta penalty:", paste(format(x$alpha_penalty), collapse = ", "),
        "\n")
  }
  cat("Timing rule:", x$time_rule, "\n")
  invisible(x)
}
