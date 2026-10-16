#' Summarise a decay-adjusted intervention effect over time
#'
#' Calculate the remaining proportional reduction and corresponding
#' one-intervention prevalence multiplier at non-negative lags. When supplied
#' an `spmodels_fit` with estimated intervention parameters and usable
#' covariance, pointwise delta-method intervals are added.
#'
#' @param object An intervention specification or fitted model containing one.
#' @param lag Finite non-negative lags in the intervention time units.
#' @param level Pointwise confidence level strictly between zero and one.
#'
#' @return An object of class `spmodels_intervention_curve` containing an
#'   unrounded `values` data frame and intervention metadata.
#' @export
sp_intervention_curve <- function(object, lag, level = 0.95) {
  fit <- if (inherits(object, "spmodels_fit")) object else NULL
  intervention <- if (is.null(fit)) object else fit$intervention
  if (!inherits(intervention, "spmodels_intervention")) {
    stop("`object` must be an intervention specification or a fitted model containing one.",
         call. = FALSE)
  }
  if (!is.numeric(lag) || !length(lag) || anyNA(lag) ||
      any(!is.finite(lag)) || any(lag < 0)) {
    stop("`lag` must contain finite non-negative values.", call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      !is.finite(level) || level <= 0 || level >= 1) {
    stop("`level` must be one number strictly between zero and one.",
         call. = FALSE)
  }
  reduction <- intervention$alpha * exp(-lag / intervention$decay)
  standard_error <- lower <- upper <- rep(NA_real_, length(lag))
  if (!is.null(fit) && isTRUE(fit$intervention_estimated) &&
      !is.null(fit$intervention_vcov) &&
      all(is.finite(fit$intervention_vcov))) {
    gradient <- cbind(
      exp(-lag / intervention$decay),
      reduction * lag / intervention$decay^2
    )
    variance <- rowSums((gradient %*% fit$intervention_vcov) * gradient)
    standard_error <- sqrt(pmax(0, variance))
    z <- stats::qnorm((1 + level) / 2)
    lower <- pmax(0, reduction - z * standard_error)
    upper <- pmin(1, reduction + z * standard_error)
  }
  structure(list(
    values = data.frame(
      lag = as.numeric(lag), reduction = reduction,
      multiplier = 1 - reduction, std_error = standard_error,
      lower_reduction = lower, upper_reduction = upper,
      lower_multiplier = 1 - upper, upper_multiplier = 1 - lower
    ),
    alpha = intervention$alpha, decay = intervention$decay,
    time_units = intervention$time_units, level = level,
    estimated = !is.null(fit) && isTRUE(fit$intervention_estimated),
    call = match.call()
  ), class = "spmodels_intervention_curve")
}

#' @export
print.spmodels_intervention_curve <- function(x, ...) {
  cat("SPmodels intervention-effect curve\n")
  cat("Alpha:", format(x$alpha), " Decay:", format(x$decay),
      x$time_units, "\n")
  print(x$values, row.names = FALSE)
  invisible(x)
}

#' @export
plot.spmodels_intervention_curve <- function(x, ...) {
  values <- x$values
  graphics::plot(values$lag, values$reduction, type = "l",
                 xlab = paste("Lag (", x$time_units, ")", sep = ""),
                 ylab = "Remaining proportional reduction", ylim = c(0, 1),
                 ...)
  if (all(is.finite(values$lower_reduction)) && nrow(values) > 1L) {
    graphics::lines(values$lag, values$lower_reduction, lty = 2)
    graphics::lines(values$lag, values$upper_reduction, lty = 2)
  }
  invisible(x)
}
