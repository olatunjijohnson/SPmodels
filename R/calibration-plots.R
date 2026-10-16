#' Assess predictive interval calibration
#'
#' Compare central predictive intervals from joint prediction draws with known
#' or subsequently observed values over a declared set of nominal levels.
#'
#' @param object An object returned by [sp_predict_draws()].
#' @param truth One finite value for every prediction location.
#' @param levels Numeric central interval probabilities strictly between zero
#'   and one.
#' @param confidence Confidence level for Wilson uncertainty bounds around the
#'   empirical coverage proportions.
#'
#' @return An object of class `spmodels_calibration` containing the unrounded
#'   calibration curve and its pointwise Wilson bounds.
#' @export
sp_calibration <- function(object, truth,
                           levels = c(0.5, 0.8, 0.9, 0.95),
                           confidence = 0.95) {
  if (!inherits(object, "spmodels_prediction_draws")) {
    stop("`object` must be returned by `sp_predict_draws()`.", call. = FALSE)
  }
  if (!is.numeric(truth) || length(truth) != ncol(object$draws) ||
      anyNA(truth) || any(!is.finite(truth))) {
    stop("`truth` must contain one finite value per prediction location.",
         call. = FALSE)
  }
  if (!is.numeric(levels) || !length(levels) || anyNA(levels) ||
      any(!is.finite(levels)) || any(levels <= 0 | levels >= 1)) {
    stop("`levels` must contain probabilities strictly between zero and one.",
         call. = FALSE)
  }
  if (anyDuplicated(levels)) {
    stop("`levels` must not contain duplicates.", call. = FALSE)
  }
  if (!is.numeric(confidence) || length(confidence) != 1L ||
      is.na(confidence) || !is.finite(confidence) || confidence <= 0 ||
      confidence >= 1) {
    stop("`confidence` must be one probability strictly between zero and one.",
         call. = FALSE)
  }
  levels <- sort(as.numeric(levels))
  n <- length(truth)
  z <- stats::qnorm((1 + confidence) / 2)
  curve <- do.call(rbind, lapply(levels, function(nominal) {
    alpha <- (1 - nominal) / 2
    limits <- t(apply(object$draws, 2L, stats::quantile,
                      probs = c(alpha, 1 - alpha), names = FALSE, type = 8))
    covered <- truth >= limits[, 1L] & truth <= limits[, 2L]
    observed <- mean(covered)
    denominator <- 1 + z^2 / n
    centre <- (observed + z^2 / (2 * n)) / denominator
    half_width <- z * sqrt(observed * (1 - observed) / n +
                             z^2 / (4 * n^2)) / denominator
    data.frame(
      nominal = nominal, observed = observed,
      lower = max(0, centre - half_width),
      upper = min(1, centre + half_width),
      locations = n, stringsAsFactors = FALSE
    )
  }))
  structure(list(
    calibration = curve, confidence = confidence, truth = truth,
    family = object$family, target = object$target, type = object$type,
    call = match.call()
  ), class = "spmodels_calibration")
}

#' @export
print.spmodels_calibration <- function(x, ...) {
  cat("SPmodels predictive interval calibration\n")
  print(x$calibration, row.names = FALSE)
  invisible(x)
}

#' @export
plot.spmodels_calibration <- function(x, xlab = "Nominal coverage",
                                      ylab = "Observed coverage",
                                      main = "Predictive interval calibration",
                                      ...) {
  values <- x$calibration
  graphics::plot(values$nominal, values$observed, xlim = c(0, 1),
                 ylim = c(0, 1), xlab = xlab, ylab = ylab,
                 main = main, ...)
  graphics::segments(values$nominal, values$lower,
                     values$nominal, values$upper)
  graphics::abline(a = 0, b = 1, lty = 2)
  invisible(x)
}

#' Display cross-validation predictions or errors
#'
#' @param x An object returned by [sp_cross_validate()].
#' @param type Display observed against predicted values, errors against
#'   predicted values, or error distributions by fold.
#' @param xlab,ylab,main Plot labels.
#' @param ... Further arguments passed to the principal base-graphics call.
#' @export
plot.spmodels_cross_validation <- function(
    x, type = c("observed", "error", "fold"), xlab = NULL, ylab = NULL,
    main = NULL, ...) {
  type <- match.arg(type)
  values <- x$predictions
  if (type == "observed") {
    graphics::plot(values$predicted, values$observed,
      xlab = if (is.null(xlab)) "Predicted" else xlab,
      ylab = if (is.null(ylab)) "Observed" else ylab,
      main = if (is.null(main)) "Cross-validated predictions" else main, ...)
    graphics::abline(a = 0, b = 1, lty = 2)
  } else if (type == "error") {
    graphics::plot(values$predicted, values$error,
      xlab = if (is.null(xlab)) "Predicted" else xlab,
      ylab = if (is.null(ylab)) "Observed minus predicted" else ylab,
      main = if (is.null(main)) "Cross-validation errors" else main, ...)
    graphics::abline(h = 0, lty = 2)
  } else {
    graphics::boxplot(values$error ~ factor(values$fold),
      xlab = if (is.null(xlab)) "Fold" else xlab,
      ylab = if (is.null(ylab)) "Observed minus predicted" else ylab,
      main = if (is.null(main)) "Cross-validation errors by fold" else main,
      ...)
    graphics::abline(h = 0, lty = 2)
  }
  invisible(x)
}

#' Display prediction recovery from a draw assessment
#'
#' @param x An object returned by [sp_simulation_assess()].
#' @param type Display predicted against true values or errors against truth.
#' @param xlab,ylab,main Plot labels.
#' @param ... Further arguments passed to the principal base-graphics call.
#' @export
plot.spmodels_simulation_assessment <- function(
    x, type = c("prediction", "error"), xlab = NULL, ylab = NULL,
    main = NULL, ...) {
  type <- match.arg(type)
  values <- x$locations
  if (type == "prediction") {
    graphics::plot(values$truth, values$mean,
      xlab = if (is.null(xlab)) "Truth" else xlab,
      ylab = if (is.null(ylab)) "Predictive mean" else ylab,
      main = if (is.null(main)) "Prediction recovery" else main, ...)
    graphics::segments(values$truth, values$lower,
                       values$truth, values$upper)
    graphics::abline(a = 0, b = 1, lty = 2)
  } else {
    graphics::plot(values$truth, values$error,
      xlab = if (is.null(xlab)) "Truth" else xlab,
      ylab = if (is.null(ylab)) "Predictive mean minus truth" else ylab,
      main = if (is.null(main)) "Prediction errors" else main, ...)
    graphics::abline(h = 0, lty = 2)
  }
  invisible(x)
}

#' Display replicate estimates from simulation recovery
#'
#' @param x An object returned by [sp_simulation_recovery()].
#' @param parameter Name or index of one fitted parameter.
#' @param intervals Include reported Wald intervals when finite standard errors
#'   are available.
#' @param xlab,ylab,main Plot labels.
#' @param ... Further arguments passed to [graphics::plot()].
#' @export
plot.spmodels_simulation_recovery <- function(
    x, parameter = 1L, intervals = TRUE, xlab = "Simulation replicate",
    ylab = NULL, main = NULL, ...) {
  if (is.character(parameter)) {
    if (length(parameter) != 1L || !parameter %in% colnames(x$estimates)) {
      stop("`parameter` must name one recorded parameter.", call. = FALSE)
    }
    j <- match(parameter, colnames(x$estimates))
  } else {
    if (!is.numeric(parameter) || length(parameter) != 1L ||
        is.na(parameter) || parameter != as.integer(parameter) ||
        parameter < 1L || parameter > ncol(x$estimates)) {
      stop("`parameter` must index one recorded parameter.", call. = FALSE)
    }
    j <- as.integer(parameter)
  }
  parameter_name <- colnames(x$estimates)[j]
  estimate <- x$estimates[, j]
  valid <- is.finite(estimate)
  if (!any(valid)) stop("No successful estimate is available.", call. = FALSE)
  replicate <- seq_len(nrow(x$estimates))
  graphics::plot(replicate[valid], estimate[valid],
    xlab = xlab, ylab = if (is.null(ylab)) parameter_name else ylab,
    main = if (is.null(main)) "Parameter recovery" else main, ...)
  if (isTRUE(intervals)) {
    se <- x$standard_errors[, j]
    with_se <- valid & is.finite(se)
    z <- stats::qnorm((1 + x$level) / 2)
    graphics::segments(replicate[with_se], estimate[with_se] - z * se[with_se],
                       replicate[with_se], estimate[with_se] + z * se[with_se])
  }
  graphics::abline(h = x$truth[j], lty = 2)
  invisible(x)
}
