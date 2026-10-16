#' Structured model-estimate table
#'
#' Convert an SPmodels fit into a stable, machine-readable table without
#' requiring downstream code to inspect internal fit fields.
#'
#' @param object A fitted `spmodels_fit` object.
#' @param include_smoothness Include the fixed Matérn smoothness in the table?
#'
#' @return A data frame with one row per mean or covariance parameter and
#'   columns `component`, `term`, `estimate`, `std_error`, `statistic`,
#'   `p_value`, `status`, and `units`.
#' @export
sp_model_table <- function(object, include_smoothness = TRUE) {
  if (!inherits(object, "spmodels_fit") &&
      !inherits(object, "spmodels_neighbor_fit")) {
    stop("`object` must be a fitted SPmodels model.", call. = FALSE)
  }
  if (!is.logical(include_smoothness) || length(include_smoothness) != 1L ||
      is.na(include_smoothness)) {
    stop("`include_smoothness` must be TRUE or FALSE.", call. = FALSE)
  }
  coefficient_se <- sqrt(diag(object$coefficient_vcov))
  coefficient_statistic <- object$coefficients / coefficient_se
  mean_table <- data.frame(
    component = "mean", term = names(object$coefficients),
    estimate = unname(object$coefficients),
    std_error = unname(coefficient_se),
    statistic = unname(coefficient_statistic),
    p_value = 2 * stats::pnorm(-abs(unname(coefficient_statistic))),
    status = "estimated", units = NA_character_,
    stringsAsFactors = FALSE
  )
  covariance_names <- if (isTRUE(object$spatial)) names(object$covariance) else
    "nugget_variance"
  covariance_values <- object$covariance[covariance_names]
  covariance_se <- sqrt(diag(object$covariance_vcov))[covariance_names]
  covariance_status <- if (isTRUE(object$covariance_estimated))
    "estimated" else "fixed"
  covariance_units <- rep(NA_character_, length(covariance_names))
  covariance_units[covariance_names == "scale"] <- object$distance_units
  covariance_table <- data.frame(
    component = "covariance", term = covariance_names,
    estimate = unname(covariance_values),
    std_error = unname(covariance_se), statistic = NA_real_,
    p_value = NA_real_, status = covariance_status,
    units = covariance_units, stringsAsFactors = FALSE
  )
  if (include_smoothness && isTRUE(object$spatial)) {
    covariance_table <- rbind(covariance_table, data.frame(
      component = "covariance", term = "smoothness",
      estimate = unname(object$smoothness), std_error = NA_real_,
      statistic = NA_real_, p_value = NA_real_, status = "fixed",
      units = NA_character_, stringsAsFactors = FALSE
    ))
  }
  if (length(object$random_variance)) {
    random_se <- sqrt(diag(object$random_variance_vcov))
    covariance_table <- rbind(covariance_table, data.frame(
      component = "random", term = names(object$random_variance),
      estimate = unname(object$random_variance),
      std_error = unname(random_se), statistic = NA_real_,
      p_value = NA_real_,
      status = if (isTRUE(object$covariance_estimated)) "estimated" else "fixed",
      units = NA_character_, stringsAsFactors = FALSE
    ))
  }
  random_slope_values <- .sp_random_slope_parameters(
    object$random_slope_specification
  )
  if (length(random_slope_values)) {
    random_slope_se <- if (is.null(object$random_slope_vcov)) {
      rep(NA_real_, length(random_slope_values))
    } else sqrt(diag(object$random_slope_vcov))
    slope_status <- unlist(lapply(object$random_slope_specification,
      function(component) {
        rep(if (isTRUE(component$estimate)) "estimated" else "fixed",
            length(.sp_random_slope_parameters(list(component))))
      }), use.names = FALSE)
    covariance_table <- rbind(covariance_table, data.frame(
      component = "random_slope", term = names(random_slope_values),
      estimate = unname(random_slope_values),
      std_error = unname(random_slope_se),
      statistic = NA_real_, p_value = NA_real_, status = slope_status,
      units = NA_character_, stringsAsFactors = FALSE
    ))
  }
  if (!is.null(object$temporal)) {
    covariance_table <- rbind(covariance_table, data.frame(
      component = "temporal", term = "scale",
      estimate = object$temporal_scale,
      std_error = if (isTRUE(object$temporal_scale_estimated))
        object$temporal_scale_se else NA_real_,
      statistic = NA_real_, p_value = NA_real_,
      status = if (isTRUE(object$temporal_scale_estimated)) "estimated" else
        "fixed",
      units = object$temporal$time_units, stringsAsFactors = FALSE
    ))
  }
  if (!is.null(object$intervention)) {
    intervention_se <- if (isTRUE(object$intervention_estimated)) {
      sqrt(diag(object$intervention_vcov))
    } else c(NA_real_, NA_real_)
    covariance_table <- rbind(covariance_table, data.frame(
      component = "intervention", term = c("alpha", "decay"),
      estimate = c(object$intervention$alpha, object$intervention$decay),
      std_error = unname(intervention_se), statistic = NA_real_,
      p_value = NA_real_,
      status = if (isTRUE(object$intervention_estimated)) "estimated" else
        "fixed",
      units = c(NA_character_, object$intervention$time_units),
      stringsAsFactors = FALSE
    ))
  }
  row.names(mean_table) <- row.names(covariance_table) <- NULL
  rbind(mean_table, covariance_table)
}

#' Structured prediction-summary table
#'
#' Summarise joint prediction draws at each location, or return the existing
#' grouped-target summary. The calculations preserve the draw dependence used
#' by [sp_predict_target()].
#'
#' @param object An object returned by [sp_predict_draws()],
#'   [sp_predict_draws_neighbor()], [sp_predict_target()] or
#'   [sp_predict_target_neighbor()].
#' @param threshold Optional finite exceedance threshold. It may only be
#'   supplied for prediction draws; target objects retain their own threshold.
#' @param probs Numeric quantile probabilities between zero and one. It is used
#'   only for prediction draws.
#'
#' @return A data frame containing identifiers, means, standard deviations,
#'   requested quantiles, optional exceedance probabilities, and prediction
#'   metadata.
#' @export
sp_prediction_table <- function(object, threshold = NULL,
                                probs = c(0.025, 0.5, 0.975)) {
  if (inherits(object, "spmodels_prediction_draws")) {
    target_object <- sp_predict_target(object, threshold = threshold,
                                       probs = probs)
  } else if (inherits(object, "spmodels_prediction_target")) {
    if (!is.null(threshold)) {
      stop("`threshold` is already defined by a prediction-target object.",
           call. = FALSE)
    }
    target_object <- object
  } else {
    stop("`object` must be a fitted prediction-draw or prediction-target result.",
         call. = FALSE)
  }
  out <- target_object$summary
  out$family <- target_object$family
  out$prediction_target <- target_object$target
  out$scale <- target_object$type
  row.names(out) <- NULL
  out
}

.sp_metric_table <- function(data, assessment, scope, id = NULL) {
  if (!is.data.frame(data)) data <- as.data.frame(data)
  if (is.null(id)) {
    id <- rep("all", nrow(data))
  } else {
    id <- as.character(id)
  }
  metric_columns <- names(data)[vapply(data, is.numeric, logical(1))]
  if (!length(metric_columns)) {
    stop("The validation object contains no numeric metrics.", call. = FALSE)
  }
  do.call(rbind, lapply(seq_len(nrow(data)), function(i) {
    data.frame(
      assessment = assessment, scope = scope, id = id[i],
      metric = metric_columns,
      value = vapply(data[metric_columns], function(column) column[i],
                     numeric(1)),
      stringsAsFactors = FALSE, row.names = NULL
    )
  }))
}

#' Structured validation table
#'
#' Return a common long-form metric schema for cross-validation, predictive-
#' draw assessment, or repeated simulation recovery.
#'
#' @param object An SPmodels cross-validation, simulation-assessment, or
#'   simulation-recovery object.
#' @param detail For cross-validation, return overall or fold metrics. Other
#'   supported objects currently expose their natural summary level.
#'
#' @return A data frame with columns `assessment`, `scope`, `id`, `metric`, and
#'   `value`.
#' @export
sp_validation_table <- function(object, detail = c("summary", "fold")) {
  detail <- match.arg(detail)
  if (inherits(object, "spmodels_cross_validation")) {
    if (detail == "fold") {
      return(.sp_metric_table(
        object$fold_summary[setdiff(names(object$fold_summary), "fold")],
        "cross_validation", "fold", object$fold_summary$fold
      ))
    }
    return(.sp_metric_table(object$summary, "cross_validation", "overall"))
  }
  if (detail == "fold") {
    stop("`detail = \"fold\"` is only available for cross-validation objects.",
         call. = FALSE)
  }
  if (inherits(object, "spmodels_simulation_assessment")) {
    return(.sp_metric_table(object$summary, "predictive_assessment", "overall"))
  }
  if (inherits(object, "spmodels_simulation_recovery")) {
    metrics <- object$summary[setdiff(names(object$summary), "parameter")]
    return(.sp_metric_table(metrics, "simulation_recovery", "parameter",
                            object$summary$parameter))
  }
  if (inherits(object, "spmodels_simulation_scenarios")) {
    scenario_order <- unique(object$summary$scenario)
    pieces <- lapply(scenario_order, function(scenario) {
      x <- object$summary[object$summary$scenario == scenario, , drop = FALSE]
      metrics <- x[setdiff(names(x), c("scenario", "parameter"))]
      table <- .sp_metric_table(metrics, "simulation_recovery", "parameter",
                                x$parameter)
      table$scope <- paste0("scenario:", x$scenario[[1L]])
      table
    })
    answer <- do.call(rbind, pieces)
    row.names(answer) <- NULL
    return(answer)
  }
  stop("`object` must be an SPmodels validation or simulation-assessment result.",
       call. = FALSE)
}
