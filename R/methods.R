#' @export
coef.spmodels_fit <- function(object, ...) object$coefficients

#' @export
vcov.spmodels_fit <- function(object, ...) object$coefficient_vcov

#' @export
logLik.spmodels_fit <- function(object, ...) {
  value <- object$log_likelihood
  attr(value, "nobs") <- object$nobs
  attr(value, "df") <- object$rank +
    length(.sp_fitted_covariance_names(object)) +
    if (isTRUE(object$covariance_estimated))
      length(object$random_variance) else 0L
  attr(value, "df") <- attr(value, "df") +
    if (isTRUE(object$intervention_estimated)) 2L else 0L
  attr(value, "df") <- attr(value, "df") +
    if (isTRUE(object$temporal_scale_estimated)) 1L else 0L
  attr(value, "df") <- attr(value, "df") +
    length(if (is.null(object$random_slope_parameterization)) numeric() else
      object$random_slope_parameterization$values)
  class(value) <- "logLik"
  value
}

.sp_prediction_components <- function(object, newdata, unseen) {
  if (missing(newdata) || is.null(newdata)) {
    newdata <- object$model
    x_new <- object$x
    coordinates_new <- object$coordinates
    offset_new <- object$offset
    random_z_new <- if (is.null(object$random_specification)) NULL else
      object$random_specification$z
    random_slope_new <- object$random_slope_specification
    temporal_values_new <- object$temporal_values
  } else {
    if (!is.data.frame(newdata)) stop("`newdata` must be a data frame.", call. = FALSE)
    terms_new <- stats::delete.response(object$terms)
    model_newdata <- if (inherits(newdata, "sf")) {
      if (!requireNamespace("sf", quietly = TRUE)) {
        stop("Package `sf` is required for `sf` input.", call. = FALSE)
      }
      sf::st_drop_geometry(newdata)
    } else newdata
    mf_new <- stats::model.frame(terms_new, model_newdata, na.action = stats::na.fail,
                                 xlev = object$xlevels)
    x_new <- stats::model.matrix(terms_new, mf_new,
                                 contrasts.arg = object$contrasts)
    if (identical(object$coordinate_source, "sf") && !inherits(newdata, "sf")) {
      stop("A model fitted from `sf` data requires `sf` prediction data.",
           call. = FALSE)
    }
    if (isTRUE(object$spatial)) {
      coordinate_info <- .sp_coordinates(object$coords, newdata,
                                         expected_crs = object$crs)
      coordinate_multiplier <- object$coordinate_multiplier
      if (is.null(coordinate_multiplier)) coordinate_multiplier <- 1
      coordinates_new <- coordinate_info$values * coordinate_multiplier
      if (nrow(coordinates_new) != nrow(x_new)) {
        stop("Coordinate and predictor rows in `newdata` do not align.", call. = FALSE)
      }
    } else {
      coordinates_new <- NULL
    }
    offset_new <- stats::model.offset(mf_new)
    if (is.null(offset_new)) offset_new <- rep(0, nrow(x_new))
    random_z_new <- .sp_random_new(object$random_specification, newdata,
                                   unseen = unseen)
    random_slope_new <- .sp_random_slope_new(
      object$random_slope_specification, newdata, unseen = unseen
    )
    temporal_values_new <- .sp_temporal_values(object$temporal, newdata)
  }
  if (!isTRUE(object$spatial)) {
    fit <- as.vector(x_new %*% object$coefficients + offset_new)
    return(list(fit = fit, x = x_new, cross_covariance = NULL))
  }
  cross_distance <- .sp_dist(object$coordinates, coordinates_new)
  cross_correlation <- if (is.null(object$temporal)) {
    matern_cor(cross_distance, object$covariance[["scale"]], object$smoothness)
  } else {
    spacetime_cor(cross_distance,
      .sp_temporal_lag(object$temporal_values, temporal_values_new),
      object$covariance[["scale"]], object$smoothness,
      object$temporal_scale)
  }
  cross_covariance <- object$covariance[["spatial_variance"]] *
    cross_correlation
  if (!is.null(object$random_specification)) {
    cross_covariance <- cross_covariance + .sp_random_covariance(
      object$random_specification$z, random_z_new,
      object$random_specification$column_group, object$random_variance
    )
  }
  if (!is.null(object$random_slope_specification)) {
    cross_covariance <- cross_covariance + .sp_random_slope_covariance(
      object$random_slope_specification, random_slope_new
    )
  }
  alpha <- .sp_solve_chol(object$fitted_state$chol_v,
                          object$residual_mean)
  fit <- as.vector(x_new %*% object$coefficients + offset_new +
                     crossprod(cross_covariance, alpha))
  list(fit = fit, x = x_new, cross_covariance = cross_covariance,
       random_z = random_z_new, random_slope = random_slope_new,
       temporal_values = temporal_values_new)
}

.sp_logistic_normal_moments <- function(mean, variance) {
  if (variance <= 100 * .Machine$double.eps) {
    probability <- stats::plogis(mean)
    return(c(first = probability, second = probability^2))
  }
  standard_deviation <- sqrt(variance)
  first <- stats::integrate(
    function(z) stats::plogis(mean + standard_deviation * z) *
      stats::dnorm(z), -Inf, Inf, rel.tol = 1e-9
  )$value
  second <- stats::integrate(
    function(z) stats::plogis(mean + standard_deviation * z)^2 *
      stats::dnorm(z), -Inf, Inf, rel.tol = 1e-9
  )$value
  c(first = first, second = second)
}

.sp_predict_nongaussian <- function(object, newdata, se.fit, target, type,
                                    unseen) {
  if (is.null(newdata)) {
    x_new <- object$x
    coordinates_new <- object$coordinates
    offset_new <- object$offset
    probability_multiplier <- object$intervention_multiplier
    temporal_values_new <- object$temporal_values
    random_z_new <- if (is.null(object$random_specification)) NULL else
      object$random_specification$z
    random_slope_new <- object$random_slope_specification
  } else {
    if (!is.data.frame(newdata)) stop("`newdata` must be a data frame.", call. = FALSE)
    terms_new <- stats::delete.response(object$terms)
    model_newdata <- if (inherits(newdata, "sf")) {
      if (!requireNamespace("sf", quietly = TRUE)) {
        stop("Package `sf` is required for `sf` input.", call. = FALSE)
      }
      sf::st_drop_geometry(newdata)
    } else newdata
    mf_new <- stats::model.frame(terms_new, model_newdata,
      na.action = stats::na.fail, xlev = object$xlevels)
    x_new <- stats::model.matrix(terms_new, mf_new,
                                 contrasts.arg = object$contrasts)
    if (identical(object$coordinate_source, "sf") && !inherits(newdata, "sf")) {
      stop("A model fitted from `sf` data requires `sf` prediction data.",
           call. = FALSE)
    }
    coordinate_info <- .sp_coordinates(object$coords, newdata,
                                        expected_crs = object$crs)
    coordinate_multiplier <- object$coordinate_multiplier
    if (is.null(coordinate_multiplier)) coordinate_multiplier <- 1
    coordinates_new <- coordinate_info$values * coordinate_multiplier
    if (nrow(coordinates_new) != nrow(x_new)) {
      stop("Coordinate and predictor rows in `newdata` do not align.",
           call. = FALSE)
    }
    offset_new <- stats::model.offset(mf_new)
    if (is.null(offset_new)) offset_new <- rep(0, nrow(x_new))
    probability_multiplier <- .sp_intervention_multiplier(
      object$intervention, newdata
    )
    temporal_values_new <- .sp_temporal_values(object$temporal, newdata)
    random_z_new <- .sp_random_new(object$random_specification, newdata,
                                   unseen = unseen)
    random_slope_new <- .sp_random_slope_new(
      object$random_slope_specification, newdata, unseen = unseen
    )
  }
  cross_distance <- .sp_dist(object$coordinates, coordinates_new)
  cross_correlation <- if (is.null(object$temporal)) {
    matern_cor(cross_distance, object$covariance[["scale"]],
               object$smoothness)
  } else {
    spacetime_cor(cross_distance,
      .sp_temporal_lag(object$temporal_values, temporal_values_new),
      object$covariance[["scale"]], object$smoothness,
      object$temporal_scale)
  }
  cross_covariance <- object$covariance[["spatial_variance"]] *
    cross_correlation
  if (!is.null(object$random_specification)) {
    cross_covariance <- cross_covariance + .sp_random_covariance(
      object$random_specification$z, random_z_new,
      object$random_specification$column_group, object$random_variance
    )
  }
  if (!is.null(object$random_slope_specification)) {
    cross_covariance <- cross_covariance + .sp_random_slope_covariance(
      object$random_slope_specification, random_slope_new
    )
  }
  precision_cross <- .sp_solve_chol(object$fitted_state$chol_sigma,
                                    cross_covariance)
  latent_mean <- as.vector(x_new %*% object$coefficients + offset_new +
    crossprod(cross_covariance,
      .sp_solve_chol(object$fitted_state$chol_sigma,
                     object$latent_mode$mode)))
  prior_variance <- object$covariance[["spatial_variance"]] +
    if (target == "observation") object$covariance[["nugget_variance"]] else 0
  if (!is.null(object$random_specification)) {
    random_new_covariance <- .sp_random_new_covariance(
      random_z_new, object$random_specification, object$random_variance
    )
    prior_variance <- prior_variance + diag(random_new_covariance)
  }
  if (!is.null(object$random_slope_specification)) {
    prior_variance <- prior_variance + diag(
      .sp_random_slope_covariance(random_slope_new)
    )
  }
  conditional_variance <- prior_variance -
    colSums(cross_covariance * precision_cross)
  mode_uncertainty <- colSums(precision_cross *
    (object$latent_mode$covariance %*% precision_cross))
  latent_variance <- conditional_variance + mode_uncertainty
  tolerance <- 100 * .Machine$double.eps * max(1, prior_variance)
  if (any(latent_variance < -tolerance)) {
    warning("A negative predictive variance was produced; check conditioning.",
            call. = FALSE)
  }
  latent_variance <- pmax(0, latent_variance)
  if (type == "link") {
    if (!se.fit) return(latent_mean)
    return(list(fit = latent_mean, se.fit = sqrt(latent_variance),
                target = target, type = type))
  }
  if (object$family$family == "poisson") {
    fit <- exp(latent_mean + 0.5 * latent_variance)
    response_variance <- fit^2 * expm1(latent_variance)
    if (target == "observation") response_variance <- response_variance + fit
  } else {
    moments <- vapply(seq_along(latent_mean), function(i) {
      .sp_logistic_normal_moments(latent_mean[i], latent_variance[i])
    }, numeric(2))
    if (is.null(probability_multiplier)) {
      probability_multiplier <- rep(1, length(latent_mean))
    }
    fit <- probability_multiplier * moments["first", ]
    second <- probability_multiplier^2 * moments["second", ]
    response_variance <- pmax(0, second - fit^2)
    if (target == "observation") response_variance <- fit * (1 - fit)
  }
  if (!se.fit) return(unname(fit))
  list(fit = unname(fit), se.fit = sqrt(unname(response_variance)),
       target = target, type = type)
}

#' Predict from a fitted SPmodels model
#'
#' Return conditional Gaussian predictions or Laplace-approximated binomial and
#' Poisson point predictions at new coordinates. Non-Gaussian response means
#' integrate over approximate latent uncertainty rather than merely applying
#' the inverse link to the predictive linear-predictor mean.
#'
#' @param object A fitted `spmodels_fit` object.
#' @param newdata Optional data frame containing predictors and coordinates.
#' @param se.fit Return a list containing predictions and standard errors?
#' @param target Predict the latent spatial process or a future noisy
#'   observation. Observation uncertainty includes the nugget variance.
#' @param type Return predictions on the response or linear-predictor scale.
#'   These scales coincide for the Gaussian identity-link model.
#' @param unseen Policy for grouping levels absent during fitting. `"error"`
#'   rejects them. `"marginal"` integrates new zero-mean grouped coefficients
#'   with their declared covariance and preserves dependence between rows
#'   sharing the same new level.
#' @param ... Reserved for method compatibility.
#'
#' @return A numeric prediction vector, or a list with `fit`, `se.fit`,
#'   `target`, and, for non-Gaussian fits, `type`.
#' @export
predict.spmodels_fit <- function(object, newdata = NULL, se.fit = FALSE,
                                 target = c("latent", "observation"),
                                 type = c("response", "link"),
                                 unseen = c("error", "marginal"), ...) {
  target <- match.arg(target)
  type <- match.arg(type)
  unseen <- match.arg(unseen)
  if (!is.logical(se.fit) || length(se.fit) != 1L || is.na(se.fit)) {
    stop("`se.fit` must be TRUE or FALSE.", call. = FALSE)
  }
  if (object$family$family != "gaussian") {
    return(.sp_predict_nongaussian(object, newdata, se.fit, target, type,
                                   unseen))
  }
  pieces <- .sp_prediction_components(object, newdata, unseen)
  if (!se.fit) return(pieces$fit)
  if (!isTRUE(object$spatial)) {
    prediction_variance <- rowSums(
      (pieces$x %*% object$coefficient_vcov) * pieces$x
    )
    if (target == "observation") {
      prediction_variance <- prediction_variance +
        object$covariance[["nugget_variance"]]
    }
    return(list(fit = pieces$fit,
                se.fit = sqrt(pmax(0, unname(prediction_variance))),
                target = target))
  }
  vinv_c <- .sp_solve_chol(object$fitted_state$chol_v,
                           pieces$cross_covariance)
  process_reduction <- colSums(pieces$cross_covariance * vinv_c)
  mean_difference <- pieces$x -
    crossprod(pieces$cross_covariance, object$fitted_state$vinv_x)
  mean_variance <- rowSums((mean_difference %*% object$coefficient_vcov) *
                             mean_difference)
  prior_variance <- rep(object$covariance[["spatial_variance"]],
                        nrow(pieces$x))
  if (!is.null(object$random_specification)) {
    random_new_covariance <- .sp_random_new_covariance(
      pieces$random_z, object$random_specification, object$random_variance
    )
    prior_variance <- prior_variance + diag(random_new_covariance)
  }
  if (!is.null(object$random_slope_specification)) {
    prior_variance <- prior_variance + diag(
      .sp_random_slope_covariance(pieces$random_slope)
    )
  }
  prediction_variance <- prior_variance -
    process_reduction + mean_variance
  if (target == "observation") {
    prediction_variance <- prediction_variance +
      object$covariance[["nugget_variance"]]
  }
  tolerance <- 100 * .Machine$double.eps *
    max(1, object$covariance[["spatial_variance"]])
  if (any(prediction_variance < -tolerance)) {
    warning("A negative prediction variance was produced; check conditioning.",
            call. = FALSE)
  }
  prediction_variance <- unname(prediction_variance)
  list(fit = pieces$fit, se.fit = sqrt(pmax(0, prediction_variance)),
       target = target)
}

#' @export
fitted.spmodels_fit <- function(object, ...) {
  if (object$family$family != "gaussian") {
    return(as.vector(object$latent_mode$conditional_mean))
  }
  stats::predict(object)
}

#' @export
residuals.spmodels_fit <- function(object, ...) object$y - stats::fitted(object)

#' @export
print.spmodels_fit <- function(x, ...) {
  label <- if (x$family$family == "gaussian") "Gaussian" else
    paste0(tools::toTitleCase(x$family$family), "-", x$family$link)
  cat(if (isTRUE(x$spatial)) paste(label, "geostatistical regression\n") else
        paste("Independent", label, "regression\n"))
  cat("Formula: "); print(x$formula)
  cat("Observations:", x$nobs, "\n")
  if (isTRUE(x$spatial)) cat("Distance units:", x$distance_units, "\n")
  if (!is.null(x$temporal)) {
    cat("Temporal scale:", format(x$temporal_scale),
        x$temporal$time_units,
        if (isTRUE(x$temporal_scale_estimated)) "(estimated)" else "(fixed)",
        "\n")
  }
  method <- if (!is.null(x$likelihood_method)) x$likelihood_method else "ML"
  cat("Likelihood method:", method, "\n")
  covariance_names <- if (isTRUE(x$spatial)) names(x$covariance) else
    "nugget_variance"
  estimated_covariance <- .sp_fitted_covariance_names(x)
  covariance_status <- ifelse(covariance_names %in% estimated_covariance,
                              "estimated", "fixed")
  cat("Covariance parameters:\n")
  if (isTRUE(x$spatial)) {
    print(c(x$covariance, smoothness = x$smoothness))
  } else {
    print(x$covariance["nugget_variance"])
  }
  cat("Status:", paste(paste0(covariance_names, "=", covariance_status),
                        collapse = ", "), "\n")
  if (length(x$random_variance)) {
    cat("Random-intercept variances:\n")
    print(x$random_variance)
  }
  random_slope_values <- .sp_random_slope_parameters(
    x$random_slope_specification
  )
  if (length(random_slope_values)) {
    cat("Grouped random-slope covariance parameters:\n")
    print(random_slope_values)
  }
  if (!is.null(x$intervention)) {
    cat(if (isTRUE(x$intervention_estimated))
      "Estimated intervention parameters:\n" else
      "Fixed intervention parameters:\n")
    print(c(alpha = x$intervention$alpha, decay = x$intervention$decay))
  }
  if (!is.null(x$mcml)) {
    if (!is.null(x$mcml$adaptation)) {
      cat("Adaptive pilot ESS:",
          format(x$mcml$adaptation$pilot_diagnostics$effective_sample_size),
          "of", x$mcml$control$pilot_nsim, "\n")
      cat("Adaptive moments:",
          if (isTRUE(x$mcml$adaptation$used)) "used" else "not used",
          "(", x$mcml$adaptation$reason, ")\n")
    }
    cat("Importance ESS:", format(x$mcml$effective_sample_size), "of",
        x$mcml$control$nsim, "\n")
    cat("Log-likelihood MCSE:", format(x$mcml$log_likelihood_mcse), "\n")
  }
  invisible(x)
}

#' @export
summary.spmodels_fit <- function(object, ...) {
  standard_error <- sqrt(diag(stats::vcov(object)))
  table <- cbind(Estimate = stats::coef(object), `Std. Error` = standard_error,
                 `z value` = stats::coef(object) / standard_error)
  covariance_standard_error <- sqrt(diag(object$covariance_vcov))
  covariance_table <- cbind(Estimate = object$covariance,
                            `Std. Error` = covariance_standard_error)
  if (!isTRUE(object$spatial)) {
    covariance_table <- covariance_table["nugget_variance", , drop = FALSE]
  }
  random_table <- NULL
  if (length(object$random_variance)) {
    random_table <- cbind(
      Estimate = object$random_variance,
      `Std. Error` = sqrt(diag(object$random_variance_vcov))
    )
  }
  random_slope_values <- .sp_random_slope_parameters(
    object$random_slope_specification
  )
  random_slope_table <- if (!length(random_slope_values)) NULL else cbind(
    Estimate = random_slope_values,
    `Std. Error` = if (is.null(object$random_slope_vcov)) NA_real_ else
      sqrt(diag(object$random_slope_vcov))
  )
  intervention_table <- NULL
  if (!is.null(object$intervention)) {
    intervention_values <- c(alpha = object$intervention$alpha,
                             decay = object$intervention$decay)
    intervention_table <- cbind(
      Estimate = intervention_values,
      `Std. Error` = if (isTRUE(object$intervention_estimated))
        sqrt(diag(object$intervention_vcov)) else c(NA_real_, NA_real_)
    )
  }
  temporal_table <- if (is.null(object$temporal)) NULL else cbind(
    Estimate = object$temporal_scale,
    `Std. Error` = if (isTRUE(object$temporal_scale_estimated))
      object$temporal_scale_se else NA_real_
  )
  if (!is.null(temporal_table)) rownames(temporal_table) <- "scale"
  structure(list(call = object$call, spatial = object$spatial,
                 coefficients = table,
                 covariance = object$covariance,
                 covariance_table = covariance_table,
                 random = object$random_variance,
                 random_table = random_table,
                 random_slope = random_slope_values,
                 random_slope_table = random_slope_table,
                 intervention = object$intervention,
                 intervention_table = intervention_table,
                 intervention_estimated = object$intervention_estimated,
                 temporal = object$temporal,
                 temporal_table = temporal_table,
                 temporal_scale_estimated = object$temporal_scale_estimated,
                 smoothness = object$smoothness,
                 covariance_estimated = object$covariance_estimated,
                 covariance_estimated_names =
                   .sp_fitted_covariance_names(object),
                 log_likelihood = stats::logLik(object), nobs = object$nobs,
                 distance_units = object$distance_units,
                 family = object$family,
                 likelihood_method = object$likelihood_method,
                 mcml_diagnostics = if (is.null(object$mcml)) NULL else
                   object$mcml[c("effective_sample_size",
                     "effective_sample_fraction", "maximum_weight",
                     "log_likelihood_mcse")],
                 mcml_pilot_diagnostics = if (is.null(object$mcml$adaptation))
                   NULL else object$mcml$adaptation$pilot_diagnostics,
                 mcml_adaptation_used = if (is.null(object$mcml$adaptation))
                   NULL else object$mcml$adaptation$used,
                 mcml_adaptation_reason = if (is.null(object$mcml$adaptation))
                   NULL else object$mcml$adaptation$reason,
                 convergence = if (is.null(object$optimisation)) NA_integer_ else
                   object$optimisation$convergence,
                 convergence_message = if (is.null(object$optimisation)) NULL else
                   object$optimisation$message),
            class = "summary.spmodels_fit")
}

#' @export
print.summary.spmodels_fit <- function(x, ...) {
  label <- if (x$family$family == "gaussian") "Gaussian" else
    paste0(tools::toTitleCase(x$family$family), "-", x$family$link)
  cat(if (isTRUE(x$spatial)) paste(label, "geostatistical regression\n\nCall:\n") else
        paste("Independent", label, "regression\n\nCall:\n"))
  print(x$call)
  cat("\nMean coefficients:\n")
  stats::printCoefmat(x$coefficients)
  method <- if (!is.null(x$likelihood_method)) x$likelihood_method else "ML"
  cat("\nLikelihood method:", method, "\n")
  cat("Covariance parameters:\n")
  print(x$covariance_table)
  if (!is.null(x$random_table)) {
    cat("Random-intercept variances:\n")
    print(x$random_table)
  }
  if (!is.null(x$random_slope_table)) {
    cat("Fixed grouped random-slope covariance parameters:\n")
    print(x$random_slope_table)
  }
  if (!is.null(x$intervention)) {
    cat(if (isTRUE(x$intervention_estimated))
      "Estimated intervention parameters:\n" else
      "Fixed intervention parameters:\n")
    print(x$intervention_table)
  }
  if (!is.null(x$temporal_table)) {
    cat(if (isTRUE(x$temporal_scale_estimated))
      "Estimated temporal scale:\n" else "Fixed temporal scale:\n")
    print(x$temporal_table)
    cat("Time units:", x$temporal$time_units, "\n")
  }
  if (isTRUE(x$spatial)) cat("Fixed smoothness:", format(x$smoothness), "\n")
  if (x$covariance_estimated && isTRUE(x$spatial)) {
    cat("Optimiser convergence code:", x$convergence, "\n")
  }
  cat("\nLog-likelihood:", format(as.numeric(x$log_likelihood)), "\n")
  if (!is.null(x$mcml_diagnostics)) {
    if (!is.null(x$mcml_pilot_diagnostics)) {
      cat("Adaptive pilot ESS:",
          format(x$mcml_pilot_diagnostics$effective_sample_size),
          "(fraction",
          format(x$mcml_pilot_diagnostics$effective_sample_fraction),
          ")\n")
      cat("Adaptive pilot maximum normalized weight:",
          format(x$mcml_pilot_diagnostics$maximum_weight), "\n")
      cat("Adaptive moments:",
          if (isTRUE(x$mcml_adaptation_used)) "used" else "not used",
          "(", x$mcml_adaptation_reason, ")\n")
    }
    cat("Importance ESS:",
        format(x$mcml_diagnostics$effective_sample_size),
        "(fraction", format(x$mcml_diagnostics$effective_sample_fraction),
        ")\n")
    cat("Maximum normalized weight:",
        format(x$mcml_diagnostics$maximum_weight), "\n")
    cat("Log-likelihood MCSE:",
        format(x$mcml_diagnostics$log_likelihood_mcse), "\n")
  }
  if (isTRUE(x$spatial)) cat("Distance units:", x$distance_units, "\n")
  invisible(x)
}
