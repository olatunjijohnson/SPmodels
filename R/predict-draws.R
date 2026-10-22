.sp_prediction_design <- function(object, newdata, unseen = "error") {
  if (is.null(newdata)) {
    return(list(x = object$x, coordinates = object$coordinates,
                offset = object$offset, data = object$data_used,
                intervention_multiplier = object$intervention_multiplier,
                temporal_values = object$temporal_values,
                random_z = if (is.null(object$random_specification)) NULL else
                  object$random_specification$z,
                random_slope = object$random_slope_specification))
  }
  if (!is.data.frame(newdata)) {
    stop("`newdata` must be a data frame.", call. = FALSE)
  }
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
  offset_new <- stats::model.offset(mf_new)
  if (is.null(offset_new)) offset_new <- rep(0, nrow(x_new))
  if (!isTRUE(object$spatial)) {
    return(list(x = x_new, coordinates = NULL, offset = offset_new,
                random_z = NULL))
  }
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
  random_z <- .sp_random_new(object$random_specification, newdata,
                             unseen = unseen)
  random_slope <- .sp_random_slope_new(
    object$random_slope_specification, newdata, unseen = unseen
  )
  list(x = x_new, coordinates = coordinates_new, offset = offset_new,
       data = newdata,
       random_z = random_z,
       random_slope = random_slope,
       temporal_values = .sp_temporal_values(object$temporal, newdata),
       intervention_multiplier = .sp_intervention_multiplier(
         object$intervention, newdata
       ))
}

.sp_covariance_root <- function(covariance) {
  covariance <- (covariance + t(covariance)) / 2
  scale <- max(1, max(abs(diag(covariance))))
  tolerance <- 1000 * .Machine$double.eps * scale * nrow(covariance)
  chol_covariance <- tryCatch(chol(covariance), error = function(e) NULL)
  if (!is.null(chol_covariance)) {
    return(list(root = chol_covariance, orientation = "upper"))
  }
  decomposition <- eigen(covariance, symmetric = TRUE)
  if (min(decomposition$values) < -tolerance) {
    stop("The joint predictive covariance is not positive semidefinite; check conditioning.",
         call. = FALSE)
  }
  values <- pmax(0, decomposition$values)
  root <- sweep(decomposition$vectors, 2L, sqrt(values), "*")
  list(root = root, orientation = "lower")
}

.sp_prediction_random_slope <- function(object, design) {
  prediction <- design$random_slope
  fitted <- object$random_slope_specification
  if (is.null(prediction) || is.null(fitted)) return(prediction)
  if (length(prediction) != length(fitted) ||
      !identical(vapply(prediction, `[[`, character(1L), "label"),
                 vapply(fitted, `[[`, character(1L), "label"))) {
    stop("Fitted and prediction grouped random-slope components do not align.",
         call. = FALSE)
  }
  for (j in seq_along(prediction)) {
    prediction[[j]]$covariance <- fitted[[j]]$covariance
  }
  prediction
}

.sp_joint_gaussian <- function(object, design, target) {
  mean <- as.vector(design$x %*% object$coefficients + design$offset)
  if (!isTRUE(object$spatial)) {
    covariance <- design$x %*% object$coefficient_vcov %*% t(design$x)
    if (target == "observation") {
      diag(covariance) <- diag(covariance) +
        object$covariance[["nugget_variance"]]
    }
    return(list(mean = mean, covariance = covariance))
  }
  cross_distance <- .sp_dist(object$coordinates, design$coordinates)
  cross_correlation <- if (is.null(object$temporal)) {
    matern_cor(cross_distance, object$covariance[["scale"]],
               object$smoothness)
  } else {
    spacetime_cor(cross_distance,
      .sp_temporal_lag(object$temporal_values, design$temporal_values),
      object$covariance[["scale"]], object$smoothness,
      object$temporal_scale)
  }
  cross_covariance <- object$covariance[["spatial_variance"]] *
    cross_correlation
  prediction_random_slope <- .sp_prediction_random_slope(object, design)
  if (!is.null(object$random_specification)) {
    cross_covariance <- cross_covariance + .sp_random_covariance(
      object$random_specification$z, design$random_z,
      object$random_specification$column_group, object$random_variance
    )
  }
  if (!is.null(object$random_slope_specification)) {
    cross_covariance <- cross_covariance + .sp_random_slope_covariance(
      object$random_slope_specification, prediction_random_slope
    )
  }
  vinv_c <- .sp_solve_chol(object$fitted_state$chol_v, cross_covariance)
  mean <- mean + as.vector(crossprod(
    cross_covariance,
    .sp_solve_chol(object$fitted_state$chol_v, object$residual_mean)
  ))
  new_distance <- .sp_dist(design$coordinates)
  new_correlation <- if (is.null(object$temporal)) {
    matern_cor(new_distance, object$covariance[["scale"]], object$smoothness)
  } else {
    spacetime_cor(new_distance, .sp_temporal_lag(design$temporal_values),
      object$covariance[["scale"]], object$smoothness,
      object$temporal_scale)
  }
  covariance <- object$covariance[["spatial_variance"]] * new_correlation -
    crossprod(cross_covariance, vinv_c)
  if (!is.null(object$random_specification)) {
    covariance <- covariance + .sp_random_new_covariance(
      design$random_z, object$random_specification, object$random_variance
    )
  }
  if (!is.null(object$random_slope_specification)) {
    covariance <- covariance + .sp_random_slope_covariance(
      prediction_random_slope
    )
  }
  mean_difference <- design$x -
    crossprod(cross_covariance, object$fitted_state$vinv_x)
  covariance <- covariance + mean_difference %*%
    object$coefficient_vcov %*% t(mean_difference)
  if (target == "observation") {
    diag(covariance) <- diag(covariance) +
      object$covariance[["nugget_variance"]]
  }
  list(mean = mean, covariance = covariance)
}

.sp_joint_nongaussian <- function(object, design, target) {
  cross_distance <- .sp_dist(object$coordinates, design$coordinates)
  cross_correlation <- if (is.null(object$temporal)) {
    matern_cor(cross_distance, object$covariance[["scale"]],
               object$smoothness)
  } else {
    spacetime_cor(cross_distance,
      .sp_temporal_lag(object$temporal_values, design$temporal_values),
      object$covariance[["scale"]], object$smoothness,
      object$temporal_scale)
  }
  cross_covariance <- object$covariance[["spatial_variance"]] *
    cross_correlation
  prediction_random_slope <- .sp_prediction_random_slope(object, design)
  if (!is.null(object$random_specification)) {
    cross_covariance <- cross_covariance + .sp_random_covariance(
      object$random_specification$z, design$random_z,
      object$random_specification$column_group, object$random_variance
    )
  }
  if (!is.null(object$random_slope_specification)) {
    cross_covariance <- cross_covariance + .sp_random_slope_covariance(
      object$random_slope_specification, prediction_random_slope
    )
  }
  precision_cross <- .sp_solve_chol(object$fitted_state$chol_sigma,
                                    cross_covariance)
  mean <- as.vector(design$x %*% object$coefficients + design$offset +
    crossprod(cross_covariance,
      .sp_solve_chol(object$fitted_state$chol_sigma,
                     object$latent_mode$mode)))
  new_distance <- .sp_dist(design$coordinates)
  new_correlation <- if (is.null(object$temporal)) {
    matern_cor(new_distance, object$covariance[["scale"]], object$smoothness)
  } else {
    spacetime_cor(new_distance, .sp_temporal_lag(design$temporal_values),
      object$covariance[["scale"]], object$smoothness,
      object$temporal_scale)
  }
  covariance <- object$covariance[["spatial_variance"]] * new_correlation -
    crossprod(cross_covariance, precision_cross)
  if (!is.null(object$random_specification)) {
    covariance <- covariance + .sp_random_new_covariance(
      design$random_z, object$random_specification, object$random_variance
    )
  }
  if (!is.null(object$random_slope_specification)) {
    covariance <- covariance + .sp_random_slope_covariance(
      prediction_random_slope
    )
  }
  covariance <- covariance + crossprod(
    precision_cross,
    object$latent_mode$covariance %*% precision_cross
  )
  if (target == "observation") {
    diag(covariance) <- diag(covariance) +
      object$covariance[["nugget_variance"]]
  }
  list(mean = mean, covariance = covariance)
}

.sp_parameter_specification <- function(object) {
  beta <- object$coefficients
  p <- length(beta)
  covariance_names <- .sp_fitted_covariance_names(object)
  covariance_count <- length(covariance_names)
  if (object$family$family != "gaussian" &&
      !is.null(object$parameter_vcov)) {
    covariance <- object$parameter_vcov
    if (anyNA(covariance) || any(!is.finite(covariance))) {
      stop("Parameter uncertainty is unavailable because the fitted Hessian covariance is unusable.",
           call. = FALSE)
    }
    intervention_count <- if (isTRUE(object$intervention_estimated)) 2L else 0L
    random_count <- if (isTRUE(object$covariance_estimated))
      length(object$random_variance) else 0L
    random_slope_parameterization <- .sp_random_slope_parameterization(
      object$random_slope_specification
    )
    random_slope_count <- length(random_slope_parameterization$values)
    temporal_count <- if (isTRUE(object$temporal_scale_estimated)) 1L else 0L
    values <- c(beta, if (covariance_count)
      log(object$covariance[covariance_names]) else NULL,
      if (random_count) log(object$random_variance) else NULL,
      random_slope_parameterization$values,
      if (temporal_count) log(object$temporal_scale) else NULL,
      if (intervention_count) c(
        stats::qlogis(object$intervention$alpha),
        log(object$intervention$decay)
      ) else NULL)
    if (length(values) != nrow(covariance)) {
      stop("Fitted parameter values and Hessian covariance do not align.",
           call. = FALSE)
    }
    return(list(values = values, covariance = covariance, beta_count = p,
                covariance_count = covariance_count,
                covariance_names = covariance_names,
                random_count = random_count,
                random_slope_count = random_slope_count,
                random_slope_parameterization = random_slope_parameterization,
                temporal_count = temporal_count,
                intervention_count = intervention_count))
  }
  temporal_count <- if (isTRUE(object$temporal_scale_estimated)) 1L else 0L
  random_count <- if (isTRUE(object$covariance_estimated))
    length(object$random_variance) else 0L
  random_slope_parameterization <- .sp_random_slope_parameterization(
    object$random_slope_specification
  )
  random_slope_count <- length(random_slope_parameterization$values)
  values <- beta
  covariance <- object$coefficient_vcov
  if (covariance_count + random_count + temporal_count +
      random_slope_count > 0L) {
    random_names <- if (random_count)
      paste0("random::", names(object$random_variance)) else character()
    variance_names <- c(covariance_names, random_names,
      if (temporal_count) "temporal_scale" else character())
    covariance_values <- c(object$covariance[covariance_names],
      if (random_count) stats::setNames(object$random_variance,
                                        random_names) else numeric(),
      if (temporal_count) c(temporal_scale = object$temporal_scale) else numeric())
    covariance_vcov <- object$variance_parameter_vcov[
      variance_names, variance_names, drop = FALSE
    ]
    if (anyNA(covariance_vcov) || any(!is.finite(covariance_vcov))) {
      stop("Parameter uncertainty is unavailable because the fitted covariance Hessian is unusable.",
           call. = FALSE)
    }
    base_count <- covariance_count + random_count + temporal_count
    total_count <- base_count + random_slope_count
    transformed_vcov <- if (random_slope_count &&
      !is.null(object$transformed_parameter_vcov)) {
      transformed_names <- c(
        paste0("log_", covariance_names),
        if (random_count)
          paste0("log_random::", names(object$random_variance)) else
            character(),
        random_slope_parameterization$names,
        if (temporal_count) "log_temporal_scale" else character()
      )
      object$transformed_parameter_vcov[
        transformed_names, transformed_names, drop = FALSE
      ]
    } else {
      inverse_jacobian <- diag(1 / unname(covariance_values), base_count)
      inverse_jacobian %*% covariance_vcov %*% inverse_jacobian
    }
    if (!identical(dim(transformed_vcov), c(total_count, total_count))) {
      stop("Fitted transformed covariance parameters do not align.",
           call. = FALSE)
    }
    covariance <- rbind(
      cbind(covariance, matrix(0, p, total_count)),
      cbind(matrix(0, total_count, p), transformed_vcov)
    )
    values <- c(values, log(covariance_values),
                random_slope_parameterization$values)
  }
  if (anyNA(covariance) || any(!is.finite(covariance))) {
    stop("Parameter uncertainty is unavailable because its covariance is unusable.",
         call. = FALSE)
  }
  list(values = values, covariance = covariance, beta_count = p,
       covariance_count = covariance_count,
       covariance_names = covariance_names, random_count = random_count,
       random_slope_count = random_slope_count,
       random_slope_parameterization = random_slope_parameterization,
       temporal_count = temporal_count,
       intervention_count = 0L)
}

.sp_object_at_parameters <- function(object, parameters, specification) {
  answer <- object
  p <- specification$beta_count
  answer$coefficients <- stats::setNames(parameters[seq_len(p)],
                                         names(object$coefficients))
  if (specification$covariance_count > 0L) {
    answer$covariance[specification$covariance_names] <- exp(
      parameters[p + seq_len(specification$covariance_count)]
    )
  }
  if (specification$random_count > 0L) {
    random_start <- p + specification$covariance_count
    answer$random_variance <- stats::setNames(exp(parameters[
      random_start + seq_len(specification$random_count)
    ]), names(object$random_variance))
  }
  if (specification$random_slope_count > 0L) {
    slope_start <- p + specification$covariance_count +
      specification$random_count
    answer$random_slope_specification <- .sp_random_slope_update(
      object$random_slope_specification,
      parameters[slope_start + seq_len(specification$random_slope_count)],
      specification$random_slope_parameterization
    )
    answer$random_slope <- .sp_random_slope_export(
      answer$random_slope_specification
    )
  }
  if (specification$temporal_count > 0L) {
    answer$temporal_scale <- exp(parameters[
      p + specification$covariance_count + specification$random_count +
        specification$random_slope_count + 1L
    ])
    answer$temporal$scale <- answer$temporal_scale
  }
  if (specification$intervention_count > 0L) {
    intervention_start <- p + specification$covariance_count +
      specification$random_count + specification$random_slope_count +
      specification$temporal_count
    answer$intervention$alpha <- unname(stats::plogis(
      parameters[intervention_start + 1L]
    ))
    answer$intervention$decay <- unname(exp(
      parameters[intervention_start + 2L]
    ))
    answer$intervention_multiplier <- .sp_intervention_multiplier(
      answer$intervention, object$data_used
    )
  }
  answer$coefficient_vcov[,] <- 0
  if (object$family$family == "gaussian") {
    if (isTRUE(object$spatial)) {
      distance <- .sp_dist(object$coordinates)
      v <- .sp_covariance(
        distance, answer$covariance, object$smoothness, diagonal = TRUE,
        random_z = if (is.null(object$random_specification)) NULL else
          object$random_specification$z,
        random_variance = answer$random_variance,
        random_column_group = if (is.null(object$random_specification)) NULL else
          object$random_specification$column_group,
        random_slope_specification = answer$random_slope_specification,
        temporal_lag = .sp_temporal_lag(object$temporal_values),
        temporal_scale = answer$temporal_scale
      )
      chol_v <- tryCatch(chol(v), error = function(e) NULL)
      if (is.null(chol_v)) return(NULL)
      answer$fitted_state$chol_v <- chol_v
      answer$fitted_state$vinv_x <- .sp_solve_chol(chol_v, object$x)
      answer$residual_mean <- as.vector(
        object$y - object$offset - object$x %*% answer$coefficients
      )
    }
  } else {
    distance <- .sp_dist(object$coordinates)
    state <- .sp_laplace_state(
      answer$coefficients, answer$covariance, distance, object$x,
      object$offset, object$y, object$trials, object$family,
      object$smoothness,
      probability_multiplier = answer$intervention_multiplier,
      temporal_lag = .sp_temporal_lag(object$temporal_values),
      temporal_scale = answer$temporal_scale,
      random_z = if (is.null(object$random_specification)) NULL else
        object$random_specification$z,
      random_variance = answer$random_variance,
      random_column_group = if (is.null(object$random_specification)) NULL else
        object$random_specification$column_group,
      random_slope_specification = answer$random_slope_specification,
      initial_mode = object$latent_mode$mode
    )
    if (is.null(state)) return(NULL)
    answer$fitted_state <- state
    answer$latent_mode <- state$mode
  }
  answer
}

#' Draw jointly from a fitted model's point predictions
#'
#' Preserve dependence among prediction locations while drawing from the
#' exact Gaussian conditional distribution or the local Gaussian predictive
#' approximation for a binomial-logit or Poisson-log fit.
#'
#' @param object A fitted `spmodels_fit` object.
#' @param newdata Optional prediction data containing all predictors and point
#'   coordinates. The fitted locations are used when this is `NULL`.
#' @param nsim Positive number of joint draws.
#' @param target Draw the smooth latent target or a future observation-site
#'   target. The latter includes the independent variance component.
#' @param type Draw on the linear-predictor or response scale. For binomial and
#'   Poisson fits, a latent response draw is a probability or intensity; an
#'   observation response draw is a Bernoulli outcome or Poisson count.
#' @param seed Optional integer seed local to this call.
#' @param parameter_uncertainty Include a local Gaussian approximation to
#'   fitted regression and covariance-parameter uncertainty? This requires a
#'   usable fitted Hessian covariance. Gaussian covariance and regression
#'   blocks are treated as asymptotically independent. Fitted intervention
#'   parameters rebuild the multiplier and latent state for every draw.
#' @param unseen Policy for grouping levels absent during fitting. `"error"`
#'   rejects them; `"marginal"` draws new zero-mean random intercepts with the
#'   fitted variance, shared by rows carrying the same new level.
#'
#' @return An object of class `spmodels_prediction_draws` containing an
#'   `nsim` by number-of-locations draw matrix and its defining metadata.
#'
#' @details Joint draws require a dense covariance over all prediction rows and
#'   become expensive as the grid grows. At 1,000 or more locations and at
#'   least 100 draws, this function warns about the cost. For pointwise means
#'   and standard errors on large grids, use [sp_predict_chunked()], which is
#'   exact but does not retain cross-location dependence. Use joint draws when
#'   spatially dependent uncertainty or draw-based target summaries are needed.
#' @export
sp_predict_draws <- function(object, newdata = NULL, nsim = 1000L,
                             target = c("latent", "observation"),
                             type = c("response", "link"), seed = NULL,
                             parameter_uncertainty = FALSE,
                             unseen = c("error", "marginal")) {
  if (!inherits(object, "spmodels_fit")) {
    stop("`object` must be a fitted `spmodels_fit` model.", call. = FALSE)
  }
  if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) ||
      !is.finite(nsim) || nsim < 1 || nsim != as.integer(nsim)) {
    stop("`nsim` must be one positive integer.", call. = FALSE)
  }
  target <- match.arg(target)
  type <- match.arg(type)
  unseen <- match.arg(unseen)
  if (!is.logical(parameter_uncertainty) ||
      length(parameter_uncertainty) != 1L || is.na(parameter_uncertainty)) {
    stop("`parameter_uncertainty` must be TRUE or FALSE.", call. = FALSE)
  }
  design <- .sp_prediction_design(object, newdata, unseen)
  if (length(design$coordinates) && nrow(design$coordinates) >= 1000L &&
      as.integer(nsim) >= 100L) {
    warning(paste0(
      "Joint prediction draws over at least 1,000 locations can be slow and ",
      "memory-intensive. Use `sp_predict_chunked()` for exact pointwise ",
      "means/SEs; it does not provide joint draws or target summaries."
    ), call. = FALSE)
  }
  distribution <- if (object$family$family == "gaussian") {
    .sp_joint_gaussian(object, design, target)
  } else {
    .sp_joint_nongaussian(object, design, target)
  }
  seed_scope <- .sp_seed_scope(seed)
  if (!is.null(seed_scope)) {
    on.exit(.sp_restore_seed(seed_scope), add = TRUE)
    set.seed(seed_scope$seed)
  }
  parameter_draws <- NULL
  intervention_multiplier_draws <- NULL
  if (parameter_uncertainty) {
    specification <- .sp_parameter_specification(object)
    parameter_root <- .sp_covariance_root(specification$covariance)
    parameter_z <- matrix(stats::rnorm(as.integer(nsim) *
      length(specification$values)), nrow = as.integer(nsim))
    parameter_draws <- if (parameter_root$orientation == "upper") {
      sweep(parameter_z %*% parameter_root$root, 2L,
            specification$values, "+")
    } else {
      sweep(parameter_z %*% t(parameter_root$root), 2L,
            specification$values, "+")
    }
    link_draws <- matrix(NA_real_, nrow = as.integer(nsim),
                         ncol = length(distribution$mean))
    if (specification$intervention_count > 0L) {
      intervention_multiplier_draws <- matrix(
        NA_real_, nrow = as.integer(nsim), ncol = length(distribution$mean)
      )
    }
    for (draw_index in seq_len(as.integer(nsim))) {
      draw_object <- .sp_object_at_parameters(
        object, parameter_draws[draw_index, ], specification
      )
      if (is.null(draw_object)) {
        stop(sprintf("Parameter draw %d produced an unusable fitted state.",
                     draw_index), call. = FALSE)
      }
      draw_distribution <- if (object$family$family == "gaussian") {
        .sp_joint_gaussian(draw_object, design, target)
      } else {
        .sp_joint_nongaussian(draw_object, design, target)
      }
      draw_root <- .sp_covariance_root(draw_distribution$covariance)
      prediction_z <- stats::rnorm(length(draw_distribution$mean))
      link_draws[draw_index, ] <- draw_distribution$mean +
        if (draw_root$orientation == "upper") {
          as.vector(prediction_z %*% draw_root$root)
        } else {
          as.vector(prediction_z %*% t(draw_root$root))
        }
      if (!is.null(intervention_multiplier_draws)) {
        intervention_multiplier_draws[draw_index, ] <-
          .sp_intervention_multiplier(draw_object$intervention, design$data)
      }
    }
    distribution <- list(
      mean = colMeans(link_draws),
      covariance = if (nrow(link_draws) > 1L) stats::cov(link_draws) else
        matrix(NA_real_, ncol(link_draws), ncol(link_draws))
    )
  } else {
    root <- .sp_covariance_root(distribution$covariance)
    z <- matrix(stats::rnorm(as.integer(nsim) * length(distribution$mean)),
                nrow = as.integer(nsim))
    link_draws <- if (root$orientation == "upper") {
      sweep(z %*% root$root, 2L, distribution$mean, "+")
    } else {
      sweep(z %*% t(root$root), 2L, distribution$mean, "+")
    }
  }
  draws <- link_draws
  if (object$family$family == "binomial" && type == "response") {
    probability <- stats::plogis(link_draws)
    multiplier <- if (!is.null(intervention_multiplier_draws))
      intervention_multiplier_draws else design$intervention_multiplier
    if (is.matrix(multiplier)) {
      probability <- probability * multiplier
    } else if (!is.null(multiplier)) {
      probability <- sweep(probability, 2L, multiplier, "*")
    }
    draws <- if (target == "latent") probability else
      matrix(stats::rbinom(length(probability), 1L, probability),
             nrow = nrow(probability))
  } else if (object$family$family == "poisson" && type == "response") {
    intensity <- exp(link_draws)
    draws <- if (target == "latent") intensity else
      matrix(stats::rpois(length(intensity), intensity),
             nrow = nrow(intensity))
  }
  colnames(draws) <- row.names(if (is.null(newdata)) object$model else newdata)
  if (is.null(colnames(draws))) colnames(draws) <- paste0("location_", seq_len(ncol(draws)))
  structure(list(
    draws = draws, predictive_mean = distribution$mean,
    predictive_covariance = distribution$covariance,
    parameter_draws = parameter_draws,
    intervention_multiplier_draws = intervention_multiplier_draws,
    parameter_uncertainty = parameter_uncertainty,
    family = object$family$family, target = target, type = type,
    seed = seed, nsim = as.integer(nsim), call = match.call()
  ), class = "spmodels_prediction_draws")
}

#' @export
print.spmodels_prediction_draws <- function(x, ...) {
  cat("Joint", x$family, x$type, "prediction draws\n")
  cat("Draws:", x$nsim, " Locations:", ncol(x$draws), "\n")
  cat("Target:", x$target, "\n")
  if (!is.null(x$approximation)) {
    cat("Approximation:", x$approximation, "with at most", x$neighbors,
        "predecessors\n")
  }
  invisible(x)
}
