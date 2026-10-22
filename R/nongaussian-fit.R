.sp_response_data <- function(mf, family_name, weights) {
  response <- stats::model.response(mf)
  n <- nrow(mf)
  if (family_name == "binomial") {
    model_weights <- stats::model.weights(mf)
    if (is.matrix(response)) {
      if (ncol(response) != 2L || !is.null(model_weights)) {
        stop("A two-column binomial response cannot also use `weights`.",
             call. = FALSE)
      }
      successes <- response[, 1L]
      failures <- response[, 2L]
      if (anyNA(response) || any(!is.finite(response)) ||
          any(response < 0) || any(response != floor(response))) {
        stop("Binomial successes and failures must be non-negative integer counts.",
             call. = FALSE)
      }
      trials <- successes + failures
    } else {
      if (!is.numeric(response) || anyNA(response) ||
          any(!is.finite(response)) || any(response < 0) ||
          any(response > 1)) {
        stop("A binomial response must be a zero-one value, a proportion, or a two-column count response.",
             call. = FALSE)
      }
      trials <- if (is.null(model_weights)) rep(1, n) else model_weights
      if (!is.numeric(trials) || length(trials) != n || anyNA(trials) ||
          any(!is.finite(trials)) || any(trials < 0) ||
          any(trials != floor(trials))) {
        stop("Binomial `weights` must be non-negative integer trial counts.",
             call. = FALSE)
      }
      successes <- response * trials
      if (any(abs(successes - round(successes)) > 1e-8)) {
        stop("Binomial proportions multiplied by `weights` must give integer successes.",
             call. = FALSE)
      }
      successes <- round(successes)
    }
    encoding <- if (is.matrix(response)) "counts" else if (is.null(model_weights))
      "bernoulli" else "proportion"
    return(list(y = as.numeric(successes), trials = as.numeric(trials),
                response = response, encoding = encoding))
  }
  if (!is.null(stats::model.weights(mf))) {
    stop("Poisson `weights` are not implemented; use an offset for exposure.",
         call. = FALSE)
  }
  if (!is.numeric(response) || is.matrix(response) || anyNA(response) ||
      any(!is.finite(response)) || any(response < 0) ||
      any(response != floor(response))) {
    stop("The Poisson response must contain non-negative integer counts.",
         call. = FALSE)
  }
  list(y = as.numeric(response), trials = NULL, response = response,
       encoding = "count")
}

.sp_laplace_state <- function(beta, covariance, distance, x, offset, y,
                               trials, family, smoothness,
                               probability_multiplier = NULL,
                               temporal_lag = NULL,
                               temporal_scale = NULL,
                               random_z = NULL,
                               random_variance = numeric(),
                               random_column_group = NULL,
                               random_slope_specification = NULL,
                               initial_mode = NULL) {
  sigma <- .sp_covariance(distance, covariance, smoothness, diagonal = TRUE,
                          temporal_lag = temporal_lag,
                          temporal_scale = temporal_scale,
                          random_z = random_z,
                          random_variance = random_variance,
                          random_column_group = random_column_group)
  if (!is.null(random_slope_specification)) {
    sigma <- sigma + .sp_random_slope_covariance(
      random_slope_specification
    )
  }
  chol_sigma <- tryCatch(chol(sigma), error = function(e) NULL)
  if (is.null(chol_sigma)) return(NULL)
  mode <- tryCatch(
    sp_latent_mode(y, as.vector(x %*% beta + offset), sigma,
                   family = family, trials = trials,
                   probability_multiplier = probability_multiplier,
                   initial = initial_mode),
    error = function(e) NULL
  )
  if (is.null(mode) || !mode$converged) return(NULL)
  chol_mode_covariance <- tryCatch(chol(mode$covariance),
                                   error = function(e) NULL)
  if (is.null(chol_mode_covariance)) return(NULL)
  log_determinant_sigma <- 2 * sum(log(diag(chol_sigma)))
  log_determinant_mode_covariance <-
    2 * sum(log(diag(chol_mode_covariance)))
  log_likelihood <- mode$log_posterior - 0.5 * log_determinant_sigma +
    0.5 * log_determinant_mode_covariance
  if (!is.finite(log_likelihood)) return(NULL)
  list(beta = beta, covariance = covariance, sigma = sigma,
       chol_sigma = chol_sigma, mode = mode,
       log_likelihood = log_likelihood)
}

.sp_fit_nongaussian <- function(call, formula, data, family, coords, spatial,
                                 weights, subset, na.action, offset, control,
                                 verbose, intervention = NULL, time = NULL,
                                 random = NULL, random_slope = list()) {
  if (!spatial) {
    stop("Independent non-Gaussian fitting is not yet implemented.",
         call. = FALSE)
  }
  model_data <- if (inherits(data, "sf")) {
    if (!requireNamespace("sf", quietly = TRUE)) {
      stop("Package `sf` is required for `sf` input.", call. = FALSE)
    }
    sf::st_drop_geometry(data)
  } else data
  mf_call <- call[c(1L, match(c("formula", "data", "subset", "na.action",
                               "offset", "weights"), names(call), 0L))]
  mf_call[[1L]] <- quote(stats::model.frame)
  mf_call[["formula"]] <- quote(formula)
  mf_call[["data"]] <- quote(model_data)
  if (!is.null(weights)) {
    # Insert the already data-masked vector itself.  model.frame evaluates
    # model arguments in the formula environment, where the internal name
    # `weights` could otherwise resolve to stats::weights().
    mf_call[["weights"]] <- weights
  } else if ("weights" %in% names(mf_call)) {
    mf_call[["weights"]] <- NULL
  }
  mf <- eval(mf_call, environment())
  terms_object <- attr(mf, "terms")
  x <- stats::model.matrix(terms_object, mf)
  if (qr(x)$rank < ncol(x)) {
    stop("The mean-model design matrix is rank deficient.", call. = FALSE)
  }
  response <- .sp_response_data(mf, family$family, weights)
  y <- response$y
  trials <- response$trials
  model_offset <- stats::model.offset(mf)
  if (is.null(model_offset)) model_offset <- rep(0, length(y))
  row_index <- match(row.names(mf), row.names(model_data))
  if (anyNA(row_index)) stop("Could not align model rows with `data`.", call. = FALSE)
  if (!is.null(intervention) && family$family != "binomial") {
    stop("Decay-adjusted intervention histories require a binomial-logit model.",
         call. = FALSE)
  }
  probability_multiplier <- .sp_intervention_multiplier(
    intervention, data, row_index
  )
  temporal_values <- .sp_temporal_values(time, data, row_index)
  temporal_lag <- .sp_temporal_lag(temporal_values)
  if (!is.null(time) && isTRUE(time$estimate) &&
      !any(temporal_lag[upper.tri(temporal_lag)] > 0)) {
    stop("Estimated temporal scale requires at least two distinct observed times.",
         call. = FALSE)
  }
  random_specification <- .sp_random_training(
    random, data[row_index, , drop = FALSE]
  )
  random_slope_specification <- .sp_random_slope_training(
    random_slope, data[row_index, , drop = FALSE]
  )
  .sp_random_slope_identifiability(random_slope_specification,
                                   random_specification)
  coordinate_info <- .sp_coordinates(coords, data, row_index)
  coordinates <- coordinate_info$values
  distance_units <- control$distance_units
  if (is.null(distance_units)) distance_units <- coordinate_info$inferred_units
  if (is.null(distance_units)) {
    stop("Supply explicit coordinate units with `sp_control(distance_units = ...)`.",
         call. = FALSE)
  }
  coordinate_multiplier <- .sp_distance_unit_multiplier(
    coordinate_info$inferred_units, distance_units
  )
  coordinates <- coordinates * coordinate_multiplier
  distance <- .sp_dist(coordinates)
  positive_distances <- distance[upper.tri(distance) & distance > 0]
  if (!length(positive_distances)) {
    stop("At least two distinct coordinate locations are required.", call. = FALSE)
  }
  if (!control$estimate && any(vapply(
    list(control$spatial_variance, control$scale, control$nugget_variance),
    is.null, logical(1)))) {
    stop("Fixed spatial fitting requires `spatial_variance`, `scale`, and `nugget_variance`.",
         call. = FALSE)
  }

  glm_y <- if (family$family == "binomial") y / pmax(trials, 1) else y
  glm_weights <- if (family$family == "binomial") trials else rep(1, length(y))
  initial_glm <- suppressWarnings(stats::glm.fit(
    x, glm_y, weights = glm_weights, offset = model_offset, family = family
  ))
  initial_beta <- initial_glm$coefficients
  initial_beta[!is.finite(initial_beta)] <- 0
  initial_covariance <- c(
    spatial_variance = if (is.null(control$spatial_variance)) 0.8 else
      control$spatial_variance,
    scale = if (is.null(control$scale)) stats::median(positive_distances) else
      control$scale,
    nugget_variance = if (is.null(control$nugget_variance)) 0.2 else
      control$nugget_variance
  )
  covariance_estimated_names <- .sp_estimated_covariance_names(control)
  random_initial <- .sp_resolve_random_variance(
    control$random_variance,
    if (is.null(random_specification)) character() else
      random_specification$variables,
    control$estimate, default = 0.1
  )
  p <- ncol(x)
  covariance_count <- length(covariance_estimated_names)
  random_count <- if (control$estimate) length(random_initial) else 0L
  random_slope_parameterization <- .sp_random_slope_parameterization(
    random_slope_specification
  )
  random_slope_count <- length(random_slope_parameterization$values)
  temporal_estimated <- !is.null(time) && isTRUE(time$estimate)
  temporal_count <- if (temporal_estimated) 1L else 0L
  intervention_estimated <- !is.null(intervention) &&
    isTRUE(intervention$estimate)
  intervention_count <- if (intervention_estimated) 2L else 0L
  parameter_start <- c(
    initial_beta,
    if (covariance_count)
      log(initial_covariance[covariance_estimated_names]) else NULL,
    if (random_count) log(random_initial) else NULL,
    random_slope_parameterization$values,
    if (temporal_estimated) log(time$scale) else NULL,
    if (intervention_estimated)
      c(stats::qlogis(intervention$alpha), log(intervention$decay)) else NULL
  )
  names(parameter_start) <- c(
    colnames(x),
    if (covariance_count) paste0("log_", covariance_estimated_names) else NULL,
    if (random_count) paste0("log_random::", names(random_initial)) else NULL,
    random_slope_parameterization$names,
    if (temporal_estimated) "log_temporal_scale" else NULL,
    if (intervention_estimated) c("logit_alpha", "log_decay") else NULL
  )
  intervention_indices <- if (intervention_estimated) {
    p + covariance_count + random_count + random_slope_count +
      temporal_count + seq_len(2L)
  } else integer()
  resolve_intervention <- function(parameters) {
    if (!intervention_estimated) return(intervention)
    current <- intervention
    current$alpha <- unname(stats::plogis(
      parameters[intervention_indices[[1L]]]
    ))
    current$decay <- unname(exp(parameters[intervention_indices[[2L]]]))
    current
  }
  objective <- function(parameters) {
    value <- tryCatch({
      beta <- parameters[seq_len(p)]
      covariance <- initial_covariance
      if (covariance_count) {
        covariance[covariance_estimated_names] <-
          exp(parameters[p + seq_len(covariance_count)])
      }
      random_variance <- if (random_count) {
        answer <- exp(parameters[p + covariance_count + seq_len(random_count)])
        stats::setNames(answer, names(random_initial))
      } else random_initial
      random_slope_start <- p + covariance_count + random_count
      current_random_slope <- .sp_random_slope_update(
        random_slope_specification,
        parameters[random_slope_start + seq_len(random_slope_count)],
        random_slope_parameterization
      )
      current_temporal_scale <- if (temporal_estimated) {
        exp(parameters[random_slope_start + random_slope_count + 1L])
      } else if (!is.null(time)) time$scale else NULL
      current_intervention <- resolve_intervention(parameters)
      current_multiplier <- .sp_intervention_multiplier(
        current_intervention, data, row_index
      )
      state <- .sp_laplace_state(beta, covariance, distance, x, model_offset,
                                 y, trials, family, control$smoothness,
                                 current_multiplier, temporal_lag,
                                 current_temporal_scale,
                                 random_z = if (is.null(random_specification))
                                   NULL else random_specification$z,
                                 random_variance = random_variance,
                                 random_column_group = if (is.null(random_specification))
                                   NULL else random_specification$column_group,
                                 random_slope_specification = current_random_slope)
      if (is.null(state)) NULL else {
        penalty <- if (intervention_estimated &&
          !is.null(current_intervention$alpha_penalty)) {
          stats::dbeta(
            current_intervention$alpha,
            current_intervention$alpha_penalty[["shape1"]],
            current_intervention$alpha_penalty[["shape2"]], log = TRUE
          )
        } else 0
        -(state$log_likelihood + penalty)
      }
    }, error = function(e) NULL)
    if (is.null(value) || !is.finite(value)) 1e20 else value
  }
  lower <- c(
    rep(-Inf, p),
    if (covariance_count) log(c(
      spatial_variance = 1e-7,
      scale = min(positive_distances) * 1e-4,
      nugget_variance = 1e-7
    )[covariance_estimated_names]) else NULL,
    if (random_count) rep(log(1e-7), random_count) else NULL,
    unlist(lapply(random_slope_parameterization$components, function(details) {
      index <- details$index
      vapply(seq_len(nrow(index)), function(k) {
        if (index[k, 1L] == index[k, 2L]) log(sqrt(1e-7)) else -sqrt(1e5)
      }, numeric(1L))
    }), use.names = FALSE),
    if (temporal_estimated) {
      positive_time <- temporal_lag[upper.tri(temporal_lag) &
        temporal_lag > 0]
      log(if (length(positive_time)) min(positive_time) * 1e-4 else
        time$scale * 1e-4)
    } else NULL
  )
  upper <- c(
    rep(Inf, p),
    if (covariance_count) log(c(
      spatial_variance = 1e5,
      scale = max(positive_distances) * 1e4,
      nugget_variance = 1e5
    )[covariance_estimated_names]) else NULL,
    if (random_count) rep(log(1e5), random_count) else NULL,
    unlist(lapply(random_slope_parameterization$components, function(details) {
      index <- details$index
      vapply(seq_len(nrow(index)), function(k) {
        if (index[k, 1L] == index[k, 2L]) log(sqrt(1e5)) else sqrt(1e5)
      }, numeric(1L))
    }), use.names = FALSE),
    if (temporal_estimated) {
      positive_time <- temporal_lag[upper.tri(temporal_lag) &
        temporal_lag > 0]
      log(if (length(positive_time)) max(positive_time) * 1e4 else
        time$scale * 1e4)
    } else NULL
  )
  if (intervention_estimated) {
    survey_times <- model_data[[intervention$time]][row_index]
    time_span <- diff(range(c(survey_times, intervention$history$time)))
    if (!is.finite(time_span) || time_span <= 0) time_span <- 1
    lower <- c(lower, stats::qlogis(1e-6), log(time_span * 1e-6))
    upper <- c(upper, stats::qlogis(1 - 1e-6), log(time_span * 1e4))
  }
  optimisation <- stats::optim(
    parameter_start, objective, method = "L-BFGS-B", lower = lower,
    upper = upper, control = list(maxit = control$maxit,
      factr = control$reltol / .Machine$double.eps,
      trace = if (verbose) 1 else 0), hessian = TRUE
  )
  beta <- optimisation$par[seq_len(p)]
  names(beta) <- colnames(x)
  covariance <- initial_covariance
  if (covariance_count) {
    covariance[covariance_estimated_names] <-
      exp(optimisation$par[p + seq_len(covariance_count)])
  }
  random_variance <- if (random_count) {
    answer <- exp(optimisation$par[p + covariance_count + seq_len(random_count)])
    stats::setNames(answer, names(random_initial))
  } else random_initial
  random_slope_start <- p + covariance_count + random_count
  random_slope_specification <- .sp_random_slope_update(
    random_slope_specification,
    optimisation$par[random_slope_start + seq_len(random_slope_count)],
    random_slope_parameterization
  )
  temporal_scale <- if (temporal_estimated) {
    unname(exp(optimisation$par[
      random_slope_start + random_slope_count + 1L
    ]))
  } else if (!is.null(time)) time$scale else NULL
  fitted_intervention <- resolve_intervention(optimisation$par)
  probability_multiplier <- .sp_intervention_multiplier(
    fitted_intervention, data, row_index
  )
  state <- .sp_laplace_state(beta, covariance, distance, x, model_offset,
                             y, trials, family, control$smoothness,
                             probability_multiplier, temporal_lag,
                             temporal_scale,
                             random_z = if (is.null(random_specification))
                               NULL else random_specification$z,
                             random_variance = random_variance,
                             random_column_group = if (is.null(random_specification))
                               NULL else random_specification$column_group,
                             random_slope_specification = random_slope_specification)
  if (is.null(state)) stop("The final Laplace approximation could not be evaluated.",
                           call. = FALSE)
  if (optimisation$convergence != 0L) {
    warning(sprintf("Laplace optimisation ended with convergence code %d: %s",
                    optimisation$convergence, optimisation$message),
            call. = FALSE)
  }
  parameter_vcov <- matrix(NA_real_, length(optimisation$par),
                           length(optimisation$par))
  if (all(is.finite(optimisation$hessian))) {
    chol_hessian <- tryCatch(chol(optimisation$hessian),
                             error = function(e) NULL)
    if (!is.null(chol_hessian)) parameter_vcov <- chol2inv(chol_hessian)
  }
  dimnames(parameter_vcov) <- list(names(parameter_start), names(parameter_start))
  coefficient_vcov <- parameter_vcov[seq_len(p), seq_len(p), drop = FALSE]
  dimnames(coefficient_vcov) <- list(colnames(x), colnames(x))
  covariance_vcov <- matrix(NA_real_, 3L, 3L,
    dimnames = list(names(covariance), names(covariance)))
  if (covariance_count && all(is.finite(parameter_vcov))) {
    log_indices <- p + seq_len(covariance_count)
    jacobian <- diag(unname(covariance[covariance_estimated_names]),
                     covariance_count)
    covariance_vcov[covariance_estimated_names,
      covariance_estimated_names] <- jacobian %*%
      parameter_vcov[log_indices, log_indices, drop = FALSE] %*% jacobian
  }
  random_variance_vcov <- matrix(
    NA_real_, length(random_variance), length(random_variance),
    dimnames = list(names(random_variance), names(random_variance))
  )
  if (random_count && all(is.finite(parameter_vcov))) {
    random_indices <- p + covariance_count + seq_len(random_count)
    jacobian <- diag(unname(random_variance), random_count)
    random_variance_vcov <- jacobian %*%
      parameter_vcov[random_indices, random_indices, drop = FALSE] %*%
      jacobian
    dimnames(random_variance_vcov) <- list(names(random_variance),
                                           names(random_variance))
  }
  random_slope_values <- .sp_random_slope_parameters(
    random_slope_specification
  )
  random_slope_vcov <- matrix(
    NA_real_, length(random_slope_values), length(random_slope_values),
    dimnames = list(names(random_slope_values), names(random_slope_values))
  )
  if (random_slope_count && all(is.finite(parameter_vcov))) {
    slope_indices <- random_slope_start + seq_len(random_slope_count)
    slope_jacobian <- .sp_random_slope_natural_jacobian(
      random_slope_specification
    )
    estimated_rows <- rowSums(abs(slope_jacobian)) > 0
    random_slope_vcov[estimated_rows, estimated_rows] <-
      slope_jacobian[estimated_rows, , drop = FALSE] %*%
      parameter_vcov[slope_indices, slope_indices, drop = FALSE] %*%
      t(slope_jacobian[estimated_rows, , drop = FALSE])
  }
  temporal_scale_se <- NA_real_
  if (temporal_estimated && all(is.finite(parameter_vcov))) {
    temporal_index <- random_slope_start + random_slope_count + 1L
    temporal_variance <- parameter_vcov[temporal_index, temporal_index]
    if (temporal_variance >= 0) {
      temporal_scale_se <- temporal_scale * sqrt(temporal_variance)
    }
  }
  intervention_vcov <- matrix(
    NA_real_, 2L, 2L,
    dimnames = list(c("alpha", "decay"), c("alpha", "decay"))
  )
  if (intervention_estimated && all(is.finite(parameter_vcov))) {
    jacobian <- diag(c(
      fitted_intervention$alpha * (1 - fitted_intervention$alpha),
      fitted_intervention$decay
    ), 2L)
    intervention_vcov <- jacobian %*%
      parameter_vcov[intervention_indices, intervention_indices,
                     drop = FALSE] %*% jacobian
    dimnames(intervention_vcov) <- list(c("alpha", "decay"),
                                        c("alpha", "decay"))
  }
  penalty_value <- if (intervention_estimated &&
    !is.null(fitted_intervention$alpha_penalty)) {
    stats::dbeta(fitted_intervention$alpha,
      fitted_intervention$alpha_penalty[["shape1"]],
      fitted_intervention$alpha_penalty[["shape2"]], log = TRUE)
  } else 0
  structure(list(
    call = call, formula = formula, terms = terms_object, spatial = TRUE,
    contrasts = attr(x, "contrasts"),
    xlevels = stats::.getXlevels(terms_object, mf), coords = coords,
    coordinate_names = colnames(coordinates),
    coordinate_source = coordinate_info$source, crs = coordinate_info$crs,
    distance_units = distance_units,
    coordinate_multiplier = coordinate_multiplier,
    family = family, model = mf, x = x,
    y = y, trials = trials, offset = model_offset, coordinates = coordinates,
    data_used = data[row_index, , drop = FALSE], control = control,
    response_encoding = response$encoding,
    intervention = fitted_intervention,
    intervention_multiplier = probability_multiplier,
    random = random, random_specification = random_specification,
    random_variance = random_variance,
    random_slope = .sp_random_slope_export(random_slope_specification),
    random_slope_specification = random_slope_specification,
    random_variance_vcov = random_variance_vcov,
    random_slope_vcov = random_slope_vcov,
    random_slope_parameterization = random_slope_parameterization,
    temporal = time, temporal_values = temporal_values,
    temporal_scale = temporal_scale,
    temporal_scale_estimated = temporal_estimated,
    temporal_scale_se = temporal_scale_se,
    intervention_estimated = intervention_estimated,
    intervention_vcov = intervention_vcov,
    external_offset = "offset" %in% names(call),
    covariance = covariance, smoothness = control$smoothness,
    covariance_estimated = covariance_count > 0L,
    covariance_estimated_names = covariance_estimated_names,
    covariance_vcov = covariance_vcov, coefficients = beta,
    coefficient_vcov = coefficient_vcov,
    parameter_vcov = parameter_vcov,
    residual_mean = y - state$mode$conditional_mean,
    fitted_state = state, latent_mode = state$mode,
    log_likelihood = state$log_likelihood + penalty_value,
    likelihood_method = if (intervention_estimated &&
      !is.null(fitted_intervention$alpha_penalty)) {
      "Penalized Laplace"
    } else "Laplace",
    optimisation = optimisation, nobs = length(y), rank = p
  ), class = "spmodels_fit")
}
