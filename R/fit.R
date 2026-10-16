.sp_dist <- function(a, b = a) {
  aa <- rowSums(a * a)
  bb <- rowSums(b * b)
  squared <- outer(aa, bb, "+") - 2 * tcrossprod(a, b)
  sqrt(pmax(squared, 0))
}

.sp_solve_chol <- function(chol_v, rhs) {
  backsolve(chol_v, forwardsolve(t(chol_v), rhs))
}

.sp_coordinates <- function(coords, data, rows = NULL, expected_crs = NULL) {
  if (inherits(data, "sf")) {
    if (!requireNamespace("sf", quietly = TRUE)) {
      stop("Package `sf` is required for `sf` input.", call. = FALSE)
    }
    if (!is.null(coords)) {
      stop("Do not supply `coords` when `data` is an `sf` object.", call. = FALSE)
    }
    geometry_types <- unique(as.character(sf::st_geometry_type(data)))
    if (!length(geometry_types) || any(!geometry_types %in% "POINT")) {
      stop("`sf` input must contain only non-empty POINT geometries.", call. = FALSE)
    }
    if (any(sf::st_is_empty(data))) {
      stop("`sf` input must contain only non-empty POINT geometries.", call. = FALSE)
    }
    crs <- sf::st_crs(data)
    if (is.na(crs)) {
      stop("`sf` input must have a known projected coordinate reference system.",
           call. = FALSE)
    }
    if (isTRUE(sf::st_is_longlat(data))) {
      stop("Longitude/latitude coordinates must be projected before fitting.",
           call. = FALSE)
    }
    if (!is.null(expected_crs) && !isTRUE(crs == expected_crs)) {
      stop("`newdata` must use the same coordinate reference system as the fitted data.",
           call. = FALSE)
    }
    answer <- sf::st_coordinates(data)
    if (ncol(answer) != 2L) {
      stop("Only two-dimensional POINT geometry is currently supported.",
           call. = FALSE)
    }
    if (!is.null(rows)) answer <- answer[rows, , drop = FALSE]
    units <- crs$units_gdal
    if (is.null(units) || is.na(units) || !nzchar(units)) units <- NULL
    return(list(values = answer, crs = crs, inferred_units = units,
                source = "sf"))
  }
  if (!inherits(coords, "formula") || length(coords) != 2L) {
    stop("`coords` must be a one-sided formula containing two variables, such as `~ x + y`.",
         call. = FALSE)
  }
  coordinate_frame <- stats::model.frame(coords, data = data,
                                          na.action = stats::na.pass)
  if (ncol(coordinate_frame) != 2L ||
      !all(vapply(coordinate_frame, is.numeric, logical(1)))) {
    stop("`coords` must select exactly two numeric coordinate variables.",
         call. = FALSE)
  }
  answer <- as.matrix(coordinate_frame)
  if (!is.null(rows)) answer <- answer[rows, , drop = FALSE]
  if (anyNA(answer) || any(!is.finite(answer))) {
    stop("Coordinate values must be finite and non-missing.", call. = FALSE)
  }
  list(values = answer, crs = NULL, inferred_units = NULL,
       source = "columns")
}

.sp_covariance <- function(distance, pars, smoothness, diagonal = FALSE,
                           random_z = NULL, random_variance = numeric(),
                           random_column_group = NULL,
                           random_slope_specification = NULL,
                           temporal_lag = NULL, temporal_scale = NULL) {
  correlation <- if (is.null(temporal_lag)) {
    .sp_matern_cor_unchecked(distance, pars[["scale"]], smoothness)
  } else {
    .sp_matern_cor_unchecked(distance, pars[["scale"]], smoothness) *
      exp(-temporal_lag / temporal_scale)
  }
  answer <- pars[["spatial_variance"]] * correlation
  if (!is.null(random_z) && length(random_variance)) {
    answer <- answer + .sp_random_covariance(
      random_z, random_z, random_column_group, random_variance
    )
  }
  if (!is.null(random_slope_specification)) {
    answer <- answer + .sp_random_slope_covariance(
      random_slope_specification
    )
  }
  if (diagonal) diag(answer) <- diag(answer) + pars[["nugget_variance"]]
  answer
}

.sp_gaussian_state <- function(pars, distance, x, y_centered, smoothness,
                               random_z = NULL,
                               random_variance = numeric(),
                               random_column_group = NULL,
                               random_slope_specification = NULL,
                               temporal_lag = NULL, temporal_scale = NULL) {
  v <- .sp_covariance(distance, pars, smoothness, diagonal = TRUE,
                      random_z = random_z,
                      random_variance = random_variance,
                      random_column_group = random_column_group,
                      random_slope_specification = random_slope_specification,
                      temporal_lag = temporal_lag,
                      temporal_scale = temporal_scale)
  .sp_gaussian_state_v(v, x, y_centered)
}

.sp_gaussian_state_v <- function(v, x, y_centered) {
  chol_v <- tryCatch(chol(v), error = function(e) NULL)
  if (is.null(chol_v)) return(NULL)
  vinv_x <- .sp_solve_chol(chol_v, x)
  xt_vinv_x <- crossprod(x, vinv_x)
  chol_information <- tryCatch(chol(xt_vinv_x), error = function(e) NULL)
  if (is.null(chol_information)) return(NULL)
  vinv_y <- .sp_solve_chol(chol_v, y_centered)
  beta <- .sp_solve_chol(chol_information, crossprod(x, vinv_y))
  residual <- as.vector(y_centered - x %*% beta)
  quadratic <- sum(residual * .sp_solve_chol(chol_v, residual))
  n <- length(y_centered)
  log_likelihood <- -0.5 * (n * log(2 * pi) +
    2 * sum(log(diag(chol_v))) + quadratic)
  list(v = v, chol_v = chol_v, vinv_x = vinv_x,
       information = xt_vinv_x, chol_information = chol_information,
       beta = as.vector(beta), residual = residual,
       log_likelihood = log_likelihood)
}

#' Fit a geostatistical regression
#'
#' Fit Gaussian identity-link regression by maximum likelihood, or
#' binomial-logit and Poisson-log latent Gaussian regression by a Laplace
#' approximation. Models use a linear mean and stationary isotropic Matérn
#' covariance with fixed smoothness.
#'
#' @param formula,data A model formula and data frame.
#' @param family A Gaussian identity, binomial logit, or Poisson log family.
#' @param coords A one-sided formula selecting exactly two numeric coordinate
#'   columns.
#' @param spatial Fit spatial dependence? `FALSE` fits the independent Gaussian
#'   maximum-likelihood limiting model without coordinates.
#' @param random Optional one-sided additive formula of grouping variables for
#'   independent latent Gaussian random intercepts in any supported spatial
#'   family, such as `~ clinic + survey`.
#' @param random_slope Optional grouped random-slope specification made by
#'   [sp_random_slope()], or a list of additive specifications. Each covariance
#'   may be fixed or estimated through its constrained Cholesky parameters.
#' @param time Optional separable temporal specification created by
#'   [sp_temporal()].
#' @param intervention Optional fixed decay-adjusted intervention history
#'   created by [sp_intervention()] for a binomial-logit fit.
#' @param weights Integer binomial trial counts for a proportion response.
#' @param subset Optional row subset, following `model.frame` conventions.
#' @param na.action Missing-value action. The default fails explicitly.
#' @param offset Optional known linear-predictor offset.
#' @param control An object returned by [sp_control()].
#' @param verbose Logical scalar controlling optimisation tracing.
#'
#' @return An object of class `spmodels_fit`.
#' @export
sp_fit <- function(formula, data, family = stats::gaussian(), coords = NULL,
                   spatial = TRUE, time = NULL, random = NULL,
                   intervention = NULL, weights = NULL, subset = NULL,
                   na.action = stats::na.fail, offset = NULL,
                   control = sp_control(), verbose = FALSE,
                   random_slope = NULL) {
  if (!is.data.frame(data)) {
    stop("This implementation requires `data` to be a data frame.", call. = FALSE)
  }
  if (!inherits(formula, "formula")) stop("`formula` must be a formula.", call. = FALSE)
  if (is.function(family)) family <- family()
  if (!inherits(family, "family") ||
      !((family$family == "gaussian" && family$link == "identity") ||
        (family$family == "binomial" && family$link == "logit") ||
        (family$family == "poisson" && family$link == "log"))) {
    stop("Supported families are Gaussian identity, binomial logit, and Poisson log.",
         call. = FALSE)
  }
  if (!is.logical(spatial) || length(spatial) != 1L || is.na(spatial)) {
    stop("`spatial` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(time) && !inherits(time, "spmodels_temporal")) {
    stop("`time` must be created by `sp_temporal()`.", call. = FALSE)
  }
  if (!is.null(time) && !spatial) {
    stop("A temporal covariance currently requires `spatial = TRUE`.",
         call. = FALSE)
  }
  if (!inherits(control, "spmodels_control")) {
    stop("`control` must be created by `sp_control()`.", call. = FALSE)
  }
  if (!is.null(time) && isTRUE(time$estimate) && !isTRUE(control$estimate)) {
    stop("Estimated temporal scale requires `sp_control(estimate = TRUE)`.",
         call. = FALSE)
  }
  if (!is.null(intervention) &&
      !inherits(intervention, "spmodels_intervention")) {
    stop("`intervention` must be created by `sp_intervention()`.",
         call. = FALSE)
  }
  if (is.null(random) && !is.null(control$random_variance)) {
    stop("`random_variance` was supplied but no `random` grouping formula was declared.",
         call. = FALSE)
  }
  random_slope <- .sp_random_slope_list(random_slope)
  if (!isTRUE(control$estimate) && any(vapply(random_slope, function(component) {
    isTRUE(component$estimate)
  }, logical(1)))) {
    stop("Estimated random-slope covariance requires `sp_control(estimate = TRUE)`.",
         call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }

  # Resolve a supplied weights expression with data-mask semantics before the
  # internal fitting call is reconstructed.  This supports both column names
  # (the usual model-frame interface) and programmatically supplied vectors.
  if (!missing(weights)) {
    weights <- eval(substitute(weights), data, parent.frame())
  }

  call <- match.call()
  if (family$family != "gaussian") {
    return(.sp_fit_nongaussian(
      call, formula, data, family, coords, spatial, weights, subset,
      na.action, offset, control, verbose, intervention, time, random,
      random_slope
    ))
  }
  if (!is.null(intervention)) {
    stop("Decay-adjusted intervention histories require a binomial-logit model.",
         call. = FALSE)
  }
  if (!is.null(weights)) {
    stop("Gaussian `weights` are not implemented.", call. = FALSE)
  }
  if ((!is.null(random) || length(random_slope)) && !isTRUE(spatial)) {
    stop("Grouped random effects currently require `spatial = TRUE`.",
         call. = FALSE)
  }
  model_data <- if (inherits(data, "sf")) {
    if (!requireNamespace("sf", quietly = TRUE)) {
      stop("Package `sf` is required for `sf` input.", call. = FALSE)
    }
    sf::st_drop_geometry(data)
  } else data
  mf_call <- call[c(1L, match(c("formula", "data", "subset", "na.action", "offset"),
                              names(call), 0L))]
  mf_call[[1L]] <- quote(stats::model.frame)
  mf_call[["formula"]] <- quote(formula)
  mf_call[["data"]] <- quote(model_data)
  mf <- eval(mf_call, environment())
  y <- stats::model.response(mf, "numeric")
  if (!is.numeric(y) || is.matrix(y)) {
    stop("The Gaussian response must be one numeric vector.", call. = FALSE)
  }
  terms_object <- attr(mf, "terms")
  x <- stats::model.matrix(terms_object, mf)
  if (qr(x)$rank < ncol(x)) {
    stop("The mean-model design matrix is rank deficient.", call. = FALSE)
  }
  row_index <- match(row.names(mf), row.names(model_data))
  if (anyNA(row_index)) stop("Could not align model rows with `data`.", call. = FALSE)
  model_offset <- stats::model.offset(mf)
  if (is.null(model_offset)) model_offset <- rep(0, length(y))
  y_centered <- y - model_offset

  if (!spatial) {
    ordinary <- stats::lm.fit(x, y_centered)
    residual <- as.vector(ordinary$residuals)
    residual_variance <- if (control$estimate) {
      mean(residual^2)
    } else {
      if (is.null(control$nugget_variance)) {
        stop("Fixed non-spatial fitting requires `nugget_variance`.", call. = FALSE)
      }
      control$nugget_variance
    }
    if (!is.finite(residual_variance) || residual_variance <= 0) {
      stop("The independent Gaussian residual variance is not positive.",
           call. = FALSE)
    }
    state <- .sp_gaussian_state_v(diag(residual_variance, length(y)),
                                  x, y_centered)
    names(state$beta) <- colnames(x)
    beta_covariance <- chol2inv(state$chol_information)
    dimnames(beta_covariance) <- list(colnames(x), colnames(x))
    covariance <- c(spatial_variance = 0, scale = NA_real_,
                    nugget_variance = residual_variance)
    covariance_vcov <- matrix(NA_real_, 3L, 3L,
                              dimnames = list(names(covariance), names(covariance)))
    if (control$estimate) {
      covariance_vcov["nugget_variance", "nugget_variance"] <-
        2 * residual_variance^2 / length(y)
    }
    return(structure(list(
      call = call, formula = formula, terms = terms_object,
      contrasts = attr(x, "contrasts"),
      xlevels = stats::.getXlevels(terms_object, mf),
      spatial = FALSE, coords = NULL, coordinate_names = NULL,
      coordinate_source = "none", crs = NULL, distance_units = NULL,
      family = family, model = mf, x = x, y = y, offset = model_offset,
      data_used = data[row_index, , drop = FALSE], control = control,
      external_offset = "offset" %in% names(call),
      coordinates = NULL, covariance = covariance, smoothness = NA_real_,
      covariance_estimated = control$estimate,
      covariance_vcov = covariance_vcov,
      coefficients = state$beta, coefficient_vcov = beta_covariance,
      residual_mean = state$residual, fitted_state = state,
      log_likelihood = state$log_likelihood, optimisation = NULL,
      nobs = length(y), rank = ncol(x)
    ), class = "spmodels_fit"))
  }

  coordinate_info <- .sp_coordinates(coords, data, row_index)
  coordinates <- coordinate_info$values
  distance_units <- control$distance_units
  if (is.null(distance_units)) distance_units <- coordinate_info$inferred_units
  if (is.null(distance_units)) {
    stop("Supply explicit coordinate units with `sp_control(distance_units = ...)`.",
         call. = FALSE)
  }
  if (nrow(coordinates) != length(y)) stop("Coordinate and model rows do not align.", call. = FALSE)
  random_specification <- .sp_random_training(
    random, data[row_index, , drop = FALSE]
  )
  random_slope_specification <- .sp_random_slope_training(
    random_slope, data[row_index, , drop = FALSE]
  )
  .sp_random_slope_identifiability(random_slope_specification,
                                   random_specification)
  temporal_values <- .sp_temporal_values(time, data, row_index)
  temporal_lag <- .sp_temporal_lag(temporal_values)
  if (!is.null(time) && isTRUE(time$estimate) &&
      !any(temporal_lag[upper.tri(temporal_lag)] > 0)) {
    stop("Estimated temporal scale requires at least two distinct observed times.",
         call. = FALSE)
  }
  if (!control$estimate && any(vapply(
    list(control$spatial_variance, control$scale, control$nugget_variance),
    is.null, logical(1)))) {
    stop("Fixed spatial fitting requires `spatial_variance`, `scale`, and `nugget_variance`.",
         call. = FALSE)
  }

  distance <- .sp_dist(coordinates)
  positive_distances <- distance[upper.tri(distance) & distance > 0]
  if (!length(positive_distances)) {
    stop("At least two distinct coordinate locations are required.", call. = FALSE)
  }

  ordinary_residual <- stats::lm.fit(x, y_centered)$residuals
  response_variance <- mean((ordinary_residual - mean(ordinary_residual))^2)
  if (!is.finite(response_variance) || response_variance <= 0) {
    response_variance <- max(mean((y_centered - mean(y_centered))^2), 1)
  }
  initial <- c(
    spatial_variance = if (is.null(control$spatial_variance)) 0.8 * response_variance else control$spatial_variance,
    scale = if (is.null(control$scale)) stats::median(positive_distances) else control$scale,
    nugget_variance = if (is.null(control$nugget_variance)) 0.2 * response_variance else control$nugget_variance
  )
  random_initial <- .sp_resolve_random_variance(
    control$random_variance,
    if (is.null(random_specification)) character() else
      random_specification$variables,
    control$estimate,
    default = 0.1 * max(response_variance, 1e-8)
  )
  random_slope_parameterization <- .sp_random_slope_parameterization(
    random_slope_specification
  )
  parameter_initial <- c(
    initial,
    if (length(random_initial)) stats::setNames(
      random_initial, paste0("random::", names(random_initial))
    ) else numeric(),
    if (!is.null(time) && isTRUE(time$estimate))
      c(temporal_scale = time$scale) else numeric()
  )
  optimisation_initial <- c(
    stats::setNames(log(parameter_initial), paste0("log_", names(parameter_initial))),
    random_slope_parameterization$values
  )

  optimisation <- NULL
  if (control$estimate) {
    d_min <- min(positive_distances)
    d_max <- max(positive_distances)
    variance_floor <- max(response_variance, 1) * 1e-10
    lower <- c(log(c(variance_floor, d_min * 1e-4, variance_floor,
                   rep(variance_floor, length(random_initial)),
                   if (!is.null(time) && isTRUE(time$estimate)) {
                     positive_time <- temporal_lag[upper.tri(temporal_lag) &
                       temporal_lag > 0]
                     if (length(positive_time)) min(positive_time) * 1e-4 else
                       time$scale * 1e-4
                   } else numeric())),
      unlist(lapply(random_slope_parameterization$components, function(details) {
        index <- details$index
        vapply(seq_len(nrow(index)), function(k) {
          if (index[k, 1L] == index[k, 2L]) log(sqrt(variance_floor)) else
            -sqrt(max(response_variance, 1) * 1e10)
        }, numeric(1L))
      }), use.names = FALSE))
    upper <- c(log(c(max(response_variance, 1) * 1e10, d_max * 1e4,
                   max(response_variance, 1) * 1e10,
                   rep(max(response_variance, 1) * 1e10,
                       length(random_initial)),
                   if (!is.null(time) && isTRUE(time$estimate)) {
                     positive_time <- temporal_lag[upper.tri(temporal_lag) &
                       temporal_lag > 0]
                     if (length(positive_time)) max(positive_time) * 1e4 else
                       time$scale * 1e4
                   } else numeric())),
      unlist(lapply(random_slope_parameterization$components, function(details) {
        index <- details$index
        vapply(seq_len(nrow(index)), function(k) {
          if (index[k, 1L] == index[k, 2L]) {
            log(sqrt(max(response_variance, 1) * 1e10))
          } else sqrt(max(response_variance, 1) * 1e10)
        }, numeric(1L))
      }), use.names = FALSE))
    base_count <- length(parameter_initial)
    objective <- function(transformed_pars) {
      values <- stats::setNames(exp(transformed_pars[seq_len(base_count)]),
                                names(parameter_initial))
      pars <- values[names(initial)]
      random_values <- if (length(random_initial)) {
        answer <- values[paste0("random::", names(random_initial))]
        names(answer) <- names(random_initial)
        answer
      } else numeric()
      temporal_scale <- if (!is.null(time) && isTRUE(time$estimate)) {
        values[["temporal_scale"]]
      } else if (!is.null(time)) time$scale else NULL
      current_random_slope <- .sp_random_slope_update(
        random_slope_specification,
        transformed_pars[base_count + seq_along(random_slope_parameterization$values)],
        random_slope_parameterization
      )
      state <- .sp_gaussian_state(pars, distance, x, y_centered,
                                  control$smoothness,
                                  random_z = if (is.null(random_specification))
                                    NULL else random_specification$z,
                                  random_variance = random_values,
                                  random_column_group = if (is.null(random_specification))
                                    NULL else random_specification$column_group,
                                  random_slope_specification = current_random_slope,
                                  temporal_lag = temporal_lag,
                                  temporal_scale = temporal_scale)
      if (is.null(state) || !is.finite(state$log_likelihood)) return(1e20)
      -state$log_likelihood
    }
    optimisation <- stats::optim(
      optimisation_initial, objective, method = "L-BFGS-B", lower = lower,
      upper = upper, control = list(maxit = control$maxit,
                                    factr = control$reltol / .Machine$double.eps,
                                    trace = if (verbose) 1 else 0),
      hessian = TRUE
    )
    fitted_parameters <- stats::setNames(
      exp(optimisation$par[seq_len(length(parameter_initial))]),
                                         names(parameter_initial))
    covariance <- fitted_parameters[names(initial)]
    random_variance <- if (length(random_initial)) {
      answer <- fitted_parameters[paste0("random::", names(random_initial))]
      names(answer) <- names(random_initial)
      answer
    } else numeric()
    temporal_scale <- if (!is.null(time) && isTRUE(time$estimate)) {
      fitted_parameters[["temporal_scale"]]
    } else if (!is.null(time)) time$scale else NULL
    random_slope_specification <- .sp_random_slope_update(
      random_slope_specification,
      optimisation$par[length(parameter_initial) +
        seq_along(random_slope_parameterization$values)],
      random_slope_parameterization
    )
    if (optimisation$convergence != 0L) {
      warning(sprintf("Covariance optimisation ended with convergence code %d: %s",
                      optimisation$convergence, optimisation$message),
              call. = FALSE)
    }
  } else {
    covariance <- initial
    random_variance <- random_initial
    temporal_scale <- if (!is.null(time)) time$scale else NULL
  }
  state <- .sp_gaussian_state(covariance, distance, x, y_centered,
                              control$smoothness,
                              random_z = if (is.null(random_specification))
                                NULL else random_specification$z,
                              random_variance = random_variance,
                              random_column_group = if (is.null(random_specification))
                                NULL else random_specification$column_group,
                              random_slope_specification = random_slope_specification,
                              temporal_lag = temporal_lag,
                              temporal_scale = temporal_scale)
  if (is.null(state)) stop("The fitted covariance or mean system is not positive definite.", call. = FALSE)
  names(state$beta) <- colnames(x)
  beta_covariance <- chol2inv(state$chol_information)
  dimnames(beta_covariance) <- list(colnames(x), colnames(x))
  covariance_vcov <- matrix(NA_real_, 3L, 3L,
                            dimnames = list(names(covariance), names(covariance)))
  random_variance_vcov <- matrix(
    NA_real_, length(random_variance), length(random_variance),
    dimnames = list(names(random_variance), names(random_variance))
  )
  variance_parameter_vcov <- matrix(
    NA_real_, length(parameter_initial), length(parameter_initial),
    dimnames = list(names(parameter_initial), names(parameter_initial))
  )
  random_slope_values <- .sp_random_slope_parameters(
    random_slope_specification
  )
  random_slope_vcov <- matrix(
    NA_real_, length(random_slope_values), length(random_slope_values),
    dimnames = list(names(random_slope_values), names(random_slope_values))
  )
  transformed_parameter_vcov <- matrix(
    NA_real_, length(optimisation_initial), length(optimisation_initial),
    dimnames = list(names(optimisation_initial), names(optimisation_initial))
  )
  if (control$estimate && !is.null(optimisation$hessian) &&
      all(is.finite(optimisation$hessian))) {
    chol_hessian <- tryCatch(chol(optimisation$hessian),
                             error = function(e) NULL)
    if (!is.null(chol_hessian)) {
      log_vcov <- chol2inv(chol_hessian)
      dimnames(log_vcov) <- list(names(optimisation_initial),
                                 names(optimisation_initial))
      transformed_parameter_vcov <- log_vcov
      fitted_values <- c(covariance, if (length(random_variance))
        stats::setNames(random_variance,
                        paste0("random::", names(random_variance))) else
          numeric(), if (!is.null(time) && isTRUE(time$estimate))
            c(temporal_scale = temporal_scale) else numeric())
      jacobian <- diag(unname(fitted_values), nrow = length(fitted_values))
      base_log_vcov <- log_vcov[seq_along(fitted_values),
                                seq_along(fitted_values), drop = FALSE]
      variance_parameter_vcov <- jacobian %*% base_log_vcov %*% jacobian
      dimnames(variance_parameter_vcov) <- list(names(fitted_values),
                                                names(fitted_values))
      covariance_vcov <- variance_parameter_vcov[names(covariance),
                                                  names(covariance), drop = FALSE]
      dimnames(covariance_vcov) <- list(names(covariance), names(covariance))
      if (length(random_variance)) {
        random_names <- paste0("random::", names(random_variance))
        random_variance_vcov <- variance_parameter_vcov[random_names,
          random_names, drop = FALSE]
        dimnames(random_variance_vcov) <- list(names(random_variance),
                                               names(random_variance))
      }
      if (length(random_slope_parameterization$values)) {
        slope_indices <- length(parameter_initial) +
          seq_along(random_slope_parameterization$values)
        slope_jacobian <- .sp_random_slope_natural_jacobian(
          random_slope_specification
        )
        estimated_rows <- rowSums(abs(slope_jacobian)) > 0
        random_slope_vcov[estimated_rows, estimated_rows] <-
          slope_jacobian[estimated_rows, , drop = FALSE] %*%
          log_vcov[slope_indices, slope_indices, drop = FALSE] %*%
          t(slope_jacobian[estimated_rows, , drop = FALSE])
      }
    }
  }

  structure(list(
    call = call, formula = formula, terms = terms_object, spatial = TRUE,
    contrasts = attr(x, "contrasts"), xlevels = stats::.getXlevels(terms_object, mf),
    coords = coords, coordinate_names = colnames(coordinates),
    coordinate_source = coordinate_info$source, crs = coordinate_info$crs,
    distance_units = distance_units, family = family,
    model = mf, x = x, y = y, offset = model_offset,
    data_used = data[row_index, , drop = FALSE], control = control,
    external_offset = "offset" %in% names(call),
    coordinates = coordinates, covariance = covariance,
    random = random, random_specification = random_specification,
    random_variance = random_variance,
    random_slope = .sp_random_slope_export(random_slope_specification),
    random_slope_specification = random_slope_specification,
    random_variance_vcov = random_variance_vcov,
    random_slope_vcov = random_slope_vcov,
    random_slope_parameterization = random_slope_parameterization,
    temporal = time, temporal_values = temporal_values,
    temporal_scale = temporal_scale,
    temporal_scale_estimated = !is.null(time) && isTRUE(time$estimate),
    temporal_scale_se = if (!is.null(time) && isTRUE(time$estimate) &&
      is.finite(variance_parameter_vcov["temporal_scale", "temporal_scale"]) &&
      variance_parameter_vcov["temporal_scale", "temporal_scale"] >= 0)
      sqrt(variance_parameter_vcov["temporal_scale", "temporal_scale"]) else
        NA_real_,
    variance_parameter_vcov = variance_parameter_vcov,
    transformed_parameter_vcov = transformed_parameter_vcov,
    smoothness = control$smoothness, covariance_estimated = control$estimate,
    covariance_vcov = covariance_vcov,
    coefficients = state$beta, coefficient_vcov = beta_covariance,
    residual_mean = state$residual, fitted_state = state,
    log_likelihood = state$log_likelihood, optimisation = optimisation,
    nobs = length(y), rank = ncol(x)
  ), class = "spmodels_fit")
}
