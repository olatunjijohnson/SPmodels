.sp_neighbor_graph <- function(coordinates, neighbors,
                               ordering = c("input", "maximin"),
                               observation_order = NULL) {
  ordering <- match.arg(ordering)
  count <- nrow(coordinates)
  if (!is.numeric(neighbors) || length(neighbors) != 1L || is.na(neighbors) ||
      !is.finite(neighbors) || neighbors < 1 || neighbors != as.integer(neighbors)) {
    stop("`neighbors` must be one positive integer.", call. = FALSE)
  }
  neighbors <- as.integer(neighbors)
  if (is.null(observation_order)) {
    observation_order <- if (ordering == "input") seq_len(count) else
      .sp_maximin_prediction_order(coordinates)
  } else {
    if (!is.numeric(observation_order) || length(observation_order) != count ||
        anyNA(observation_order) || any(!is.finite(observation_order)) ||
        any(observation_order != as.integer(observation_order)) ||
        !setequal(as.integer(observation_order), seq_len(count))) {
      stop("`observation_order` must be a permutation of observation rows.",
           call. = FALSE)
    }
    observation_order <- as.integer(observation_order)
    ordering <- "user"
  }
  graph <- vector("list", count)
  for (step in seq_len(count)) {
    current <- observation_order[[step]]
    predecessors <- if (step == 1L) integer() else
      observation_order[seq_len(step - 1L)]
    if (length(predecessors) > neighbors) {
      distance <- sqrt(rowSums(sweep(
        coordinates[predecessors, , drop = FALSE], 2L,
        coordinates[current, ], "-"
      )^2))
      predecessors <- predecessors[order(distance, seq_along(distance))[
        seq_len(neighbors)
      ]]
    }
    graph[[current]] <- predecessors
  }
  list(neighbors = graph, order = observation_order, ordering = ordering,
       budget = min(neighbors, max(0L, count - 1L)))
}

.sp_neighbor_gaussian_state <- function(covariance, coordinates, graph,
                                         x, y_centered, smoothness) {
  count <- length(y_centered)
  transformed_y <- numeric(count)
  transformed_x <- matrix(0, count, ncol(x))
  conditional_variance <- numeric(count)
  coefficients <- vector("list", count)
  for (current in graph$order) {
    selected <- graph$neighbors[[current]]
    if (!length(selected)) {
      coefficient <- numeric()
      variance <- covariance[["spatial_variance"]] +
        covariance[["nugget_variance"]]
    } else {
      neighbor_distance <- .sp_dist(coordinates[selected, , drop = FALSE])
      neighbor_covariance <- covariance[["spatial_variance"]] *
        .sp_matern_cor_unchecked(neighbor_distance, covariance[["scale"]],
                                 smoothness)
      diag(neighbor_covariance) <- diag(neighbor_covariance) +
        covariance[["nugget_variance"]]
      cross_distance <- .sp_dist(
        coordinates[selected, , drop = FALSE],
        coordinates[current, , drop = FALSE]
      )
      cross <- as.vector(covariance[["spatial_variance"]] *
        .sp_matern_cor_unchecked(cross_distance, covariance[["scale"]],
                                 smoothness))
      coefficient <- tryCatch(
        .sp_psd_solve(neighbor_covariance, cross),
        error = function(e) NULL
      )
      if (is.null(coefficient)) return(NULL)
      variance <- covariance[["spatial_variance"]] +
        covariance[["nugget_variance"]] - sum(cross * coefficient)
    }
    scale <- max(1, sum(covariance[c("spatial_variance",
                                     "nugget_variance")]))
    tolerance <- 1000 * .Machine$double.eps * scale *
      max(1L, length(selected))
    if (!is.finite(variance) || variance <= tolerance) return(NULL)
    coefficients[[current]] <- coefficient
    conditional_variance[[current]] <- variance
    transformed_y[[current]] <- y_centered[[current]] -
      if (length(selected)) sum(coefficient * y_centered[selected]) else 0
    transformed_x[current, ] <- x[current, ] -
      if (length(selected)) colSums(x[selected, , drop = FALSE] * coefficient) else 0
  }
  weighted_x <- transformed_x / sqrt(conditional_variance)
  information <- crossprod(weighted_x)
  chol_information <- tryCatch(chol(information), error = function(e) NULL)
  if (is.null(chol_information)) return(NULL)
  weighted_y <- transformed_y / conditional_variance
  beta <- .sp_solve_chol(chol_information,
                         crossprod(transformed_x, weighted_y))
  transformed_residual <- as.vector(transformed_y - transformed_x %*% beta)
  log_likelihood <- -0.5 * sum(
    log(2 * pi) + log(conditional_variance) +
      transformed_residual^2 / conditional_variance
  )
  fitted_values <- numeric(count)
  residual <- as.vector(y_centered - x %*% beta)
  for (current in graph$order) {
    selected <- graph$neighbors[[current]]
    fitted_values[[current]] <- as.vector(x[current, , drop = FALSE] %*% beta) +
      if (length(selected)) sum(coefficients[[current]] * residual[selected]) else 0
  }
  list(beta = as.vector(beta), information = information,
       chol_information = chol_information,
       coefficient_vcov = chol2inv(chol_information),
       conditional_coefficients = coefficients,
       conditional_variance = conditional_variance,
       transformed_x = transformed_x, transformed_y = transformed_y,
       transformed_residual = transformed_residual,
       residual = residual, fitted_centered = fitted_values,
       log_likelihood = log_likelihood)
}

.sp_coordinate_distance_range <- function(coordinates) {
  minimum <- Inf
  maximum <- 0
  if (nrow(coordinates) > 1L) {
    for (index in 2:nrow(coordinates)) {
      distance <- sqrt(rowSums(sweep(
        coordinates[seq_len(index - 1L), , drop = FALSE], 2L,
        coordinates[index, ], "-"
      )^2))
      positive <- distance[distance > 0]
      if (length(positive)) {
        minimum <- min(minimum, positive)
        maximum <- max(maximum, positive)
      }
    }
  }
  if (!is.finite(minimum)) {
    stop("At least two distinct coordinate locations are required.",
         call. = FALSE)
  }
  c(minimum = minimum, maximum = maximum)
}

#' Fit an approximate Gaussian nearest-predecessor model
#'
#' Fit a Gaussian identity-link geostatistical regression by multiplying
#' ordered conditional densities that retain at most `neighbors` earlier
#' observations. The calculation uses bounded covariance blocks. Keeping all
#' predecessors recovers the ordinary dense Gaussian likelihood up to numerical
#' precision. This explicit opt-in fitter does not silently replace [sp_fit()].
#'
#' @param formula,data A model formula and data frame.
#' @param coords A one-sided formula selecting two numeric coordinate columns.
#' @param neighbors Positive maximum predecessor count.
#' @param ordering Input-row or deterministic maximin order.
#' @param observation_order Optional explicit permutation of observation rows.
#' @param subset,na.action,offset Standard model-frame arguments.
#' @param control Settings from [sp_control()]. Grouped and temporal covariance
#'   components are not supported by this initial approximation.
#' @param verbose Show optimization tracing?
#'
#' @return An object of class `spmodels_neighbor_fit` with the graph, local
#'   conditional coefficients and variances retained for audit.
#' @export
sp_fit_neighbor <- function(
    formula, data, coords, neighbors = 30L,
    ordering = c("input", "maximin"), observation_order = NULL,
    subset = NULL, na.action = stats::na.fail, offset = NULL,
    control = sp_control(), verbose = FALSE) {
  if (!is.data.frame(data)) {
    stop("This implementation requires `data` to be a data frame.", call. = FALSE)
  }
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula.", call. = FALSE)
  }
  if (!inherits(control, "spmodels_control")) {
    stop("`control` must be created by `sp_control()`.", call. = FALSE)
  }
  if (!is.null(control$random_variance)) {
    stop("Nearest-predecessor fitting does not yet support grouped covariance components.",
         call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }
  ordering <- match.arg(ordering)
  call <- match.call()
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
  coordinate_info <- .sp_coordinates(coords, data, row_index)
  coordinates <- coordinate_info$values
  distance_units <- control$distance_units
  if (is.null(distance_units)) distance_units <- coordinate_info$inferred_units
  if (is.null(distance_units)) {
    stop("Supply explicit coordinate units with `sp_control(distance_units = ...)`.",
         call. = FALSE)
  }
  distance_range <- .sp_coordinate_distance_range(coordinates)
  graph <- .sp_neighbor_graph(coordinates, neighbors, ordering,
                              observation_order)
  ordinary_residual <- stats::lm.fit(x, y_centered)$residuals
  response_variance <- mean((ordinary_residual - mean(ordinary_residual))^2)
  if (!is.finite(response_variance) || response_variance <= 0) {
    response_variance <- max(mean((y_centered - mean(y_centered))^2), 1)
  }
  covariance_initial <- c(
    spatial_variance = if (is.null(control$spatial_variance))
      0.8 * response_variance else control$spatial_variance,
    scale = if (is.null(control$scale))
      sqrt(distance_range[["minimum"]] * distance_range[["maximum"]]) else
        control$scale,
    nugget_variance = if (is.null(control$nugget_variance))
      0.2 * response_variance else control$nugget_variance
  )
  if (!control$estimate && any(vapply(
    list(control$spatial_variance, control$scale, control$nugget_variance),
    is.null, logical(1)))) {
    stop("Fixed nearest-predecessor fitting requires all three covariance parameters.",
         call. = FALSE)
  }
  optimisation <- NULL
  if (control$estimate) {
    variance_floor <- max(response_variance, 1) * 1e-10
    lower <- log(c(variance_floor, distance_range[["minimum"]] * 1e-4,
                   variance_floor))
    upper <- log(c(max(response_variance, 1) * 1e10,
                   distance_range[["maximum"]] * 1e4,
                   max(response_variance, 1) * 1e10))
    objective <- function(log_parameters) {
      covariance <- stats::setNames(exp(log_parameters),
                                    names(covariance_initial))
      state <- .sp_neighbor_gaussian_state(
        covariance, coordinates, graph, x, y_centered, control$smoothness
      )
      if (is.null(state) || !is.finite(state$log_likelihood)) 1e20 else
        -state$log_likelihood
    }
    optimisation <- stats::optim(
      log(covariance_initial), objective, method = "L-BFGS-B",
      lower = lower, upper = upper,
      control = list(maxit = control$maxit,
                     factr = control$reltol / .Machine$double.eps,
                     trace = if (verbose) 1 else 0), hessian = TRUE
    )
    covariance <- stats::setNames(exp(optimisation$par),
                                  names(covariance_initial))
    if (optimisation$convergence != 0L) {
      warning(sprintf("Neighbor-likelihood optimisation ended with convergence code %d: %s",
                      optimisation$convergence, optimisation$message),
              call. = FALSE)
    }
  } else covariance <- covariance_initial
  state <- .sp_neighbor_gaussian_state(
    covariance, coordinates, graph, x, y_centered, control$smoothness
  )
  if (is.null(state)) {
    stop("The fitted local covariance or mean system is not positive definite.",
         call. = FALSE)
  }
  names(state$beta) <- colnames(x)
  dimnames(state$coefficient_vcov) <- list(colnames(x), colnames(x))
  covariance_vcov <- matrix(NA_real_, 3L, 3L,
    dimnames = list(names(covariance), names(covariance)))
  if (control$estimate && !is.null(optimisation$hessian) &&
      all(is.finite(optimisation$hessian))) {
    hessian_root <- tryCatch(chol(optimisation$hessian),
                             error = function(e) NULL)
    if (!is.null(hessian_root)) {
      log_vcov <- chol2inv(hessian_root)
      covariance_vcov <- diag(covariance) %*% log_vcov %*% diag(covariance)
      dimnames(covariance_vcov) <- list(names(covariance), names(covariance))
    }
  }
  structure(list(
    call = call, formula = formula, terms = terms_object, spatial = TRUE,
    contrasts = attr(x, "contrasts"),
    xlevels = stats::.getXlevels(terms_object, mf), coords = coords,
    coordinate_names = colnames(coordinates),
    coordinate_source = coordinate_info$source, crs = coordinate_info$crs,
    distance_units = distance_units, family = stats::gaussian(), model = mf,
    x = x, y = y, offset = model_offset,
    data_used = data[row_index, , drop = FALSE], control = control,
    external_offset = "offset" %in% names(call), coordinates = coordinates,
    covariance = covariance, covariance_vcov = covariance_vcov,
    covariance_estimated = control$estimate, smoothness = control$smoothness,
    random_variance = numeric(), random_variance_vcov = matrix(numeric(), 0, 0),
    random_slope_specification = NULL, temporal = NULL, intervention = NULL,
    coefficients = state$beta, coefficient_vcov = state$coefficient_vcov,
    residual_mean = state$residual, fitted_state = state,
    fitted_values = state$fitted_centered + model_offset,
    log_likelihood = state$log_likelihood, optimisation = optimisation,
    neighbor_graph = graph$neighbors, observation_order = graph$order,
    ordering = graph$ordering, neighbors = graph$budget,
    likelihood_method = "Gaussian nearest-predecessor ML",
    nobs = length(y), rank = ncol(x)
  ), class = "spmodels_neighbor_fit")
}

#' @export
coef.spmodels_neighbor_fit <- function(object, ...) object$coefficients

#' @export
vcov.spmodels_neighbor_fit <- function(object, ...) object$coefficient_vcov

#' @export
logLik.spmodels_neighbor_fit <- function(object, ...) {
  value <- object$log_likelihood
  attr(value, "nobs") <- object$nobs
  attr(value, "df") <- object$rank + if (object$covariance_estimated) 3L else 0L
  class(value) <- "logLik"
  value
}

#' @export
fitted.spmodels_neighbor_fit <- function(object, ...) object$fitted_values

#' @export
residuals.spmodels_neighbor_fit <- function(object, ...) {
  object$y - object$fitted_values
}

#' @export
print.spmodels_neighbor_fit <- function(x, ...) {
  cat("Gaussian nearest-predecessor geostatistical regression\n")
  cat("Formula: "); print(x$formula)
  cat("Observations:", x$nobs, "\n")
  cat("Predecessors:", x$neighbors, "(", x$ordering, "order)\n")
  cat("Distance units:", x$distance_units, "\n")
  cat("Covariance parameters", if (x$covariance_estimated)
    "(estimated):\n" else "(fixed):\n")
  print(c(x$covariance, smoothness = x$smoothness))
  invisible(x)
}

#' @export
summary.spmodels_neighbor_fit <- function(object, ...) {
  coefficient_se <- sqrt(diag(object$coefficient_vcov))
  structure(list(
    call = object$call,
    coefficients = cbind(Estimate = object$coefficients,
      `Std. Error` = coefficient_se,
      `z value` = object$coefficients / coefficient_se),
    covariance_table = cbind(Estimate = object$covariance,
      `Std. Error` = sqrt(diag(object$covariance_vcov))),
    covariance_estimated = object$covariance_estimated,
    smoothness = object$smoothness, log_likelihood = stats::logLik(object),
    nobs = object$nobs, neighbors = object$neighbors,
    ordering = object$ordering, distance_units = object$distance_units,
    convergence = if (is.null(object$optimisation)) NA_integer_ else
      object$optimisation$convergence
  ), class = "summary.spmodels_neighbor_fit")
}

#' @export
print.summary.spmodels_neighbor_fit <- function(x, ...) {
  cat("Gaussian nearest-predecessor geostatistical regression\n\nCall:\n")
  print(x$call)
  cat("\nMean coefficients:\n")
  stats::printCoefmat(x$coefficients)
  cat("\nCovariance parameters", if (x$covariance_estimated)
    "(estimated):\n" else "(fixed):\n")
  print(x$covariance_table)
  cat("Fixed smoothness:", format(x$smoothness), "\n")
  cat("Predecessors:", x$neighbors, "(", x$ordering, "order)\n")
  if (x$covariance_estimated) cat("Optimiser convergence code:", x$convergence, "\n")
  cat("Log-likelihood:", format(as.numeric(x$log_likelihood)), "\n")
  cat("Distance units:", x$distance_units, "\n")
  invisible(x)
}

#' @export
predict.spmodels_neighbor_fit <- function(
    object, newdata = NULL, se.fit = FALSE,
    target = c("latent", "observation"), ...) {
  target <- match.arg(target)
  if (!is.logical(se.fit) || length(se.fit) != 1L || is.na(se.fit)) {
    stop("`se.fit` must be TRUE or FALSE.", call. = FALSE)
  }
  if (is.null(newdata)) {
    answer <- object$fitted_values
    if (!se.fit) return(answer)
    variance <- object$fitted_state$conditional_variance
    if (target == "latent") {
      variance <- variance - object$covariance[["nugget_variance"]]
    }
    variance <- variance + rowSums(
      (object$fitted_state$transformed_x %*% object$coefficient_vcov) *
        object$fitted_state$transformed_x
    )
    return(list(fit = answer,
      se.fit = sqrt(pmax(0, variance)), target = target,
      approximation = "nearest-predecessor"))
  }
  if (!is.data.frame(newdata)) stop("`newdata` must be a data frame.", call. = FALSE)
  terms_new <- stats::delete.response(object$terms)
  model_newdata <- if (inherits(newdata, "sf")) {
    if (!requireNamespace("sf", quietly = TRUE)) {
      stop("Package `sf` is required for `sf` input.", call. = FALSE)
    }
    sf::st_drop_geometry(newdata)
  } else newdata
  mf_new <- stats::model.frame(terms_new, model_newdata,
                               na.action = stats::na.fail,
                               xlev = object$xlevels)
  x_new <- stats::model.matrix(terms_new, mf_new,
                               contrasts.arg = object$contrasts)
  coordinate_info <- .sp_coordinates(object$coords, newdata,
                                      expected_crs = object$crs)
  coordinates_new <- coordinate_info$values
  if (nrow(coordinates_new) != nrow(x_new)) {
    stop("Coordinate and predictor rows in `newdata` do not align.", call. = FALSE)
  }
  offset_new <- stats::model.offset(mf_new)
  if (is.null(offset_new)) offset_new <- rep(0, nrow(x_new))
  fit <- variance <- numeric(nrow(x_new))
  budget <- if (object$neighbors >= object$nobs - 1L) object$nobs else
    object$neighbors
  for (index in seq_len(nrow(x_new))) {
    distance <- sqrt(rowSums(sweep(object$coordinates, 2L,
                                   coordinates_new[index, ], "-")^2))
    selected <- order(distance, seq_along(distance))[seq_len(budget)]
    neighbor_distance <- .sp_dist(object$coordinates[selected, , drop = FALSE])
    neighbor_covariance <- object$covariance[["spatial_variance"]] *
      .sp_matern_cor_unchecked(neighbor_distance, object$covariance[["scale"]],
                               object$smoothness)
    diag(neighbor_covariance) <- diag(neighbor_covariance) +
      object$covariance[["nugget_variance"]]
    cross <- as.vector(object$covariance[["spatial_variance"]] *
      .sp_matern_cor_unchecked(distance[selected], object$covariance[["scale"]],
                               object$smoothness))
    coefficient <- .sp_psd_solve(neighbor_covariance, cross)
    fit[[index]] <- as.vector(x_new[index, , drop = FALSE] %*%
      object$coefficients) + offset_new[[index]] +
      sum(coefficient * object$residual_mean[selected])
    mean_gradient <- as.vector(x_new[index, , drop = FALSE] -
      crossprod(coefficient, object$x[selected, , drop = FALSE]))
    variance[[index]] <- object$covariance[["spatial_variance"]] -
      sum(cross * coefficient) +
      sum((mean_gradient %*% object$coefficient_vcov) * mean_gradient)
    if (target == "observation") {
      variance[[index]] <- variance[[index]] +
        object$covariance[["nugget_variance"]]
    }
  }
  if (!se.fit) return(fit)
  tolerance <- 1000 * .Machine$double.eps *
    max(1, sum(object$covariance[c("spatial_variance", "nugget_variance")]))
  if (any(variance < -tolerance)) {
    warning("A negative local prediction variance was produced; check conditioning.",
            call. = FALSE)
  }
  list(fit = fit, se.fit = sqrt(pmax(0, variance)), target = target,
       approximation = "nearest-predecessor")
}

#' @export
confint.spmodels_neighbor_fit <- function(object, parm, level = 0.95, ...) {
  if (missing(parm)) {
    confint.spmodels_fit(object, level = level, ...)
  } else {
    confint.spmodels_fit(object, parm = parm, level = level, ...)
  }
}
