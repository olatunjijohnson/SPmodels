#' Simulate from a fitted spatial model
#'
#' Generate unconditional response replicates at the fitted observation sites,
#' using the fitted mean, family and covariance parameters.
#'
#' @param object A fitted `spmodels_fit` object.
#' @param nsim Number of simulations.
#' @param seed Optional integer seed. The caller's random-number state is
#'   restored after simulation.
#' @param ... Reserved for method compatibility.
#'
#' @return A data frame with one simulated response per column.
#' @export
simulate.spmodels_fit <- function(object, nsim = 1L, seed = NULL, ...) {
  if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) ||
      !is.finite(nsim) || nsim < 1 || nsim != as.integer(nsim)) {
    stop("`nsim` must be one positive integer.", call. = FALSE)
  }
  seed_scope <- .sp_seed_scope(seed)
  if (!is.null(seed_scope)) {
    on.exit(.sp_restore_seed(seed_scope), add = TRUE)
    set.seed(seed_scope$seed)
  }

  z <- matrix(stats::rnorm(object$nobs * as.integer(nsim)),
              nrow = object$nobs)
  if (object$family$family == "gaussian") {
    correlated_error <- crossprod(object$fitted_state$chol_v, z)
    mean_response <- as.vector(object$x %*% object$coefficients + object$offset)
    simulations <- sweep(correlated_error, 1L, mean_response, "+")
  } else {
    latent <- crossprod(object$fitted_state$chol_sigma, z)
    fixed_predictor <- as.vector(
      object$x %*% object$coefficients + object$offset
    )
    eta <- sweep(latent, 1L, fixed_predictor, "+")
    if (object$family$family == "binomial") {
      probability <- stats::plogis(eta)
      if (!is.null(object$intervention_multiplier)) {
        probability <- sweep(probability, 1L,
                             object$intervention_multiplier, "*")
      }
      simulations <- matrix(stats::rbinom(
        length(probability),
        size = rep(object$trials, as.integer(nsim)),
        prob = as.vector(probability)
      ), nrow = object$nobs)
    } else {
      intensity <- exp(eta)
      if (any(!is.finite(intensity))) {
        stop("A simulated Poisson intensity overflowed; inspect fitted parameters.",
             call. = FALSE)
      }
      simulations <- matrix(stats::rpois(length(intensity), intensity),
                            nrow = object$nobs)
    }
  }
  simulations <- as.data.frame(simulations)
  names(simulations) <- paste0("sim_", seq_len(ncol(simulations)))
  row.names(simulations) <- row.names(object$model)
  attr(simulations, "seed") <- seed
  simulations
}

#' Empirical residual variogram from a fitted spatial model
#'
#' Bin pairwise semivariances of either mean-model residuals or conditional
#' residuals. The returned object also contains the fitted model variogram at
#' each bin midpoint.
#'
#' @param object A fitted `spmodels_fit` object.
#' @param bins Positive number of equal-width distance bins.
#' @param max_distance Optional positive maximum pair distance. The default is
#'   the largest observed pair distance.
#' @param residual_type For Gaussian fits, residuals after the regression mean
#'   (`"mean"`) or conditional spatial prediction (`"conditional"`). For
#'   binomial/Poisson fits, `"pearson"`, `"deviance"`, or
#'   `"randomized_quantile"`. The family-specific first choice is the default.
#' @param seed Optional local seed for randomized-quantile residuals.
#'
#' @return An object of class `spmodels_variogram` containing a `values` data
#'   frame and metadata.
#' @export
sp_variogram <- function(object, bins = 10L, max_distance = NULL,
                         residual_type = NULL, seed = NULL) {
  if (!inherits(object, "spmodels_fit")) {
    stop("`object` must be a fitted `spmodels_fit` model.", call. = FALSE)
  }
  if (!isTRUE(object$spatial)) {
    stop("A spatial variogram requires a model fitted with `spatial = TRUE`.",
         call. = FALSE)
  }
  gaussian <- object$family$family == "gaussian"
  choices <- if (gaussian) c("mean", "conditional") else
    c("randomized_quantile", "pearson", "deviance")
  if (is.null(residual_type)) residual_type <- choices[[1L]]
  if (!is.character(residual_type) || length(residual_type) != 1L ||
      is.na(residual_type) || !residual_type %in% choices) {
    stop(sprintf("`residual_type` must be one of: %s.",
                 paste(choices, collapse = ", ")), call. = FALSE)
  }
  if (!is.numeric(bins) || length(bins) != 1L || is.na(bins) ||
      !is.finite(bins) || bins < 1 || bins != as.integer(bins)) {
    stop("`bins` must be one positive integer.", call. = FALSE)
  }
  distance_matrix <- .sp_dist(object$coordinates)
  keep_triangle <- lower.tri(distance_matrix)
  pair_distance <- distance_matrix[keep_triangle]
  if (is.null(max_distance)) max_distance <- max(pair_distance)
  if (!is.numeric(max_distance) || length(max_distance) != 1L ||
      is.na(max_distance) || !is.finite(max_distance) || max_distance <= 0) {
    stop("`max_distance` must be one finite positive number.", call. = FALSE)
  }
  residual_value <- if (gaussian && residual_type == "mean") {
    object$residual_mean
  } else if (gaussian) {
    stats::residuals(object)
  } else {
    diagnostics <- sp_diagnostics(object, seed = seed)
    switch(residual_type,
      randomized_quantile = diagnostics$randomized_quantile_residual,
      pearson = diagnostics$pearson_residual,
      deviance = diagnostics$deviance_residual)
  }
  difference <- outer(residual_value, residual_value, "-")
  pair_semivariance <- 0.5 * difference[keep_triangle]^2
  retained <- pair_distance <= max_distance
  if (!any(retained)) stop("No location pairs fall within `max_distance`.", call. = FALSE)
  breaks <- seq(0, max_distance, length.out = as.integer(bins) + 1L)
  group <- cut(pair_distance[retained], breaks = breaks,
               include.lowest = TRUE, right = TRUE)
  levels_group <- levels(group)
  count <- as.integer(table(group))
  empirical <- as.numeric(tapply(pair_semivariance[retained], group, mean))
  midpoint <- (breaks[-1L] + breaks[-length(breaks)]) / 2
  if (!gaussian) {
    fitted <- rep(NA_real_, length(midpoint))
  } else if (is.null(object$temporal)) {
    fitted <- object$covariance[["nugget_variance"]] +
      object$covariance[["spatial_variance"]] *
      (1 - matern_cor(midpoint, object$covariance[["scale"]],
                      object$smoothness))
  } else {
    pair_time_lag <- .sp_temporal_lag(object$temporal_values)[keep_triangle]
    fitted_pair <- object$covariance[["nugget_variance"]] +
      object$covariance[["spatial_variance"]] *
      (1 - spacetime_cor(pair_distance, pair_time_lag,
        object$covariance[["scale"]], object$smoothness,
        object$temporal_scale))
    fitted <- as.numeric(tapply(fitted_pair[retained], group, mean))
  }
  if (gaussian && !is.null(object$random_specification)) {
    random_covariance <- .sp_random_covariance(
      object$random_specification$z, object$random_specification$z,
      object$random_specification$column_group, object$random_variance
    )
    random_semivariance <- 0.5 * outer(diag(random_covariance),
      diag(random_covariance), "+") - random_covariance
    grouped_random <- tapply(random_semivariance[keep_triangle][retained],
                             group, mean)
    fitted <- fitted + as.numeric(grouped_random)
  }
  if (gaussian && !is.null(object$random_slope_specification)) {
    random_covariance <- .sp_random_slope_covariance(
      object$random_slope_specification
    )
    random_semivariance <- 0.5 * outer(diag(random_covariance),
      diag(random_covariance), "+") - random_covariance
    grouped_random <- tapply(random_semivariance[keep_triangle][retained],
                             group, mean)
    fitted <- fitted + as.numeric(grouped_random)
  }
  values <- data.frame(
    bin = levels_group, distance = midpoint, pairs = count,
    semivariance = empirical, fitted = unname(fitted),
    stringsAsFactors = FALSE
  )
  structure(list(values = values, residual_type = residual_type,
                 family = object$family$family, seed = seed,
                 distance_units = object$distance_units,
                 max_distance = max_distance, call = match.call()),
            class = "spmodels_variogram")
}

#' @export
plot.spmodels_variogram <- function(x, ...) {
  values <- x$values
  use <- values$pairs > 0L & is.finite(values$semivariance)
  graphics::plot(values$distance[use], values$semivariance[use],
                 xlab = paste0("Distance (", x$distance_units, ")"),
                 ylab = "Semivariance", ...)
  if (any(is.finite(values$fitted))) {
    graphics::lines(values$distance, values$fitted, col = 2L, lwd = 2L)
    graphics::legend("topleft", legend = c("Empirical", "Fitted model"),
      pch = c(1, NA), lty = c(NA, 1), col = c(1, 2), bty = "n")
  } else {
    graphics::legend("topleft", legend = "Empirical residual variogram",
                     pch = 1, bty = "n")
  }
  invisible(x)
}

.sp_binned_semivariance <- function(values, triangle, retained, group,
                                    bins) {
  difference <- outer(values, values, "-")
  pair_semivariance <- 0.5 * difference[triangle]^2
  grouped <- tapply(pair_semivariance[retained], group, mean)
  answer <- rep(NA_real_, bins)
  if (length(grouped)) {
    answer[as.integer(names(grouped))] <- as.numeric(grouped)
  }
  answer
}

#' Monte Carlo envelope for a residual variogram
#'
#' Compare an observed residual variogram with either a permutation
#' null of no location-value association or simulations from the fitted spatial
#' model. For discrete families, parametric simulation repeats the declared
#' fitting route before residual calculation. Pointwise envelopes and a global
#' curve-discrepancy p-value are returned from the same simulated curves.
#'
#' @param object A fitted spatial `spmodels_fit` object.
#' @param nsim Positive number of Monte Carlo curves.
#' @param method Permute residuals across fixed locations, or simulate from the
#'   fitted model while keeping fitted parameters fixed.
#' @param level Pointwise envelope probability strictly between zero and one.
#' @param bins,max_distance,residual_type Passed to [sp_variogram()].
#' @param seed Optional integer seed local to this call.
#'
#' @return An object of class `spmodels_variogram_envelope` containing observed
#'   and simulated curves, pointwise limits and a global Monte Carlo p-value.
#' @export
sp_variogram_envelope <- function(
    object, nsim = 999L, method = c("permutation", "parametric"),
    level = 0.95, bins = 10L, max_distance = NULL,
    residual_type = NULL, seed = NULL) {
  if (!inherits(object, "spmodels_fit") || !isTRUE(object$spatial)) {
    stop("`object` must be a fitted spatial model.", call. = FALSE)
  }
  if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) ||
      !is.finite(nsim) || nsim < 1 || nsim != as.integer(nsim)) {
    stop("`nsim` must be one positive integer.", call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      !is.finite(level) || level <= 0 || level >= 1) {
    stop("`level` must be one number strictly between zero and one.",
         call. = FALSE)
  }
  method <- match.arg(method)
  gaussian <- object$family$family == "gaussian"
  choices <- if (gaussian) c("mean", "conditional") else
    c("randomized_quantile", "pearson", "deviance")
  if (is.null(residual_type)) residual_type <- choices[[1L]]
  if (!is.character(residual_type) || length(residual_type) != 1L ||
      is.na(residual_type) || !residual_type %in% choices) {
    stop(sprintf("`residual_type` must be one of: %s.",
                 paste(choices, collapse = ", ")), call. = FALSE)
  }
  if (!gaussian && method == "parametric" && isTRUE(object$external_offset)) {
    stop("Non-Gaussian parametric envelopes require offsets to be written in the formula.",
         call. = FALSE)
  }
  observed <- sp_variogram(object, bins = bins, max_distance = max_distance,
                           residual_type = residual_type, seed = seed)
  bins <- nrow(observed$values)
  distance_matrix <- .sp_dist(object$coordinates)
  triangle <- lower.tri(distance_matrix)
  pair_distance <- distance_matrix[triangle]
  retained <- pair_distance <= observed$max_distance
  breaks <- seq(0, observed$max_distance, length.out = bins + 1L)
  group <- cut(pair_distance[retained], breaks = breaks,
               include.lowest = TRUE, right = TRUE, labels = FALSE)
  observed_residual <- if (gaussian) {
    if (residual_type == "mean") object$residual_mean else
      stats::residuals(object)
  } else {
    diagnostics <- sp_diagnostics(object, seed = seed)
    switch(residual_type,
      randomized_quantile = diagnostics$randomized_quantile_residual,
      pearson = diagnostics$pearson_residual,
      deviance = diagnostics$deviance_residual)
  }
  seed_scope <- .sp_seed_scope(seed)
  if (!is.null(seed_scope)) {
    on.exit(.sp_restore_seed(seed_scope), add = TRUE)
    set.seed(seed_scope$seed)
  }
  refit_convergence <- NULL
  refit_warnings <- NULL
  residual_simulations <- if (method == "permutation") {
    replicate(as.integer(nsim), sample(observed_residual, replace = FALSE))
  } else if (gaussian) {
    response_simulations <- as.matrix(stats::simulate(
      object, nsim = as.integer(nsim)
    ))
    mean_response <- as.vector(object$x %*% object$coefficients +
                                 object$offset)
    mean_residuals <- sweep(response_simulations, 1L, mean_response, "-")
    if (residual_type == "mean") {
      mean_residuals
    } else {
      spatial_covariance <- .sp_covariance(
        distance_matrix, object$covariance, object$smoothness,
        diagonal = FALSE,
        random_z = if (is.null(object$random_specification)) NULL else
          object$random_specification$z,
        random_variance = object$random_variance,
        random_column_group = if (is.null(object$random_specification)) NULL else
          object$random_specification$column_group,
        random_slope_specification = object$random_slope_specification,
        temporal_lag = .sp_temporal_lag(object$temporal_values),
        temporal_scale = object$temporal_scale
      )
      mean_residuals - spatial_covariance %*%
        .sp_solve_chol(object$fitted_state$chol_v, mean_residuals)
    }
  } else {
    response_simulations <- as.matrix(stats::simulate(
      object, nsim = as.integer(nsim)
    ))
    refit_convergence <- rep(NA_integer_, as.integer(nsim))
    refit_warnings <- vector("list", as.integer(nsim))
    vapply(seq_len(as.integer(nsim)), function(simulation) {
      warning_text <- character()
      fitted_simulation <- tryCatch(withCallingHandlers({
          simulated_data <- .sp_simulation_data(
            object, response_simulations[, simulation]
          )
          simulation_weights <- if (
            object$family$family == "binomial" &&
            identical(object$response_encoding, "proportion")
          ) object$trials else NULL
          refit <- sp_fit(
            object$formula, simulated_data, family = object$family,
            coords = if (inherits(simulated_data, "sf")) NULL else object$coords,
            spatial = object$spatial, time = object$temporal,
            random = object$random, random_slope = object$random_slope,
            intervention = object$intervention,
            weights = simulation_weights, control = object$control
          )
          if (identical(object$likelihood_method, "Monte Carlo ML")) {
            refit <- sp_mcml(refit, control = object$mcml$control)
          }
          refit
        }, warning = function(condition) {
          warning_text <<- c(warning_text, conditionMessage(condition))
          invokeRestart("muffleWarning")
        }), error = identity)
      if (inherits(fitted_simulation, "error")) {
        stop(sprintf("Parametric envelope refit %d failed: %s",
                     simulation, conditionMessage(fitted_simulation)),
             call. = FALSE)
      }
      refit_convergence[[simulation]] <- if (
        is.null(fitted_simulation$optimisation)
      ) 0L else fitted_simulation$optimisation$convergence
      refit_warnings[[simulation]] <- warning_text
      diagnostics <- sp_diagnostics(fitted_simulation)
      switch(residual_type,
        randomized_quantile = diagnostics$randomized_quantile_residual,
        pearson = diagnostics$pearson_residual,
        deviance = diagnostics$deviance_residual)
    }, numeric(object$nobs))
  }
  if (length(refit_convergence) &&
      any(is.na(refit_convergence) | refit_convergence != 0L)) {
    warning(sprintf(
      "%d of %d parametric envelope refits reported non-zero convergence.",
      sum(is.na(refit_convergence) | refit_convergence != 0L),
      length(refit_convergence)), call. = FALSE)
  }
  simulated <- vapply(seq_len(as.integer(nsim)), function(simulation) {
    .sp_binned_semivariance(
      residual_simulations[, simulation], triangle, retained, group, bins
    )
  }, numeric(bins))
  if (is.null(dim(simulated))) simulated <- matrix(simulated, ncol = 1L)
  alpha <- (1 - level) / 2
  curve_summary <- function(values, operation) {
    values <- values[is.finite(values)]
    if (!length(values)) return(NA_real_)
    operation(values)
  }
  lower <- apply(simulated, 1L, curve_summary, operation = function(values)
    stats::quantile(values, probs = alpha, names = FALSE, type = 8))
  upper <- apply(simulated, 1L, curve_summary, operation = function(values)
    stats::quantile(values, probs = 1 - alpha, names = FALSE, type = 8))
  centre <- apply(simulated, 1L, curve_summary, operation = mean)
  spread <- apply(simulated, 1L, curve_summary, operation = function(values)
    if (length(values) > 1L) stats::sd(values) else NA_real_)
  usable <- is.finite(observed$values$semivariance) & is.finite(centre) &
    is.finite(spread) & spread > 0
  if (!any(usable)) {
    stop("Simulated variograms have no variable, populated bins.", call. = FALSE)
  }
  observed_statistic <- sum(
    ((observed$values$semivariance[usable] - centre[usable]) /
       spread[usable])^2
  )
  simulated_statistic <- apply(simulated[usable, , drop = FALSE], 2L,
    function(curve) sum(((curve - centre[usable]) / spread[usable])^2))
  p_value <- (1 + sum(simulated_statistic >= observed_statistic)) /
    (as.integer(nsim) + 1)
  values <- observed$values
  values$lower <- lower
  values$upper <- upper
  values$simulation_mean <- centre
  structure(list(
    values = values, simulations = simulated, method = method,
    level = level, p_value = p_value,
    observed_statistic = observed_statistic,
    simulated_statistic = simulated_statistic,
    residual_type = residual_type,
    family = object$family$family,
    refit_convergence = refit_convergence,
    refit_warnings = refit_warnings,
    distance_units = object$distance_units,
    max_distance = observed$max_distance, nsim = as.integer(nsim),
    seed = seed, call = match.call()
  ), class = "spmodels_variogram_envelope")
}

#' @export
print.spmodels_variogram_envelope <- function(x, ...) {
  cat("SPmodels residual-variogram Monte Carlo envelope\n")
  cat("Method:", x$method, " Simulations:", x$nsim,
      " Pointwise level:", format(x$level), "\n")
  cat("Global curve p-value:", format.pval(x$p_value, digits = 4), "\n")
  invisible(x)
}

#' @export
plot.spmodels_variogram_envelope <- function(x, ...) {
  values <- x$values
  use <- values$pairs > 0L & is.finite(values$semivariance) &
    is.finite(values$lower) & is.finite(values$upper)
  graphics::plot(values$distance[use], values$semivariance[use], type = "n",
    ylim = range(c(values$lower[use], values$upper[use],
                   values$semivariance[use]), finite = TRUE),
    xlab = paste0("Distance (", x$distance_units, ")"),
    ylab = "Semivariance", ...)
  graphics::polygon(
    c(values$distance[use], rev(values$distance[use])),
    c(values$lower[use], rev(values$upper[use])),
    col = "grey85", border = NA
  )
  graphics::lines(values$distance[use], values$semivariance[use],
                  lwd = 2L)
  graphics::points(values$distance[use], values$semivariance[use], pch = 16)
  if (any(is.finite(values$fitted))) {
    graphics::lines(values$distance, values$fitted, col = 2L, lwd = 2L)
    graphics::legend("topleft",
      legend = c("Observed", "Pointwise envelope", "Fitted model"),
      col = c(1, "grey70", 2), lwd = c(2, 8, 2), bty = "n")
  } else {
    graphics::legend("topleft",
      legend = c("Observed", "Pointwise envelope"),
      col = c(1, "grey70"), lwd = c(2, 8), bty = "n")
  }
  invisible(x)
}

#' Fitted spatial-model diagnostics
#'
#' Calculate mean-model residuals, conditional residuals and decorrelated
#' Gaussian residuals. If the fitted covariance model is adequate, the
#' decorrelated residuals should be approximately independent standard normal
#' values, subject to parameter estimation. For binomial and Poisson fits,
#' return response, Pearson and signed deviance residuals at the fitted latent
#' mode. Randomized quantile residuals use a caller-local seed.
#'
#' @param object A fitted `spmodels_fit` object.
#' @param seed Optional integer seed for non-Gaussian randomized quantile
#'   residuals. The caller's random-number state is restored.
#'
#' @return An object of class `spmodels_diagnostics`.
#' @export
sp_diagnostics <- function(object, seed = NULL) {
  if (!inherits(object, "spmodels_fit")) {
    stop("`object` must be a fitted `spmodels_fit` model.", call. = FALSE)
  }
  if (object$family$family != "gaussian") {
    fitted_mean <- as.vector(object$latent_mode$conditional_mean)
    response_residual <- object$y - fitted_mean
    if (object$family$family == "binomial") {
      probability <- ifelse(object$trials > 0,
                            fitted_mean / object$trials, 0)
      conditional_variance <- object$trials * probability * (1 - probability)
      deviance_squared <- object$family$dev.resids(
        ifelse(object$trials > 0, object$y / object$trials, 0),
        probability, object$trials
      )
      cdf_lower <- stats::pbinom(object$y - 1, object$trials, probability)
      cdf_upper <- stats::pbinom(object$y, object$trials, probability)
    } else {
      conditional_variance <- fitted_mean
      deviance_squared <- object$family$dev.resids(
        object$y, fitted_mean, rep(1, object$nobs)
      )
      cdf_lower <- stats::ppois(object$y - 1, fitted_mean)
      cdf_upper <- stats::ppois(object$y, fitted_mean)
    }
    pearson <- rep(NA_real_, object$nobs)
    positive_variance <- conditional_variance > 0
    pearson[positive_variance] <- response_residual[positive_variance] /
      sqrt(conditional_variance[positive_variance])
    deviance <- sign(response_residual) * sqrt(pmax(0, deviance_squared))
    seed_scope <- .sp_seed_scope(seed)
    if (!is.null(seed_scope)) {
      on.exit(.sp_restore_seed(seed_scope), add = TRUE)
      set.seed(seed_scope$seed)
    }
    randomized_probability <- stats::runif(object$nobs, cdf_lower, cdf_upper)
    probability_floor <- .Machine$double.eps
    randomized_probability <- pmin(1 - probability_floor,
                                   pmax(probability_floor,
                                        randomized_probability))
    randomized_quantile <- stats::qnorm(randomized_probability)
    return(structure(list(
      response_residual = response_residual,
      pearson_residual = pearson,
      deviance_residual = deviance,
      randomized_quantile_residual = randomized_quantile,
      randomized_probability = randomized_probability,
      fitted_mean = fitted_mean,
      coordinates = object$coordinates,
      distance_units = object$distance_units,
      family = object$family$family,
      call = match.call()
    ), class = "spmodels_diagnostics"))
  }
  decorrelated <- as.vector(forwardsolve(
    t(object$fitted_state$chol_v), object$residual_mean
  ))
  structure(list(
    mean_residual = object$residual_mean,
    conditional_residual = stats::residuals(object),
    decorrelated_residual = decorrelated,
    family = "gaussian",
    coordinates = object$coordinates,
    distance_units = object$distance_units,
    call = match.call()
  ), class = "spmodels_diagnostics")
}

#' @export
plot.spmodels_diagnostics <- function(x, ...) {
  old_par <- graphics::par(mfrow = c(1, 2))
  on.exit(graphics::par(old_par), add = TRUE)
  displayed <- if (identical(x$family, "gaussian"))
    x$conditional_residual else x$pearson_residual
  graphics::plot(seq_along(displayed), displayed,
                 xlab = "Observation index",
                 ylab = if (identical(x$family, "gaussian"))
                   "Conditional residual" else "Pearson residual",
                 main = if (identical(x$family, "gaussian"))
                   "Conditional residuals" else "Pearson residuals", ...)
  graphics::abline(h = 0, lty = 2)
  qq_values <- if (identical(x$family, "gaussian")) {
    x$decorrelated_residual
  } else if (!is.null(x$randomized_quantile_residual)) {
    x$randomized_quantile_residual
  } else x$deviance_residual
  stats::qqnorm(qq_values, main = if (identical(x$family, "gaussian"))
    "Decorrelated residuals" else if (!is.null(x$randomized_quantile_residual))
      "Randomized quantile residuals" else "Deviance residuals")
  stats::qqline(qq_values, col = 2L)
  invisible(x)
}
