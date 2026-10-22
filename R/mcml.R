.sp_log_mean_exp <- function(values) {
  maximum <- max(values)
  if (!is.finite(maximum)) return(-Inf)
  maximum + log(mean(exp(values - maximum)))
}

.sp_mcml_intervention <- function(parameters, object) {
  if (!isTRUE(object$intervention_estimated)) return(object$intervention)
  p <- ncol(object$x)
  covariance_count <- length(.sp_fitted_covariance_names(object))
  random_count <- if (object$covariance_estimated)
    length(object$random_variance) else 0L
  random_slope_count <- length(
    .sp_random_slope_parameterization(object$random_slope_specification)$values
  )
  temporal_count <- if (isTRUE(object$temporal_scale_estimated)) 1L else 0L
  indices <- p + covariance_count + random_count + random_slope_count +
    temporal_count + seq_len(2L)
  current <- object$intervention
  current$alpha <- unname(stats::plogis(parameters[indices[[1L]]]))
  current$decay <- unname(exp(parameters[indices[[2L]]]))
  current
}

.sp_mcml_temporal_scale <- function(parameters, object) {
  if (!isTRUE(object$temporal_scale_estimated)) return(object$temporal_scale)
  p <- ncol(object$x)
  covariance_count <- length(.sp_fitted_covariance_names(object))
  random_count <- if (object$covariance_estimated)
    length(object$random_variance) else 0L
  random_slope_count <- length(
    .sp_random_slope_parameterization(object$random_slope_specification)$values
  )
  unname(exp(parameters[p + covariance_count + random_count +
    random_slope_count + 1L]))
}

.sp_mcml_random_slope <- function(parameters, object) {
  parameterization <- .sp_random_slope_parameterization(
    object$random_slope_specification
  )
  count <- length(parameterization$values)
  if (!count) return(object$random_slope_specification)
  p <- ncol(object$x)
  covariance_count <- length(.sp_fitted_covariance_names(object))
  random_count <- if (object$covariance_estimated)
    length(object$random_variance) else 0L
  start <- p + covariance_count + random_count
  .sp_random_slope_update(
    object$random_slope_specification,
    parameters[start + seq_len(count)], parameterization
  )
}

.sp_mcml_random_variance <- function(parameters, object) {
  if (!length(object$random_variance) || !object$covariance_estimated) {
    return(object$random_variance)
  }
  p <- ncol(object$x)
  indices <- p + length(.sp_fitted_covariance_names(object)) +
    seq_along(object$random_variance)
  stats::setNames(exp(parameters[indices]), names(object$random_variance))
}

.sp_mcml_log_weights <- function(parameters, object, samples, log_proposal,
                                 distance) {
  p <- ncol(object$x)
  beta <- parameters[seq_len(p)]
  covariance_names <- .sp_fitted_covariance_names(object)
  covariance <- if (object$covariance_estimated) {
    answer <- object$covariance
    answer[covariance_names] <-
      exp(parameters[p + seq_along(covariance_names)])
    answer
  } else object$covariance
  random_variance <- .sp_mcml_random_variance(parameters, object)
  random_slope_specification <- .sp_mcml_random_slope(parameters, object)
  temporal_scale <- .sp_mcml_temporal_scale(parameters, object)
  sigma <- tryCatch(
    .sp_covariance(distance, covariance, object$smoothness, diagonal = TRUE,
                   temporal_lag = .sp_temporal_lag(object$temporal_values),
                   temporal_scale = temporal_scale,
                   random_z = if (is.null(object$random_specification)) NULL else
                     object$random_specification$z,
                   random_variance = random_variance,
                   random_column_group = if (is.null(object$random_specification))
                     NULL else object$random_specification$column_group,
                   random_slope_specification = random_slope_specification),
    error = function(e) NULL
  )
  if (is.null(sigma)) return(NULL)
  chol_sigma <- tryCatch(chol(sigma), error = function(e) NULL)
  if (is.null(chol_sigma)) return(NULL)

  fixed_predictor <- as.vector(object$x %*% beta + object$offset)
  eta <- sweep(samples, 2L, fixed_predictor, "+")
  b <- nrow(samples)
  n <- ncol(samples)
  if (object$family$family == "binomial") {
    current_intervention <- .sp_mcml_intervention(parameters, object)
    probability_multiplier <- .sp_intervention_multiplier(
      current_intervention, object$data_used
    )
    if (is.null(probability_multiplier)) {
      probability_multiplier <- rep(1, n)
    }
    likelihood <- stats::dbinom(
      matrix(object$y, b, n, byrow = TRUE),
      size = matrix(object$trials, b, n, byrow = TRUE),
      prob = sweep(stats::plogis(eta), 2L, probability_multiplier, "*"),
      log = TRUE
    )
  } else {
    likelihood <- stats::dpois(
      matrix(object$y, b, n, byrow = TRUE), lambda = exp(eta), log = TRUE
    )
  }
  log_response <- rowSums(likelihood)
  if (isTRUE(object$intervention_estimated) &&
      !is.null(current_intervention$alpha_penalty)) {
    log_response <- log_response + stats::dbeta(
      current_intervention$alpha,
      current_intervention$alpha_penalty[["shape1"]],
      current_intervention$alpha_penalty[["shape2"]], log = TRUE
    )
  }
  precision_samples <- .sp_solve_chol(chol_sigma, t(samples))
  quadratic <- rowSums(samples * t(precision_samples))
  log_determinant <- 2 * sum(log(diag(chol_sigma)))
  log_prior <- -0.5 * (n * log(2 * pi) + log_determinant + quadratic)
  log_response + log_prior - log_proposal
}

.sp_mcml_diagnostics <- function(log_weights, component = NULL) {
  maximum <- max(log_weights)
  scaled <- exp(log_weights - maximum)
  scaled_sum <- sum(scaled)
  normalised <- scaled / scaled_sum
  relative_mcse <- if (length(scaled) <= 1L) {
    NA_real_
  } else if (is.null(component) || length(unique(component)) == 1L) {
    stats::sd(scaled) / (sqrt(length(scaled)) * mean(scaled))
  } else {
    groups <- split(scaled, component)
    probability <- lengths(groups) / length(scaled)
    variance <- sum(probability^2 * vapply(groups, function(value) {
      if (length(value) > 1L) stats::var(value) / length(value) else NA_real_
    }, numeric(1L)))
    if (is.finite(variance)) sqrt(variance) / mean(scaled) else NA_real_
  }
  list(
    effective_sample_size = 1 / sum(normalised^2),
    effective_sample_fraction = 1 / (length(normalised) * sum(normalised^2)),
    maximum_weight = max(normalised),
    relative_mcse = relative_mcse,
    log_likelihood_mcse = relative_mcse
  )
}

.sp_mcml_proposal <- function(mean, base_covariance, scales, nsim) {
  base_chol <- tryCatch(chol(base_covariance), error = function(e) NULL)
  if (is.null(base_chol)) {
    stop("The latent-posterior proposal is not positive definite.",
         call. = FALSE)
  }
  n <- length(mean)
  component_count <- length(scales)
  component <- rep(seq_len(component_count), length.out = nsim)
  counts <- tabulate(component, nbins = component_count)
  probability <- counts / nsim
  standard_draws <- matrix(stats::rnorm(nsim * n), nrow = nsim, ncol = n)
  samples <- standard_draws %*% base_chol
  samples <- samples * scales[component]
  samples <- sweep(samples, 2L, mean, "+")
  centered <- sweep(samples, 2L, mean, "-")
  base_standardized <- forwardsolve(t(base_chol), t(centered))
  base_quadratic <- colSums(base_standardized^2)
  base_log_determinant <- 2 * sum(log(diag(base_chol)))
  component_log_density <- vapply(scales, function(scale) {
    -0.5 * (n * log(2 * pi) + base_log_determinant +
      2 * n * log(scale) + base_quadratic / scale^2)
  }, numeric(nsim))
  if (component_count == 1L) {
    component_log_density <- matrix(component_log_density, ncol = 1L)
  }
  weighted_log_density <- sweep(component_log_density, 2L,
                                log(probability), "+")
  row_maximum <- apply(weighted_log_density, 1L, max)
  log_density <- row_maximum + log(rowSums(exp(
    sweep(weighted_log_density, 1L, row_maximum, "-")
  )))
  list(
    mean = mean, base_covariance = base_covariance, base_chol = base_chol,
    covariances = lapply(scales, function(scale) {
      base_covariance * scale^2
    }),
    component = component, counts = counts, probability = probability,
    samples = samples, log_density = log_density
  )
}

.sp_mcml_adapt_proposal <- function(pilot, log_weights, original_mean,
                                    original_covariance, shrinkage) {
  maximum <- max(log_weights)
  weights <- exp(log_weights - maximum)
  weights <- weights / sum(weights)
  weighted_mean <- colSums(pilot$samples * weights)
  centered <- sweep(pilot$samples, 2L, weighted_mean, "-")
  denominator <- 1 - sum(weights^2)
  if (!is.finite(denominator) || denominator <= 100 * .Machine$double.eps) {
    stop("Pilot importance weights cannot identify an adapted covariance.",
         call. = FALSE)
  }
  weighted_centered <- sweep(centered, 1L, sqrt(weights), "*")
  weighted_covariance <- crossprod(weighted_centered) / denominator
  mean <- shrinkage * original_mean + (1 - shrinkage) * weighted_mean
  covariance <- shrinkage * original_covariance +
    (1 - shrinkage) * weighted_covariance
  covariance <- (covariance + t(covariance)) / 2
  if (is.null(tryCatch(chol(covariance), error = function(e) NULL))) {
    stop("The regularized adapted proposal is not positive definite.",
         call. = FALSE)
  }
  list(mean = mean, covariance = covariance,
       weighted_mean = weighted_mean,
       weighted_covariance = weighted_covariance)
}

#' Control Monte Carlo maximum-likelihood correction
#'
#' Configure the normalized importance distribution, optional discarded-pilot
#' adaptation and deterministic optimisation used by [sp_mcml()].
#'
#' @param nsim Number of independent importance samples.
#' @param maxit Positive maximum number of optimisation iterations.
#' @param reltol Positive relative optimisation tolerance.
#' @param proposal_scale One or more distinct positive multipliers for the
#'   standard deviation of the local Gaussian latent-posterior approximation.
#'   Multiple values create a deterministic normalized mixture proposal.
#' @param ess_warning Fraction of `nsim` below which [sp_mcml()] warns that the
#'   importance distribution needs improvement.
#' @param adapt Use an independent pilot sample to estimate a regularized
#'   proposal centre and covariance before drawing the retained main sample?
#' @param pilot_nsim Number of discarded pilot importance draws when
#'   `adapt = TRUE`.
#' @param adapt_shrinkage Weight in `(0, 1]` assigned to the original Laplace
#'   proposal moments; the remaining weight is assigned to pilot weighted
#'   moments.
#'
#' @return An object of class `spmodels_mcml_control`.
#' @export
sp_mcml_control <- function(nsim = 2000L, maxit = 200L, reltol = 1e-6,
                            proposal_scale = 1,
                            ess_warning = 0.1, adapt = FALSE,
                            pilot_nsim = 500L, adapt_shrinkage = 0.5) {
  integer_scalar <- function(value, name, minimum) {
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        !is.finite(value) || value < minimum || value != as.integer(value)) {
      stop(sprintf("`%s` must be one integer not less than %d.",
                   name, minimum), call. = FALSE)
    }
    as.integer(value)
  }
  positive_scalar <- function(value, name) {
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        !is.finite(value) || value <= 0) {
      stop(sprintf("`%s` must be one finite positive number.", name),
           call. = FALSE)
    }
    as.numeric(value)
  }
  nsim <- integer_scalar(nsim, "nsim", 100L)
  maxit <- integer_scalar(maxit, "maxit", 1L)
  reltol <- positive_scalar(reltol, "reltol")
  if (!is.numeric(proposal_scale) || !length(proposal_scale) ||
      anyNA(proposal_scale) || any(!is.finite(proposal_scale)) ||
      any(proposal_scale <= 0)) {
    stop("`proposal_scale` must contain finite positive numbers.",
         call. = FALSE)
  }
  proposal_scale <- as.numeric(proposal_scale)
  if (anyDuplicated(proposal_scale)) {
    stop("`proposal_scale` values must be distinct.", call. = FALSE)
  }
  if (!is.numeric(ess_warning) || length(ess_warning) != 1L ||
      is.na(ess_warning) || !is.finite(ess_warning) || ess_warning <= 0 ||
      ess_warning >= 1) {
    stop("`ess_warning` must be one number strictly between zero and one.",
         call. = FALSE)
  }
  if (!is.logical(adapt) || length(adapt) != 1L || is.na(adapt)) {
    stop("`adapt` must be TRUE or FALSE.", call. = FALSE)
  }
  pilot_nsim <- integer_scalar(pilot_nsim, "pilot_nsim", 100L)
  if (!is.numeric(adapt_shrinkage) || length(adapt_shrinkage) != 1L ||
      is.na(adapt_shrinkage) || !is.finite(adapt_shrinkage) ||
      adapt_shrinkage <= 0 || adapt_shrinkage > 1) {
    stop("`adapt_shrinkage` must be one number greater than zero and not greater than one.",
         call. = FALSE)
  }
  structure(list(
    nsim = nsim, maxit = maxit, reltol = reltol,
    proposal_scale = proposal_scale, ess_warning = as.numeric(ess_warning),
    adapt = adapt, pilot_nsim = pilot_nsim,
    adapt_shrinkage = as.numeric(adapt_shrinkage)
  ), class = "spmodels_mcml_control")
}

#' Refine a non-Gaussian fit by Monte Carlo maximum likelihood
#'
#' Replace the Laplace objective of a fitted binomial-logit or Poisson-log
#' spatial model with an importance-sampling estimate of its marginal
#' likelihood. The proposal starts from the scaled local Gaussian
#' latent-posterior approximation at the supplied fit. Optionally, a discarded
#' independent pilot regularizes its centre and covariance before the retained
#' sample is drawn. Common random numbers make the objective deterministic
#' during a call.
#'
#' @param object A binomial-logit or Poisson-log `spmodels_fit`, initially
#'   fitted by [sp_fit()].
#' @param control Settings returned by [sp_mcml_control()].
#' @param seed Optional integer seed local to this call.
#' @param verbose Logical scalar controlling optimisation tracing.
#'
#' @return An updated `spmodels_fit`. Its `mcml` component records proposal
#'   draws, final log weights, effective sample size, maximum normalized weight
#'   and delta-method Monte Carlo standard error.
#' @export
sp_mcml <- function(object, control = sp_mcml_control(), seed = NULL,
                    verbose = FALSE) {
  if (!inherits(object, "spmodels_fit") ||
      object$family$family == "gaussian") {
    stop("`object` must be a fitted binomial or Poisson `spmodels_fit` model.",
         call. = FALSE)
  }
  if (!inherits(control, "spmodels_mcml_control")) {
    stop("`control` must be created by `sp_mcml_control()`.", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }
  seed_scope <- .sp_seed_scope(seed)
  if (!is.null(seed_scope)) {
    on.exit(.sp_restore_seed(seed_scope), add = TRUE)
    set.seed(seed_scope$seed)
  }

  proposal_mean <- object$latent_mode$mode
  proposal_base_covariance <- object$latent_mode$covariance
  original_proposal_mean <- proposal_mean
  original_proposal_covariance <- proposal_base_covariance

  p <- length(object$coefficients)
  random_slope_parameterization <- .sp_random_slope_parameterization(
    object$random_slope_specification
  )
  random_slope_count <- length(random_slope_parameterization$values)
  covariance_names <- .sp_fitted_covariance_names(object)
  covariance_count <- length(covariance_names)
  parameter_start <- c(object$coefficients,
    if (covariance_count) log(object$covariance[covariance_names]) else NULL,
    if (object$covariance_estimated && length(object$random_variance))
      log(object$random_variance) else NULL,
    random_slope_parameterization$values,
    if (isTRUE(object$temporal_scale_estimated))
      log(object$temporal_scale) else NULL,
    if (isTRUE(object$intervention_estimated))
      c(stats::qlogis(object$intervention$alpha),
        log(object$intervention$decay)) else NULL)
  names(parameter_start) <- c(names(object$coefficients),
    if (covariance_count) paste0("log_", covariance_names) else NULL,
    if (object$covariance_estimated && length(object$random_variance))
      paste0("log_random::", names(object$random_variance)) else NULL,
    random_slope_parameterization$names,
    if (isTRUE(object$temporal_scale_estimated))
      "log_temporal_scale" else NULL,
    if (isTRUE(object$intervention_estimated))
      c("logit_alpha", "log_decay") else NULL)
  distance <- .sp_dist(object$coordinates)
  pilot <- adaptation <- pilot_diagnostics <- NULL
  adaptation_reason <- NULL
  if (isTRUE(control$adapt)) {
    pilot <- .sp_mcml_proposal(
      proposal_mean, proposal_base_covariance, control$proposal_scale,
      control$pilot_nsim
    )
    pilot_log_weights <- .sp_mcml_log_weights(
      parameter_start, object, pilot$samples, pilot$log_density, distance
    )
    if (is.null(pilot_log_weights) || any(!is.finite(pilot_log_weights))) {
      stop("The MCML pilot sample did not produce finite importance weights.",
           call. = FALSE)
    }
    pilot_diagnostics <- .sp_mcml_diagnostics(
      pilot_log_weights, pilot$component
    )
    if (pilot_diagnostics$effective_sample_fraction < control$ess_warning) {
      warning(sprintf(
        "Adaptive-pilot ESS is %.1f%% of `pilot_nsim`; retaining the original proposal moments.",
        100 * pilot_diagnostics$effective_sample_fraction), call. = FALSE)
      adaptation_reason <- "pilot ESS fraction below ess_warning"
    } else {
      adaptation <- .sp_mcml_adapt_proposal(
        pilot, pilot_log_weights, proposal_mean, proposal_base_covariance,
        control$adapt_shrinkage
      )
      proposal_mean <- adaptation$mean
      proposal_base_covariance <- adaptation$covariance
      adaptation_reason <- "pilot ESS threshold met"
    }
  }
  proposal <- .sp_mcml_proposal(
    proposal_mean, proposal_base_covariance, control$proposal_scale,
    control$nsim
  )
  samples <- proposal$samples
  log_proposal <- proposal$log_density
  proposal_component <- proposal$component
  proposal_counts <- proposal$counts
  proposal_probability <- proposal$probability
  proposal_covariances <- proposal$covariances
  positive_distances <- distance[upper.tri(distance) & distance > 0]
  lower <- rep(-Inf, length(parameter_start))
  upper <- rep(Inf, length(parameter_start))
  if (covariance_count) {
    covariance_indices <- p + seq_len(covariance_count)
    lower[covariance_indices] <- log(c(
      spatial_variance = 1e-7,
      scale = min(positive_distances) * 1e-4,
      nugget_variance = 1e-7
    )[covariance_names])
    upper[covariance_indices] <- log(c(
      spatial_variance = 1e5,
      scale = max(positive_distances) * 1e4,
      nugget_variance = 1e5
    )[covariance_names])
    if (length(object$random_variance)) {
      random_indices <- p + covariance_count +
        seq_along(object$random_variance)
      lower[random_indices] <- log(1e-7)
      upper[random_indices] <- log(1e5)
    }
  }
  if (random_slope_count) {
    slope_start <- p + covariance_count +
      (if (object$covariance_estimated) length(object$random_variance) else 0L)
    cursor <- 0L
    for (details in random_slope_parameterization$components) {
      for (k in seq_len(nrow(details$index))) {
        cursor <- cursor + 1L
        diagonal <- details$index[k, 1L] == details$index[k, 2L]
        lower[slope_start + cursor] <- if (diagonal) log(sqrt(1e-7)) else
          -sqrt(1e5)
        upper[slope_start + cursor] <- if (diagonal) log(sqrt(1e5)) else
          sqrt(1e5)
      }
    }
  }
  if (isTRUE(object$temporal_scale_estimated)) {
    temporal_index <- p + covariance_count +
      (if (object$covariance_estimated) length(object$random_variance) else 0L) +
      random_slope_count + 1L
    temporal_lag <- .sp_temporal_lag(object$temporal_values)
    positive_time <- temporal_lag[upper.tri(temporal_lag) & temporal_lag > 0]
    lower[temporal_index] <- log(if (length(positive_time))
      min(positive_time) * 1e-4 else object$temporal_scale * 1e-4)
    upper[temporal_index] <- log(if (length(positive_time))
      max(positive_time) * 1e4 else object$temporal_scale * 1e4)
  }
  intervention_indices <- if (isTRUE(object$intervention_estimated)) {
    p + covariance_count +
      (if (object$covariance_estimated) length(object$random_variance) else 0L) +
      random_slope_count +
      (if (isTRUE(object$temporal_scale_estimated)) 1L else 0L) +
      seq_len(2L)
  } else integer()
  if (length(intervention_indices)) {
    survey_times <- object$data_used[[object$intervention$time]]
    time_span <- diff(range(c(survey_times, object$intervention$history$time)))
    if (!is.finite(time_span) || time_span <= 0) time_span <- 1
    lower[intervention_indices] <- c(stats::qlogis(1e-6),
                                     log(time_span * 1e-6))
    upper[intervention_indices] <- c(stats::qlogis(1 - 1e-6),
                                     log(time_span * 1e4))
  }
  objective <- function(parameters) {
    log_weights <- tryCatch(
      .sp_mcml_log_weights(parameters, object, samples, log_proposal,
                           distance),
      error = function(e) NULL
    )
    if (is.null(log_weights)) return(1e20)
    value <- .sp_log_mean_exp(log_weights)
    if (!is.finite(value)) 1e20 else -value
  }
  optimisation <- stats::optim(
    parameter_start, objective, method = "L-BFGS-B", lower = lower,
    upper = upper, hessian = TRUE,
    control = list(maxit = control$maxit,
      factr = control$reltol / .Machine$double.eps,
      trace = if (verbose) 1 else 0)
  )
  if (optimisation$convergence != 0L) {
    warning(sprintf("MCML optimisation ended with convergence code %d: %s",
                    optimisation$convergence, optimisation$message),
            call. = FALSE)
  }
  beta <- stats::setNames(optimisation$par[seq_len(p)],
                          names(object$coefficients))
  covariance <- object$covariance
  if (covariance_count) {
    covariance[covariance_names] <-
      exp(optimisation$par[p + seq_len(covariance_count)])
  }
  random_variance <- .sp_mcml_random_variance(optimisation$par, object)
  random_slope_specification <- .sp_mcml_random_slope(
    optimisation$par, object
  )
  temporal_scale <- .sp_mcml_temporal_scale(optimisation$par, object)
  fitted_intervention <- .sp_mcml_intervention(optimisation$par, object)
  fitted_multiplier <- .sp_intervention_multiplier(
    fitted_intervention, object$data_used
  )
  final_state <- .sp_laplace_state(
    beta, covariance, distance, object$x, object$offset, object$y,
    object$trials, object$family, object$smoothness,
    probability_multiplier = fitted_multiplier,
    temporal_lag = .sp_temporal_lag(object$temporal_values),
    temporal_scale = temporal_scale,
    random_z = if (is.null(object$random_specification)) NULL else
      object$random_specification$z,
    random_variance = random_variance,
    random_column_group = if (is.null(object$random_specification)) NULL else
      object$random_specification$column_group,
    random_slope_specification = random_slope_specification,
    initial_mode = object$latent_mode$mode
  )
  if (is.null(final_state)) {
    stop("The final MCML parameters did not produce a usable latent state.",
         call. = FALSE)
  }
  final_log_weights <- .sp_mcml_log_weights(
    optimisation$par, object, samples, log_proposal, distance
  )
  diagnostics <- .sp_mcml_diagnostics(final_log_weights, proposal_component)
  if (diagnostics$effective_sample_fraction < control$ess_warning) {
    warning(sprintf(
      "Importance-sampling ESS is %.1f%% of `nsim`; increase `nsim`, adjust `proposal_scale`, or refit the proposal.",
      100 * diagnostics$effective_sample_fraction), call. = FALSE)
  }

  parameter_vcov <- matrix(NA_real_, length(optimisation$par),
                           length(optimisation$par),
                           dimnames = list(names(parameter_start),
                                           names(parameter_start)))
  if (all(is.finite(optimisation$hessian))) {
    chol_hessian <- tryCatch(chol(optimisation$hessian),
                             error = function(e) NULL)
    if (!is.null(chol_hessian)) parameter_vcov <- chol2inv(chol_hessian)
  }
  coefficient_vcov <- parameter_vcov[seq_len(p), seq_len(p), drop = FALSE]
  dimnames(coefficient_vcov) <- list(names(beta), names(beta))
  covariance_vcov <- matrix(NA_real_, 3L, 3L,
    dimnames = list(names(covariance), names(covariance)))
  if (covariance_count && all(is.finite(parameter_vcov))) {
    log_indices <- p + seq_len(covariance_count)
    jacobian <- diag(unname(covariance[covariance_names]), covariance_count)
    covariance_vcov[covariance_names, covariance_names] <- jacobian %*%
      parameter_vcov[log_indices, log_indices, drop = FALSE] %*% jacobian
  }
  random_variance_vcov <- object$random_variance_vcov
  if (object$covariance_estimated && length(random_variance)) {
    random_indices <- p + covariance_count + seq_along(random_variance)
    random_variance_vcov <- matrix(
      NA_real_, length(random_variance), length(random_variance),
      dimnames = list(names(random_variance), names(random_variance))
    )
    if (all(is.finite(parameter_vcov))) {
      jacobian <- diag(unname(random_variance), length(random_variance))
      random_variance_vcov <- jacobian %*%
        parameter_vcov[random_indices, random_indices, drop = FALSE] %*%
        jacobian
      dimnames(random_variance_vcov) <- list(names(random_variance),
                                             names(random_variance))
    }
  }
  random_slope_values <- .sp_random_slope_parameters(
    random_slope_specification
  )
  random_slope_vcov <- matrix(
    NA_real_, length(random_slope_values), length(random_slope_values),
    dimnames = list(names(random_slope_values), names(random_slope_values))
  )
  if (random_slope_count && all(is.finite(parameter_vcov))) {
    slope_start <- p + covariance_count +
      (if (object$covariance_estimated) length(object$random_variance) else 0L)
    slope_indices <- slope_start + seq_len(random_slope_count)
    slope_jacobian <- .sp_random_slope_natural_jacobian(
      random_slope_specification
    )
    estimated_rows <- rowSums(abs(slope_jacobian)) > 0
    random_slope_vcov[estimated_rows, estimated_rows] <-
      slope_jacobian[estimated_rows, , drop = FALSE] %*%
      parameter_vcov[slope_indices, slope_indices, drop = FALSE] %*%
      t(slope_jacobian[estimated_rows, , drop = FALSE])
  }
  temporal_scale_se <- object$temporal_scale_se
  if (isTRUE(object$temporal_scale_estimated)) {
    temporal_index <- p + covariance_count +
      (if (object$covariance_estimated) length(object$random_variance) else 0L) +
      random_slope_count + 1L
    temporal_variance <- parameter_vcov[temporal_index, temporal_index]
    temporal_scale_se <- if (is.finite(temporal_variance) &&
      temporal_variance >= 0) temporal_scale * sqrt(temporal_variance) else
        NA_real_
  }
  intervention_vcov <- object$intervention_vcov
  if (length(intervention_indices)) {
    intervention_vcov <- matrix(NA_real_, 2L, 2L,
      dimnames = list(c("alpha", "decay"), c("alpha", "decay")))
    if (all(is.finite(parameter_vcov))) {
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
  }

  answer <- object
  answer$coefficients <- beta
  answer$covariance <- covariance
  answer$random_variance <- random_variance
  answer$random_variance_vcov <- random_variance_vcov
  answer$random_slope_specification <- random_slope_specification
  answer$random_slope <- .sp_random_slope_export(random_slope_specification)
  answer$random_slope_vcov <- random_slope_vcov
  answer$random_slope_parameterization <- .sp_random_slope_parameterization(
    random_slope_specification
  )
  answer$temporal_scale <- temporal_scale
  answer$temporal$scale <- temporal_scale
  answer$temporal_scale_se <- temporal_scale_se
  answer$intervention <- fitted_intervention
  answer$intervention_multiplier <- fitted_multiplier
  answer$intervention_vcov <- intervention_vcov
  answer$parameter_vcov <- parameter_vcov
  answer$coefficient_vcov <- coefficient_vcov
  answer$covariance_vcov <- covariance_vcov
  answer$fitted_state <- final_state
  answer$latent_mode <- final_state$mode
  answer$residual_mean <- object$y - final_state$mode$conditional_mean
  answer$log_likelihood <- .sp_log_mean_exp(final_log_weights)
  answer$likelihood_method <- "Monte Carlo ML"
  answer$optimisation <- optimisation
  answer$mcml_call <- match.call()
  answer$mcml <- c(list(
    control = control, seed = seed,
    reference_parameters = parameter_start,
    original_proposal_mean = original_proposal_mean,
    original_proposal_covariance = original_proposal_covariance,
    proposal_mean = proposal_mean,
    proposal_base_covariance = proposal_base_covariance,
    proposal_covariance = proposal_covariances[[1L]],
    proposal_covariances = proposal_covariances,
    proposal_component = proposal_component,
    proposal_probability = proposal_probability,
    log_proposal = log_proposal,
    latent_samples = samples, log_weights = final_log_weights,
    adaptation = if (is.null(pilot)) NULL else list(
      used = !is.null(adaptation), reason = adaptation_reason,
      shrinkage = control$adapt_shrinkage,
      pilot_samples = pilot$samples,
      pilot_log_proposal = pilot$log_density,
      pilot_log_weights = pilot_log_weights,
      pilot_component = pilot$component,
      pilot_probability = pilot$probability,
      pilot_diagnostics = pilot_diagnostics,
      weighted_mean = if (is.null(adaptation)) NULL else
        adaptation$weighted_mean,
      weighted_covariance = if (is.null(adaptation)) NULL else
        adaptation$weighted_covariance
    )
  ), diagnostics)
  answer
}

#' @export
print.spmodels_mcml_control <- function(x, ...) {
  cat("SPmodels MCML control\n")
  cat("Importance samples:", x$nsim, " Proposal scale:",
      paste(format(x$proposal_scale), collapse = ", "), "\n")
  if (isTRUE(x$adapt)) {
    cat("Adaptive pilot:", x$pilot_nsim, " Shrinkage:",
        format(x$adapt_shrinkage), "\n")
  }
  invisible(x)
}

.sp_mcml_seed_sequence <- function(seed, count) {
  if (is.null(seed)) return(rep(list(NULL), count))
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      !is.finite(seed) || seed != as.integer(seed)) {
    stop("`seed` must be NULL or one finite integer.", call. = FALSE)
  }
  modulus <- .Machine$integer.max
  values <- (as.double(as.integer(seed)) + seq_len(count) - 1) %% modulus
  values[values == 0] <- 1
  as.list(as.integer(values))
}

#' Repeat Monte Carlo maximum-likelihood fits independently
#'
#' Run [sp_mcml()] with independent importance samples and summarise the
#' resulting estimates and likelihood diagnostics. Replication diagnoses
#' Monte Carlo sensitivity; it does not replace scientific sampling
#' uncertainty or combine the fits into a new estimator.
#'
#' @param object A fitted binomial-logit or Poisson-log `spmodels_fit`.
#' @param nrep Number of independent MCML fits, at least two.
#' @param control Settings returned by [sp_mcml_control()].
#' @param seed Optional integer defining consecutive, call-local replicate
#'   seeds.
#' @param verbose Logical scalar passed to [sp_mcml()].
#'
#' @return An object of class `spmodels_mcml_replicates` containing every fit,
#'   replicate estimates, between-replicate summaries and MCML diagnostics.
#' @export
sp_mcml_replicates <- function(object, nrep = 3L,
                               control = sp_mcml_control(), seed = NULL,
                               verbose = FALSE) {
  if (!is.numeric(nrep) || length(nrep) != 1L || is.na(nrep) ||
      !is.finite(nrep) || nrep < 2 || nrep != as.integer(nrep)) {
    stop("`nrep` must be one integer not less than two.", call. = FALSE)
  }
  nrep <- as.integer(nrep)
  seeds <- .sp_mcml_seed_sequence(seed, nrep)
  fits <- vector("list", nrep)
  for (index in seq_len(nrep)) {
    fits[[index]] <- sp_mcml(object, control = control, seed = seeds[[index]],
                             verbose = verbose)
  }
  parameter_names <- names(fits[[1L]]$optimisation$par)
  estimates <- do.call(rbind, lapply(fits, function(fit) {
    stats::setNames(as.numeric(fit$optimisation$par), parameter_names)
  }))
  rownames(estimates) <- paste0("replicate_", seq_len(nrep))
  summary_table <- data.frame(
    term = parameter_names,
    mean = colMeans(estimates),
    between_sd = apply(estimates, 2L, stats::sd),
    mean_mc_std_error = apply(estimates, 2L, stats::sd) / sqrt(nrep),
    min = apply(estimates, 2L, min),
    max = apply(estimates, 2L, max),
    row.names = NULL,
    check.names = FALSE
  )
  mc_covariance <- stats::cov(estimates)
  dimnames(mc_covariance) <- list(parameter_names, parameter_names)
  parameter_covariances <- lapply(fits, `[[`, "parameter_vcov")
  usable_covariance <- vapply(parameter_covariances, function(value) {
    is.matrix(value) && identical(dim(value), c(length(parameter_names),
      length(parameter_names))) && all(is.finite(value))
  }, logical(1))
  sampling_covariance <- if (all(usable_covariance)) {
    Reduce(`+`, parameter_covariances) / nrep
  } else {
    matrix(NA_real_, length(parameter_names), length(parameter_names),
           dimnames = list(parameter_names, parameter_names))
  }
  dimnames(sampling_covariance) <- list(parameter_names, parameter_names)
  total_covariance <- if (all(is.finite(sampling_covariance))) {
    sampling_covariance + mc_covariance
  } else sampling_covariance
  summary_table$sampling_std_error <- sqrt(diag(sampling_covariance))
  summary_table$mc_std_error <- sqrt(diag(mc_covariance))
  summary_table$total_std_error <- sqrt(diag(total_covariance))
  natural_parameter <- function(fit) {
    slope_jacobian <- .sp_random_slope_natural_jacobian(
      fit$random_slope_specification
    )
    slope_values <- .sp_random_slope_parameters(
      fit$random_slope_specification
    )
    estimated_slope <- if (ncol(slope_jacobian))
      rowSums(abs(slope_jacobian)) > 0 else rep(FALSE, length(slope_values))
    c(stats::setNames(fit$coefficients,
                      paste0("mean::", names(fit$coefficients))),
      if (isTRUE(fit$covariance_estimated)) {
        covariance_names <- .sp_fitted_covariance_names(fit)
        stats::setNames(fit$covariance[covariance_names],
                        paste0("covariance::", covariance_names))
      } else numeric(),
      if (isTRUE(fit$covariance_estimated) && length(fit$random_variance))
        stats::setNames(fit$random_variance,
                        paste0("random::", names(fit$random_variance))) else
          numeric(),
      if (any(estimated_slope)) stats::setNames(
        slope_values[estimated_slope],
        paste0("random_slope::", names(slope_values)[estimated_slope])
      ) else numeric(),
      if (isTRUE(fit$temporal_scale_estimated))
        c(`temporal::scale` = fit$temporal_scale) else numeric(),
      if (isTRUE(fit$intervention_estimated)) c(
        `intervention::alpha` = fit$intervention$alpha,
        `intervention::decay` = fit$intervention$decay
      ) else numeric())
  }
  natural_estimates <- do.call(rbind, lapply(fits, natural_parameter))
  rownames(natural_estimates) <- rownames(estimates)
  natural_names <- colnames(natural_estimates)
  natural_covariances <- lapply(seq_along(fits), function(index) {
    fit <- fits[[index]]
    slope_jacobian <- .sp_random_slope_natural_jacobian(
      fit$random_slope_specification
    )
    estimated_slope <- if (ncol(slope_jacobian))
      rowSums(abs(slope_jacobian)) > 0 else logical(nrow(slope_jacobian))
    diagonal_derivative <- c(rep(1, length(fit$coefficients)),
      if (isTRUE(fit$covariance_estimated))
        unname(fit$covariance[.sp_fitted_covariance_names(fit)]) else
          numeric(),
      if (isTRUE(fit$covariance_estimated) && length(fit$random_variance))
        unname(fit$random_variance) else numeric())
    leading_count <- length(diagonal_derivative)
    intervention_derivative <-
      if (isTRUE(fit$intervention_estimated)) c(
        fit$intervention$alpha * (1 - fit$intervention$alpha),
        fit$intervention$decay
      ) else numeric()
    temporal_derivative <- if (isTRUE(fit$temporal_scale_estimated))
      fit$temporal_scale else numeric()
    natural_count <- leading_count + sum(estimated_slope) +
      length(temporal_derivative) + length(intervention_derivative)
    transformed_count <- leading_count + ncol(slope_jacobian) +
      length(temporal_derivative) + length(intervention_derivative)
    jacobian <- matrix(0, natural_count, transformed_count)
    if (leading_count) {
      jacobian[seq_len(leading_count), seq_len(leading_count)] <-
        diag(diagonal_derivative, leading_count)
    }
    if (ncol(slope_jacobian)) {
      natural_rows <- leading_count + seq_len(sum(estimated_slope))
      transformed_columns <- leading_count + seq_len(ncol(slope_jacobian))
      jacobian[natural_rows, transformed_columns] <-
        slope_jacobian[estimated_slope, , drop = FALSE]
    }
    if (length(temporal_derivative)) {
      natural_row <- leading_count + sum(estimated_slope) + 1L
      transformed_column <- leading_count + ncol(slope_jacobian) + 1L
      jacobian[natural_row, transformed_column] <- temporal_derivative
    }
    if (length(intervention_derivative)) {
      natural_rows <- natural_count - length(intervention_derivative) +
        seq_along(intervention_derivative)
      transformed_columns <- transformed_count -
        length(intervention_derivative) + seq_along(intervention_derivative)
      jacobian[natural_rows, transformed_columns] <-
        diag(intervention_derivative, length(intervention_derivative))
    }
    covariance <- fit$parameter_vcov
    if (!is.matrix(covariance) || any(dim(covariance) != transformed_count) ||
        any(!is.finite(covariance))) return(NULL)
    answer <- jacobian %*% covariance %*% jacobian
    dimnames(answer) <- list(natural_names, natural_names)
    answer
  })
  natural_sampling_covariance <- if (all(vapply(
    natural_covariances, is.matrix, logical(1)))) {
    Reduce(`+`, natural_covariances) / nrep
  } else matrix(NA_real_, length(natural_names), length(natural_names),
                dimnames = list(natural_names, natural_names))
  natural_mc_covariance <- stats::cov(natural_estimates)
  dimnames(natural_mc_covariance) <- list(natural_names, natural_names)
  natural_total_covariance <- if (all(is.finite(
    natural_sampling_covariance))) {
    natural_sampling_covariance + natural_mc_covariance
  } else natural_sampling_covariance
  natural_summary <- data.frame(
    term = natural_names,
    mean = colMeans(natural_estimates),
    sampling_std_error = sqrt(diag(natural_sampling_covariance)),
    mc_std_error = sqrt(diag(natural_mc_covariance)),
    total_std_error = sqrt(diag(natural_total_covariance)),
    mean_mc_std_error = sqrt(diag(natural_mc_covariance)) / sqrt(nrep),
    row.names = NULL, check.names = FALSE
  )
  diagnostics <- data.frame(
    replicate = seq_len(nrep),
    seed = if (is.null(seed)) rep(NA_integer_, nrep) else unlist(seeds),
    log_likelihood = vapply(fits, function(fit) fit$log_likelihood,
                            numeric(1L)),
    effective_sample_size = vapply(
      fits, function(fit) fit$mcml$effective_sample_size, numeric(1L)),
    effective_sample_fraction = vapply(
      fits, function(fit) fit$mcml$effective_sample_fraction, numeric(1L)),
    maximum_weight = vapply(
      fits, function(fit) fit$mcml$maximum_weight, numeric(1L)),
    relative_mcse = vapply(
      fits, function(fit) fit$mcml$relative_mcse, numeric(1L)),
    adaptation_used = vapply(fits, function(fit) {
      if (is.null(fit$mcml$adaptation)) NA else fit$mcml$adaptation$used
    }, logical(1L)),
    pilot_effective_sample_fraction = vapply(fits, function(fit) {
      if (is.null(fit$mcml$adaptation)) NA_real_ else
        fit$mcml$adaptation$pilot_diagnostics$effective_sample_fraction
    }, numeric(1L)),
    pilot_maximum_weight = vapply(fits, function(fit) {
      if (is.null(fit$mcml$adaptation)) NA_real_ else
        fit$mcml$adaptation$pilot_diagnostics$maximum_weight
    }, numeric(1L)),
    pilot_relative_mcse = vapply(fits, function(fit) {
      if (is.null(fit$mcml$adaptation)) NA_real_ else
        fit$mcml$adaptation$pilot_diagnostics$relative_mcse
    }, numeric(1L))
  )
  structure(list(
    call = match.call(), control = control, fits = fits,
    estimates = estimates, summary = summary_table,
    sampling_covariance = sampling_covariance,
    mc_covariance = mc_covariance,
    total_covariance = total_covariance,
    natural_estimates = natural_estimates,
    natural_summary = natural_summary,
    natural_sampling_covariance = natural_sampling_covariance,
    natural_mc_covariance = natural_mc_covariance,
    natural_total_covariance = natural_total_covariance,
    diagnostics = diagnostics
  ), class = "spmodels_mcml_replicates")
}

#' Escalate MCML sample size to declared precision targets
#'
#' Refit with geometrically increasing independent importance-sample sizes
#' until both the relative likelihood Monte Carlo standard error and the
#' effective-sample fraction meet their declared targets, or `max_nsim` is
#' reached. With `refresh = TRUE`, each unsuccessful fit supplies the latent
#' Gaussian proposal centre and curvature for the next attempt.
#'
#' @param object A fitted binomial-logit or Poisson-log `spmodels_fit`.
#' @param control Initial settings returned by [sp_mcml_control()].
#' @param target_relative_mcse Positive target for the relative Monte Carlo
#'   standard error of the likelihood estimate.
#' @param min_ess_fraction Required effective-sample fraction in `(0, 1)`.
#' @param max_nsim Largest permitted importance-sample count.
#' @param growth Multiplicative sample-size growth greater than one.
#' @param seed Optional integer defining consecutive, call-local attempt seeds.
#' @param refresh Whether to rebuild the next proposal from the preceding fit.
#' @param verbose Logical scalar passed to [sp_mcml()].
#'
#' @return The final `spmodels_fit`, with an `mcml_precision` component that
#'   records every attempted sample size, diagnostic and stopping reason.
#' @export
sp_mcml_precision <- function(object, control = sp_mcml_control(),
                              target_relative_mcse = 0.02,
                              min_ess_fraction = 0.1,
                              max_nsim = 32000L, growth = 2,
                              seed = NULL, refresh = TRUE,
                              verbose = FALSE) {
  if (!inherits(control, "spmodels_mcml_control")) {
    stop("`control` must be created by `sp_mcml_control()`.", call. = FALSE)
  }
  positive <- function(value, name) {
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        !is.finite(value) || value <= 0) {
      stop(sprintf("`%s` must be one finite positive number.", name),
           call. = FALSE)
    }
    as.numeric(value)
  }
  target_relative_mcse <- positive(target_relative_mcse,
                                   "target_relative_mcse")
  if (!is.numeric(min_ess_fraction) || length(min_ess_fraction) != 1L ||
      is.na(min_ess_fraction) || !is.finite(min_ess_fraction) ||
      min_ess_fraction <= 0 || min_ess_fraction >= 1) {
    stop("`min_ess_fraction` must be strictly between zero and one.",
         call. = FALSE)
  }
  if (!is.numeric(max_nsim) || length(max_nsim) != 1L || is.na(max_nsim) ||
      !is.finite(max_nsim) || max_nsim < control$nsim ||
      max_nsim != as.integer(max_nsim)) {
    stop("`max_nsim` must be an integer not less than `control$nsim`.",
         call. = FALSE)
  }
  growth <- positive(growth, "growth")
  if (growth <= 1) stop("`growth` must be greater than one.", call. = FALSE)
  if (!is.logical(refresh) || length(refresh) != 1L || is.na(refresh)) {
    stop("`refresh` must be TRUE or FALSE.", call. = FALSE)
  }
  possible_sizes <- control$nsim
  while (possible_sizes[[length(possible_sizes)]] < max_nsim) {
    current_size <- possible_sizes[[length(possible_sizes)]]
    next_size <- min(as.integer(max_nsim),
                     max(current_size + 1L,
                         ceiling(current_size * growth)))
    possible_sizes <- c(possible_sizes, next_size)
  }
  seeds <- .sp_mcml_seed_sequence(seed, length(possible_sizes))
  current_object <- object
  history <- vector("list", length(possible_sizes))
  achieved <- FALSE
  for (attempt in seq_along(possible_sizes)) {
    attempt_control <- control
    attempt_control$nsim <- possible_sizes[[attempt]]
    fit <- sp_mcml(current_object, control = attempt_control,
                   seed = seeds[[attempt]], verbose = verbose)
    achieved <- is.finite(fit$mcml$relative_mcse) &&
      fit$mcml$relative_mcse <= target_relative_mcse &&
      fit$mcml$effective_sample_fraction >= min_ess_fraction
    history[[attempt]] <- data.frame(
      attempt = attempt, nsim = possible_sizes[[attempt]],
      seed = if (is.null(seed)) NA_integer_ else seeds[[attempt]],
      relative_mcse = fit$mcml$relative_mcse,
      effective_sample_size = fit$mcml$effective_sample_size,
      effective_sample_fraction = fit$mcml$effective_sample_fraction,
      maximum_weight = fit$mcml$maximum_weight,
      adaptation_used = if (is.null(fit$mcml$adaptation)) NA else
        fit$mcml$adaptation$used,
      pilot_effective_sample_fraction =
        if (is.null(fit$mcml$adaptation)) NA_real_ else
          fit$mcml$adaptation$pilot_diagnostics$effective_sample_fraction,
      pilot_maximum_weight = if (is.null(fit$mcml$adaptation)) NA_real_ else
        fit$mcml$adaptation$pilot_diagnostics$maximum_weight,
      pilot_relative_mcse = if (is.null(fit$mcml$adaptation)) NA_real_ else
        fit$mcml$adaptation$pilot_diagnostics$relative_mcse,
      log_likelihood = fit$log_likelihood,
      targets_met = achieved
    )
    if (achieved || possible_sizes[[attempt]] == max_nsim) break
    if (refresh) current_object <- fit
  }
  history <- do.call(rbind, history[seq_len(attempt)])
  fit$mcml_precision <- list(
    target_relative_mcse = target_relative_mcse,
    min_ess_fraction = min_ess_fraction,
    max_nsim = as.integer(max_nsim), growth = growth, refresh = refresh,
    achieved = achieved,
    stopping_reason = if (achieved) "targets met" else "maximum nsim reached",
    history = history
  )
  if (!achieved) {
    warning("MCML precision targets were not met before `max_nsim`.",
            call. = FALSE)
  }
  fit
}

#' @export
print.spmodels_mcml_replicates <- function(x, ...) {
  cat("Independent SPmodels MCML replicates\n")
  cat("Replicates:", nrow(x$estimates), " Samples per replicate:",
      x$control$nsim, "\n")
  cat("Natural-scale parameter uncertainty:\n")
  print(x$natural_summary, row.names = FALSE)
  invisible(x)
}
