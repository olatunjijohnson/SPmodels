#' Control latent MCMC sampling
#'
#' @param iterations Total MALA iterations including burn-in.
#' @param burnin Number of initial iterations discarded.
#' @param thin Positive retention interval after burn-in.
#' @param step_size Positive MALA step size. The proposal is preconditioned by
#'   inverse curvature at the latent mode.
#' @param chains Positive number of independent chains.
#' @param initial_scale Non-negative multiplier for dispersed initial states
#'   drawn around the latent mode when `chains > 1` and `initial` is omitted.
#'
#' @return An object of class `spmodels_mcmc_control`.
#' @export
sp_mcmc_control <- function(iterations = 4000L, burnin = 1000L,
                            thin = 3L, step_size = 0.7, chains = 1L,
                            initial_scale = 2) {
  integer_scalar <- function(x, name, minimum) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
        !is.finite(x) || x < minimum || x != as.integer(x)) {
      stop(sprintf("`%s` must be one integer not less than %d.",
                   name, minimum), call. = FALSE)
    }
    as.integer(x)
  }
  iterations <- integer_scalar(iterations, "iterations", 1L)
  burnin <- integer_scalar(burnin, "burnin", 0L)
  thin <- integer_scalar(thin, "thin", 1L)
  chains <- integer_scalar(chains, "chains", 1L)
  if (burnin >= iterations) {
    stop("`burnin` must be smaller than `iterations`.", call. = FALSE)
  }
  if (!is.numeric(step_size) || length(step_size) != 1L ||
      is.na(step_size) || !is.finite(step_size) || step_size <= 0) {
    stop("`step_size` must be one finite positive number.", call. = FALSE)
  }
  if (!is.numeric(initial_scale) || length(initial_scale) != 1L ||
      is.na(initial_scale) || !is.finite(initial_scale) || initial_scale < 0) {
    stop("`initial_scale` must be one finite non-negative number.",
         call. = FALSE)
  }
  structure(list(iterations = iterations, burnin = burnin, thin = thin,
                 step_size = as.numeric(step_size), chains = chains,
                 initial_scale = as.numeric(initial_scale)),
            class = "spmodels_mcmc_control")
}

.sp_seed_scope <- function(seed) {
  if (is.null(seed)) return(NULL)
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      !is.finite(seed) || seed != as.integer(seed)) {
    stop("`seed` must be NULL or one finite integer.", call. = FALSE)
  }
  list(had_seed = exists(".Random.seed", envir = .GlobalEnv,
                         inherits = FALSE),
       old_seed = if (exists(".Random.seed", envir = .GlobalEnv,
                             inherits = FALSE))
         get(".Random.seed", envir = .GlobalEnv) else NULL,
       seed = as.integer(seed))
}

.sp_restore_seed <- function(scope) {
  if (is.null(scope)) return(invisible(NULL))
  if (scope$had_seed) {
    assign(".Random.seed", scope$old_seed, envir = .GlobalEnv)
  } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  invisible(NULL)
}

#' Sample a non-Gaussian latent spatial vector
#'
#' Run a Metropolis-adjusted Langevin algorithm (MALA) for a latent Gaussian
#' vector under binomial-logit or Poisson-log observations. The proposal is
#' preconditioned by inverse negative curvature at the latent mode.
#'
#' @inheritParams sp_latent_mode
#' @param control Settings from [sp_mcmc_control()].
#' @param seed Optional integer seed. A supplied seed is local to the call.
#' @param initial Optional initial latent vector, or a matrix with one row per
#'   chain and one column per latent component.
#'
#' @return An object of class `spmodels_latent_sample` containing retained
#'   samples, acceptance rate, mode and controls.
#' @export
sp_latent_sample <- function(y, linear_predictor, covariance,
                             family = stats::binomial(), trials = NULL,
                             initial = NULL,
                             control = sp_mcmc_control(), seed = NULL) {
  if (!inherits(control, "spmodels_mcmc_control")) {
    stop("`control` must be created by `sp_mcmc_control()`.", call. = FALSE)
  }
  mode_initial <- if (is.matrix(initial)) initial[1L, ] else initial
  mode <- sp_latent_mode(y, linear_predictor, covariance, family, trials,
                         initial = mode_initial)
  if (!mode$converged) stop("The latent mode did not converge.", call. = FALSE)
  family_name <- family$family
  if (family_name == "binomial" && is.null(trials)) {
    trials <- rep(1, length(y))
  }
  covariance_chol <- chol(covariance)
  prior_precision <- chol2inv(covariance_chol)
  preconditioner <- mode$covariance
  preconditioner_chol <- chol(preconditioner)
  preconditioner_precision <- chol2inv(preconditioner_chol)
  step_size <- control$step_size
  drift_factor <- step_size^2 / 2
  retain_at <- seq.int(control$burnin + control$thin,
                       control$iterations, by = control$thin)
  seed_scope <- .sp_seed_scope(seed)
  if (!is.null(seed_scope)) {
    on.exit(.sp_restore_seed(seed_scope), add = TRUE)
    set.seed(seed_scope$seed)
  }
  component_names <- paste0("s", seq_along(y))
  chain_names <- paste0("chain_", seq_len(control$chains))
  if (is.null(initial)) {
    initial_values <- matrix(mode$mode, nrow = control$chains,
                             ncol = length(y), byrow = TRUE)
    if (control$chains > 1L && control$initial_scale > 0) {
      initial_noise <- matrix(
        stats::rnorm((control$chains - 1L) * length(y)),
        nrow = control$chains - 1L
      )
      initial_values[-1L, ] <- initial_values[-1L, , drop = FALSE] +
        control$initial_scale * initial_noise %*% preconditioner_chol
    }
  } else if (is.matrix(initial)) {
    if (!is.numeric(initial) || anyNA(initial) || any(!is.finite(initial)) ||
        !identical(dim(initial), c(control$chains, length(y)))) {
      stop("Matrix `initial` must be finite with one row per chain and one column per latent component.",
           call. = FALSE)
    }
    initial_values <- initial
  } else {
    if (!is.numeric(initial) || length(initial) != length(y) ||
        anyNA(initial) || any(!is.finite(initial))) {
      stop("Vector `initial` must contain one finite value per latent component.",
           call. = FALSE)
    }
    initial_values <- matrix(initial, nrow = control$chains,
                             ncol = length(y), byrow = TRUE)
  }
  chain_samples <- array(
    NA_real_, dim = c(length(retain_at), length(y), control$chains),
    dimnames = list(NULL, component_names, chain_names)
  )
  accepted <- integer(control$chains)
  for (chain in seq_len(control$chains)) {
    current_s <- initial_values[chain, ]
    current <- .sp_latent_components(current_s, y, linear_predictor,
                                     prior_precision, family_name, trials)
    if (is.null(current) || !is.finite(current$log_posterior)) {
      stop(sprintf("Initial state for chain %d has a non-finite posterior.",
                   chain), call. = FALSE)
    }
    retained <- 0L
    for (iteration in seq_len(control$iterations)) {
      forward_mean <- current_s +
        drift_factor * as.vector(preconditioner %*% current$gradient)
      proposal_s <- forward_mean + step_size *
        as.vector(crossprod(preconditioner_chol,
                            stats::rnorm(length(y))))
      proposal <- .sp_latent_components(proposal_s, y, linear_predictor,
                                        prior_precision, family_name, trials)
      accept <- FALSE
      if (!is.null(proposal) && is.finite(proposal$log_posterior)) {
        reverse_mean <- proposal_s +
          drift_factor * as.vector(preconditioner %*% proposal$gradient)
        forward_delta <- proposal_s - forward_mean
        reverse_delta <- current_s - reverse_mean
        forward_quadratic <- sum(forward_delta *
          (preconditioner_precision %*% forward_delta)) / step_size^2
        reverse_quadratic <- sum(reverse_delta *
          (preconditioner_precision %*% reverse_delta)) / step_size^2
        log_ratio <- proposal$log_posterior - current$log_posterior -
          0.5 * reverse_quadratic + 0.5 * forward_quadratic
        accept <- log(stats::runif(1)) < min(0, log_ratio)
      }
      if (accept) {
        current_s <- proposal_s
        current <- proposal
        accepted[chain] <- accepted[chain] + 1L
      }
      if (iteration %in% retain_at) {
        retained <- retained + 1L
        chain_samples[retained, , chain] <- current_s
      }
    }
  }
  samples <- do.call(rbind, lapply(seq_len(control$chains), function(chain)
    matrix(chain_samples[, , chain], nrow = length(retain_at),
           ncol = length(y))))
  colnames(samples) <- component_names
  chain_id <- rep(chain_names, each = length(retain_at))
  acceptance_rates <- accepted / control$iterations
  structure(list(
    samples = samples, chain_samples = chain_samples, chain_id = chain_id,
    acceptance_rate = mean(acceptance_rates),
    acceptance_rates = stats::setNames(acceptance_rates, chain_names),
    accepted = stats::setNames(accepted, chain_names),
    iterations = control$iterations,
    retained_per_chain = length(retain_at), chains = control$chains,
    initial = initial_values,
    mode = mode, control = control, seed = seed,
    family = family_name, call = match.call()
  ), class = "spmodels_latent_sample")
}

.sp_effective_sample_size <- function(values, max_lag) {
  n <- length(values)
  autocorrelation <- as.numeric(stats::acf(
    values, plot = FALSE, lag.max = max_lag, demean = TRUE
  )$acf)[-1L]
  if (!length(autocorrelation)) return(c(ess = n, autocorrelation))
  usable <- length(autocorrelation) - length(autocorrelation) %% 2L
  if (usable < 2L) return(c(ess = n, autocorrelation))
  pair_sum <- autocorrelation[seq.int(1L, usable, by = 2L)] +
    autocorrelation[seq.int(2L, usable, by = 2L)]
  first_nonpositive <- which(pair_sum <= 0)[1L]
  pair_count <- if (is.na(first_nonpositive)) length(pair_sum) else
    first_nonpositive - 1L
  retained_rho <- if (pair_count > 0L) {
    autocorrelation[seq_len(2L * pair_count)]
  } else numeric()
  integrated_time <- max(1, 1 + 2 * sum(retained_rho))
  c(ess = min(n, n / integrated_time), autocorrelation)
}

.sp_split_chains <- function(chains) {
  n <- nrow(chains)
  half <- floor(n / 2L)
  if (half < 2L) return(NULL)
  cbind(chains[seq_len(half), , drop = FALSE],
        chains[seq.int(n - half + 1L, n), , drop = FALSE])
}

.sp_basic_rhat <- function(chains) {
  n <- nrow(chains)
  m <- ncol(chains)
  if (n < 2L || m < 2L) return(NA_real_)
  within <- apply(chains, 2L, stats::var)
  w <- mean(within)
  b <- n * stats::var(colMeans(chains))
  tolerance <- 100 * .Machine$double.eps *
    max(1, max(abs(chains)))^2
  if (w <= tolerance) return(if (b <= tolerance) 1 else Inf)
  variance_plus <- (n - 1) / n * w + b / n
  sqrt(max(0, variance_plus / w))
}

.sp_rank_normalize <- function(values) {
  ranks <- rank(values, ties.method = "average")
  stats::qnorm((ranks - 3 / 8) / (length(values) + 1 / 4))
}

.sp_rank_normalized_rhat <- function(chains) {
  split <- .sp_split_chains(chains)
  if (is.null(split)) return(NA_real_)
  rank_values <- matrix(.sp_rank_normalize(as.vector(split)),
                        nrow = nrow(split), ncol = ncol(split))
  folded <- abs(split - stats::median(as.vector(split)))
  folded_values <- matrix(.sp_rank_normalize(as.vector(folded)),
                          nrow = nrow(split), ncol = ncol(split))
  max(.sp_basic_rhat(rank_values), .sp_basic_rhat(folded_values))
}

.sp_multichain_ess <- function(chains, max_lag) {
  n <- nrow(chains)
  m <- ncol(chains)
  total <- n * m
  if (m == 1L) {
    return(unname(.sp_effective_sample_size(chains[, 1L], max_lag)[1L]))
  }
  chain_means <- colMeans(chains)
  w <- mean(apply(chains, 2L, stats::var))
  b <- n * stats::var(chain_means)
  variance_plus <- (n - 1) / n * w + b / n
  tolerance <- 100 * .Machine$double.eps *
    max(1, max(abs(chains)))^2
  if (!is.finite(variance_plus) || variance_plus <= tolerance) return(total)
  autocovariance <- vapply(seq_len(max_lag), function(lag) {
    mean(vapply(seq_len(m), function(chain) {
      centred <- chains[, chain] - chain_means[chain]
      sum(centred[seq_len(n - lag)] *
            centred[seq.int(lag + 1L, n)]) / n
    }, numeric(1)))
  }, numeric(1))
  rho <- 1 - (w - autocovariance) / variance_plus
  usable <- length(rho) - length(rho) %% 2L
  if (usable < 2L) return(total)
  pair_sum <- rho[seq.int(1L, usable, by = 2L)] +
    rho[seq.int(2L, usable, by = 2L)]
  positive_count <- which(pair_sum <= 0)[1L]
  positive_count <- if (is.na(positive_count)) length(pair_sum) else
    positive_count - 1L
  if (positive_count < 1L) return(total)
  pair_sum <- cummin(pair_sum[seq_len(positive_count)])
  integrated_time <- max(1, 1 + 2 * sum(pair_sum))
  min(total, total / integrated_time)
}

#' Diagnose retained latent samples
#'
#' Summarise component means, standard deviations and effective sample sizes,
#' and retain an autocorrelation sequence for a selected component.
#'
#' @param object An object returned by [sp_latent_sample()].
#' @param component Component index or name used for the trace and ACF display.
#' @param max_lag Positive maximum autocorrelation lag, smaller than the number
#'   of retained draws.
#' @param chains Optional chain indices or names to include. The default uses
#'   every chain.
#'
#' @return An object of class `spmodels_sampler_diagnostics`.
#' @export
sp_sampler_diagnostics <- function(object, component = 1L,
                                   max_lag = NULL, chains = NULL) {
  if (!inherits(object, "spmodels_latent_sample")) {
    stop("`object` must be returned by `sp_latent_sample()`.", call. = FALSE)
  }
  chain_samples <- object$chain_samples
  if (is.null(chain_samples)) {
    chain_samples <- array(object$samples,
      dim = c(nrow(object$samples), ncol(object$samples), 1L),
      dimnames = list(NULL, colnames(object$samples), "chain_1"))
  }
  chain_names <- dimnames(chain_samples)[[3L]]
  if (is.null(chains)) chains <- seq_along(chain_names)
  if (is.character(chains)) chains <- match(chains, chain_names)
  if (!is.numeric(chains) || !length(chains) || anyNA(chains) ||
      any(!is.finite(chains)) || any(chains != as.integer(chains)) ||
      any(chains < 1 | chains > length(chain_names)) || anyDuplicated(chains)) {
    stop("`chains` must identify distinct available chains.", call. = FALSE)
  }
  chains <- as.integer(chains)
  chain_samples <- chain_samples[, , chains, drop = FALSE]
  n <- dim(chain_samples)[1L]
  if (is.null(max_lag)) max_lag <- min(100L, n - 1L)
  if (!is.numeric(max_lag) || length(max_lag) != 1L || is.na(max_lag) ||
      max_lag < 1 || max_lag >= n || max_lag != as.integer(max_lag)) {
    stop("`max_lag` must be a positive integer smaller than retained sample size.",
         call. = FALSE)
  }
  component_names <- dimnames(chain_samples)[[2L]]
  if (is.character(component)) component <- match(component, component_names)
  if (!is.numeric(component) || length(component) != 1L || is.na(component) ||
      component < 1 || component > dim(chain_samples)[2L] ||
      component != as.integer(component)) {
    stop("`component` must identify one retained latent component.", call. = FALSE)
  }
  component_chains <- lapply(seq_len(dim(chain_samples)[2L]), function(j)
    matrix(chain_samples[, j, ], nrow = n, ncol = length(chains)))
  ess <- vapply(component_chains, .sp_multichain_ess, numeric(1),
                max_lag = as.integer(max_lag))
  rhat <- vapply(component_chains, .sp_rank_normalized_rhat, numeric(1))
  pooled <- lapply(component_chains, as.vector)
  standard_deviation <- vapply(pooled, stats::sd, numeric(1))
  selected_chains <- component_chains[[as.integer(component)]]
  selected_acf <- rowMeans(vapply(seq_len(ncol(selected_chains)), function(j) {
    as.numeric(stats::acf(selected_chains[, j], plot = FALSE,
      lag.max = as.integer(max_lag), demean = TRUE)$acf)[-1L]
  }, numeric(as.integer(max_lag))))
  summary <- data.frame(
    component = component_names,
    mean = vapply(pooled, mean, numeric(1)), sd = standard_deviation,
    ess = ess, mcse = standard_deviation / sqrt(ess), rhat = rhat,
    row.names = NULL
  )
  acceptance_rates <- object$acceptance_rates
  if (is.null(acceptance_rates)) acceptance_rates <- object$acceptance_rate
  acceptance_rates <- acceptance_rates[chains]
  trace <- selected_chains
  if (ncol(trace) == 1L) trace <- trace[, 1L]
  structure(list(
    summary = summary, acceptance_rate = mean(acceptance_rates),
    acceptance_rates = acceptance_rates,
    selected_chains = chain_names[chains],
    selected_component = component_names[as.integer(component)],
    trace = trace,
    lag = seq_along(selected_acf), autocorrelation = selected_acf,
    retained_draws = n * length(chains), retained_per_chain = n,
    chains = length(chains),
    converged = if (length(chains) > 1L)
      all(is.finite(rhat) & rhat <= 1.01) else NA,
    call = match.call()
  ), class = "spmodels_sampler_diagnostics")
}

#' @export
print.spmodels_sampler_diagnostics <- function(x, ...) {
  cat("SPmodels latent-sampler diagnostics\n")
  cat("Chains:", x$chains, " Retained per chain:",
      x$retained_per_chain, "\n")
  cat("Acceptance rates:",
      paste(format(x$acceptance_rates, digits = 3), collapse = ", "), "\n")
  if (x$chains > 1L) cat("All rank-normalized split R-hat <= 1.01:",
                         x$converged, "\n")
  print(x$summary, row.names = FALSE)
  invisible(x)
}

#' @export
plot.spmodels_sampler_diagnostics <- function(x, ...) {
  old_par <- graphics::par(mfrow = c(1, 2))
  on.exit(graphics::par(old_par), add = TRUE)
  if (is.matrix(x$trace)) {
    graphics::matplot(x$trace, type = "l", lty = 1,
                      xlab = "Retained draw", ylab = x$selected_component,
                      main = "Trace", ...)
  } else {
    graphics::plot(x$trace, type = "l", xlab = "Retained draw",
                   ylab = x$selected_component, main = "Trace", ...)
  }
  graphics::plot(x$lag, x$autocorrelation, type = "h",
                 xlab = "Lag", ylab = "Autocorrelation", main = "ACF")
  graphics::abline(h = 0, lty = 2)
  invisible(x)
}
