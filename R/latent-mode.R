.sp_latent_components <- function(s, y, linear_predictor, precision,
                                  family_name, trials,
                                  probability_multiplier = NULL) {
  eta <- linear_predictor + s
  if (family_name == "binomial") {
    baseline_probability <- stats::plogis(eta)
    if (is.null(probability_multiplier)) {
      probability_multiplier <- rep(1, length(y))
    }
    probability <- probability_multiplier * baseline_probability
    mean <- trials * probability
    if (all(probability_multiplier == 1)) {
      score <- y - mean
      weight <- trials * probability * (1 - probability)
    } else {
      q <- baseline_probability
      a <- probability_multiplier
      one_minus_probability <- 1 - probability
      q_derivative <- q * (1 - q)
      ratio_derivative <- a * q_derivative *
        (1 - 2 * q + a * q^2) / one_minus_probability^2
      score <- y * (1 - q) -
        (trials - y) * a * q_derivative / one_minus_probability
      weight <- y * q_derivative + (trials - y) * ratio_derivative
    }
    log_likelihood <- sum(stats::dbinom(y, size = trials,
                                        prob = probability, log = TRUE))
  } else {
    mean <- exp(eta)
    if (any(!is.finite(mean))) return(NULL)
    weight <- mean
    log_likelihood <- sum(stats::dpois(y, lambda = mean, log = TRUE))
  }
  precision_s <- as.vector(precision %*% s)
  list(
    eta = eta, mean = mean, weight = weight,
    gradient = if (family_name == "binomial") score - precision_s else
      y - mean - precision_s,
    negative_hessian = precision + diag(weight, length(weight)),
    log_posterior = log_likelihood - 0.5 * sum(s * precision_s)
  )
}

#' Calculate a latent Gaussian mode
#'
#' Find the posterior mode of a latent Gaussian vector for independent
#' binomial-logit or Poisson-log observations with a supplied prior covariance.
#' This deterministic operation is a building block for later Monte Carlo
#' maximum-likelihood fitting; it does not itself estimate model parameters.
#'
#' @param y Numeric outcome vector. Binomial values are success counts.
#' @param linear_predictor Fixed part of the linear predictor, including any
#'   offset, before adding the latent vector.
#' @param covariance Symmetric positive-definite latent covariance matrix.
#' @param family A binomial-logit or Poisson-log family object.
#' @param trials Binomial trial counts. The default is one trial per outcome.
#' @param probability_multiplier Optional known binomial prevalence multiplier
#'   in `(0, 1]`, such as that calculated from an intervention history. It is
#'   not used for Poisson outcomes.
#' @param initial Optional initial latent vector; defaults to zero.
#' @param maxit Positive maximum number of Newton iterations.
#' @param tolerance Positive convergence tolerance for the maximum step.
#'
#' @return A list containing the mode, local covariance, gradient, iteration
#'   count, convergence flag, fitted conditional means and log posterior.
#' @export
sp_latent_mode <- function(y, linear_predictor, covariance,
                           family = stats::binomial(), trials = NULL,
                           probability_multiplier = NULL,
                           initial = NULL, maxit = 100L,
                           tolerance = 1e-8) {
  if (!is.numeric(y) || !length(y) || anyNA(y) || any(!is.finite(y))) {
    stop("`y` must be a finite non-missing numeric vector.", call. = FALSE)
  }
  n <- length(y)
  if (!is.numeric(linear_predictor) || length(linear_predictor) != n ||
      anyNA(linear_predictor) || any(!is.finite(linear_predictor))) {
    stop("`linear_predictor` must be finite and have the same length as `y`.",
         call. = FALSE)
  }
  if (!is.matrix(covariance) || !is.numeric(covariance) ||
      !identical(dim(covariance), c(n, n)) || anyNA(covariance) ||
      any(!is.finite(covariance)) || !isSymmetric(covariance,
                                                   tol = 1e-12)) {
    stop("`covariance` must be a finite symmetric numeric matrix matching `y`.",
         call. = FALSE)
  }
  covariance_chol <- tryCatch(chol(covariance), error = function(e) NULL)
  if (is.null(covariance_chol)) {
    stop("`covariance` must be positive definite.", call. = FALSE)
  }
  if (!inherits(family, "family") ||
      !((family$family == "binomial" && family$link == "logit") ||
        (family$family == "poisson" && family$link == "log"))) {
    stop("Only binomial-logit and Poisson-log latent modes are implemented.",
         call. = FALSE)
  }
  family_name <- family$family
  if (family_name == "binomial") {
    if (is.null(trials)) trials <- rep(1, n)
    if (!is.numeric(trials) || length(trials) != n || anyNA(trials) ||
        any(!is.finite(trials)) || any(trials < 0) ||
        any(trials != floor(trials)) || any(y < 0) || any(y != floor(y)) ||
        any(y > trials)) {
      stop("Binomial `y` and `trials` must be integer counts with `0 <= y <= trials`.",
           call. = FALSE)
    }
    if (is.null(probability_multiplier)) {
      probability_multiplier <- rep(1, n)
    }
    if (!is.numeric(probability_multiplier) ||
        length(probability_multiplier) != n || anyNA(probability_multiplier) ||
        any(!is.finite(probability_multiplier)) ||
        any(probability_multiplier <= 0) ||
        any(probability_multiplier > 1)) {
      stop("`probability_multiplier` must contain one finite value in (0, 1] per observation.",
           call. = FALSE)
    }
  } else {
    if (!is.null(trials)) stop("`trials` is only used for binomial outcomes.", call. = FALSE)
    if (!is.null(probability_multiplier)) {
      stop("`probability_multiplier` is only used for binomial outcomes.",
           call. = FALSE)
    }
    if (any(y < 0) || any(y != floor(y))) {
      stop("Poisson `y` must contain non-negative integer counts.", call. = FALSE)
    }
  }
  if (is.null(initial)) initial <- numeric(n)
  if (!is.numeric(initial) || length(initial) != n || anyNA(initial) ||
      any(!is.finite(initial))) {
    stop("`initial` must be NULL or a finite numeric vector matching `y`.",
         call. = FALSE)
  }
  if (!is.numeric(maxit) || length(maxit) != 1L || is.na(maxit) ||
      maxit < 1 || maxit != as.integer(maxit)) {
    stop("`maxit` must be one positive integer.", call. = FALSE)
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      is.na(tolerance) || !is.finite(tolerance) || tolerance <= 0) {
    stop("`tolerance` must be one finite positive number.", call. = FALSE)
  }

  precision <- chol2inv(covariance_chol)
  s <- as.vector(initial)
  converged <- FALSE
  iteration <- 0L
  current <- .sp_latent_components(s, y, linear_predictor, precision,
                                   family_name, trials,
                                   probability_multiplier)
  for (iteration in seq_len(as.integer(maxit))) {
    chol_hessian <- tryCatch(chol(current$negative_hessian),
                             error = function(e) NULL)
    if (is.null(chol_hessian)) stop("Latent-mode curvature is not positive definite.", call. = FALSE)
    step <- .sp_solve_chol(chol_hessian, current$gradient)
    multiplier <- 1
    accepted <- FALSE
    while (multiplier >= 2^-20) {
      proposal_s <- s + multiplier * step
      proposal <- .sp_latent_components(proposal_s, y, linear_predictor,
                                        precision, family_name, trials,
                                        probability_multiplier)
      if (!is.null(proposal) && is.finite(proposal$log_posterior) &&
          proposal$log_posterior >= current$log_posterior) {
        accepted <- TRUE
        break
      }
      multiplier <- multiplier / 2
    }
    if (!accepted) break
    s <- proposal_s
    current <- proposal
    if (max(abs(multiplier * step)) <= tolerance * (1 + max(abs(s)))) {
      converged <- TRUE
      break
    }
  }
  final_chol <- chol(current$negative_hessian)
  structure(list(
    mode = s, covariance = chol2inv(final_chol),
    gradient = current$gradient, iterations = iteration,
    converged = converged, conditional_mean = current$mean,
    linear_predictor = current$eta,
    log_posterior = current$log_posterior,
    family = family_name,
    probability_multiplier = probability_multiplier
  ), class = "spmodels_latent_mode")
}
