source_testing <- file.exists(file.path("R", "mcml.R"))
if (source_testing) {
  for (file in list.files("R", full.names = TRUE)) source(file)
} else {
  library(SPmodels)
  .sp_mcml_log_weights <- getFromNamespace(".sp_mcml_log_weights", "SPmodels")
  .sp_log_mean_exp <- getFromNamespace(".sp_log_mean_exp", "SPmodels")
}

close <- function(a, b, tolerance = 1e-8) {
  isTRUE(all.equal(a, b, tolerance = tolerance, check.attributes = FALSE))
}
expect_error <- function(expression, pattern) {
  error <- tryCatch({ force(expression); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}

# In one dimension, importance sampling agrees with independent numerical
# integration of the marginal likelihood.
beta <- c("(Intercept)" = -0.35)
variance <- c(spatial_variance = 0.8, scale = 1,
              nugget_variance = 0.2)
scalar_object <- list(
  x = matrix(1, 1, 1), offset = 0, y = 1, trials = 2,
  family = binomial(), covariance = variance,
  covariance_estimated = FALSE, smoothness = 0.5
)
scalar_mode <- sp_latent_mode(1, beta, matrix(1, 1, 1),
                              family = binomial(), trials = 2)
set.seed(71)
standard <- matrix(rnorm(50000), ncol = 1)
root <- chol(scalar_mode$covariance)
samples <- sweep(standard %*% root, 2, scalar_mode$mode, "+")
log_q <- -0.5 * (log(2 * pi) + 2 * log(root[1, 1]) + standard[, 1]^2)
mc_weights <- .sp_mcml_log_weights(beta, scalar_object, samples, log_q,
                                    matrix(0, 1, 1))
mc_log_likelihood <- .sp_log_mean_exp(mc_weights)
exact_likelihood <- integrate(function(s) {
  dbinom(1, 2, plogis(beta + s)) * dnorm(s, sd = 1)
}, -Inf, Inf, rel.tol = 1e-11)$value
stopifnot(abs(mc_log_likelihood - log(exact_likelihood)) < 0.01)

# A deterministic two-scale mixture estimates the same scalar integral when
# every draw is divided by the complete normalized mixture density.
set.seed(72)
mixture_scales <- c(.6, 1.8)
mixture_component <- rep(1:2, length.out = 50000L)
mixture_standard <- rnorm(length(mixture_component))
mixture_samples <- matrix(
  scalar_mode$mode + mixture_standard * as.vector(root) *
    mixture_scales[mixture_component], ncol = 1L
)
mixture_probability <- tabulate(mixture_component) /
  length(mixture_component)
mixture_q <- vapply(mixture_scales, function(scale) {
  dnorm(mixture_samples[, 1L], scalar_mode$mode,
        as.vector(root) * scale)
}, numeric(nrow(mixture_samples))) %*% mixture_probability
mixture_weights <- .sp_mcml_log_weights(
  beta, scalar_object, mixture_samples, log(drop(mixture_q)), matrix(0, 1, 1)
)
stopifnot(abs(.sp_log_mean_exp(mixture_weights) -
  log(exact_likelihood)) < .02)

dat <- data.frame(
  y = c(0, 1, 0, 1, 1, 0, 1, 1),
  count = c(0, 1, 0, 2, 1, 1, 2, 3),
  x = seq(-1, 1, length.out = 8),
  east = rep(0:3, 2), north = rep(0:1, each = 4)
)
fit_control <- sp_control(
  spatial_variance = 0.7, scale = 1.2, nugget_variance = 0.2,
  estimate = FALSE, distance_units = "km"
)
laplace_fit <- sp_fit(y ~ x, dat, family = binomial(),
                      coords = ~ east + north, control = fit_control)
mc_control <- sp_mcml_control(nsim = 600, maxit = 60,
                              proposal_scale = 1.1)

# A supplied seed is local, creates a deterministic common-random-number
# objective, and reproduces all fitted values and importance diagnostics.
set.seed(913)
seed_before <- .Random.seed
mc_fit_one <- sp_mcml(laplace_fit, control = mc_control, seed = 44)
stopifnot(identical(.Random.seed, seed_before))
mc_fit_two <- sp_mcml(laplace_fit, control = mc_control, seed = 44)
stopifnot(close(coef(mc_fit_one), coef(mc_fit_two), 0),
          close(as.numeric(logLik(mc_fit_one)),
                as.numeric(logLik(mc_fit_two)), 0),
          close(mc_fit_one$mcml$log_weights,
                mc_fit_two$mcml$log_weights, 0),
          identical(mc_fit_one$likelihood_method, "Monte Carlo ML"),
          mc_fit_one$optimisation$convergence == 0,
          mc_fit_one$mcml$effective_sample_size > 0,
          mc_fit_one$mcml$effective_sample_size <= mc_control$nsim,
          mc_fit_one$mcml$maximum_weight > 0,
          mc_fit_one$mcml$maximum_weight <= 1,
          is.finite(mc_fit_one$mcml$relative_mcse),
          close(as.numeric(logLik(mc_fit_one)),
                .sp_log_mean_exp(mc_fit_one$mcml$log_weights), 1e-12))
mc_summary <- summary(mc_fit_one)
stopifnot(close(mc_summary$mcml_diagnostics$effective_sample_size,
                mc_fit_one$mcml$effective_sample_size, 0))

# Multiple proposal scales use balanced deterministic allocation and the exact
# normalized mixture density, not the generating-component density alone.
mixture_control <- sp_mcml_control(
  nsim = 101L, maxit = 20L, proposal_scale = c(.75, 1.5),
  ess_warning = .001
)
mixture_fit <- suppressWarnings(sp_mcml(
  laplace_fit, control = mixture_control, seed = 45L
))
component <- mixture_fit$mcml$proposal_component
stopifnot(identical(tabulate(component), c(51L, 50L)),
          close(mixture_fit$mcml$proposal_probability, c(51, 50) / 101, 0),
          length(mixture_fit$mcml$proposal_covariances) == 2L)
centered <- sweep(mixture_fit$mcml$latent_samples, 2L,
                  mixture_fit$mcml$proposal_mean, "-")
base_covariance <- laplace_fit$latent_mode$covariance
base_inverse <- solve(base_covariance)
base_log_determinant <- as.numeric(determinant(
  base_covariance, logarithm = TRUE
)$modulus)
quadratic <- rowSums((centered %*% base_inverse) * centered)
n_latent <- ncol(centered)
component_density <- vapply(mixture_control$proposal_scale, function(scale) {
  exp(-.5 * (n_latent * log(2 * pi) + base_log_determinant +
    2 * n_latent * log(scale) + quadratic / scale^2))
}, numeric(nrow(centered)))
direct_log_mixture <- log(drop(component_density %*%
  mixture_fit$mcml$proposal_probability))
stopifnot(close(mixture_fit$mcml$log_proposal, direct_log_mixture, 1e-11),
          all(is.finite(mixture_fit$mcml$log_weights)))

# Optional within-attempt adaptation uses a discarded pilot sample, shrinks
# its importance-weighted moments toward the Laplace proposal and then draws
# an independent retained sample from the fully normalized adapted mixture.
adaptive_control <- sp_mcml_control(
  nsim = 120L, maxit = 20L, proposal_scale = c(.8, 1.25),
  ess_warning = .001, adapt = TRUE, pilot_nsim = 110L,
  adapt_shrinkage = .6
)
adaptive_one <- suppressWarnings(sp_mcml(
  laplace_fit, control = adaptive_control, seed = 46L
))
adaptive_two <- suppressWarnings(sp_mcml(
  laplace_fit, control = adaptive_control, seed = 46L
))
adaptation <- adaptive_one$mcml$adaptation
stopifnot(!is.null(adaptation),
          isTRUE(adaptation$used),
          identical(adaptation$reason, "pilot ESS threshold met"),
          identical(dim(adaptation$pilot_samples), c(110L, 8L)),
          identical(dim(adaptive_one$mcml$latent_samples), c(120L, 8L)),
          all(is.finite(adaptation$pilot_log_weights)),
          !identical(adaptation$pilot_samples[seq_len(110L), ],
                     adaptive_one$mcml$latent_samples[seq_len(110L), ]),
          !is.null(chol(adaptive_one$mcml$proposal_base_covariance)),
          close(adaptive_one$mcml$latent_samples,
                adaptive_two$mcml$latent_samples, 0),
          close(adaptation$pilot_samples,
                adaptive_two$mcml$adaptation$pilot_samples, 0),
          close(adaptive_one$mcml$log_weights,
                adaptive_two$mcml$log_weights, 0))
adaptive_centered <- sweep(adaptive_one$mcml$latent_samples, 2L,
                           adaptive_one$mcml$proposal_mean, "-")
adaptive_inverse <- solve(adaptive_one$mcml$proposal_base_covariance)
adaptive_log_determinant <- as.numeric(determinant(
  adaptive_one$mcml$proposal_base_covariance, logarithm = TRUE
)$modulus)
adaptive_quadratic <- rowSums(
  (adaptive_centered %*% adaptive_inverse) * adaptive_centered
)
adaptive_component_density <- vapply(
  adaptive_control$proposal_scale, function(scale) {
    exp(-.5 * (n_latent * log(2 * pi) + adaptive_log_determinant +
      2 * n_latent * log(scale) + adaptive_quadratic / scale^2))
  }, numeric(nrow(adaptive_centered))
)
adaptive_direct_log_mixture <- log(drop(adaptive_component_density %*%
  adaptive_one$mcml$proposal_probability))
stopifnot(close(adaptive_one$mcml$log_proposal,
                adaptive_direct_log_mixture, 1e-11))
adaptive_summary <- summary(adaptive_one)
stopifnot(close(
  adaptive_summary$mcml_pilot_diagnostics$effective_sample_size,
  adaptation$pilot_diagnostics$effective_sample_size, 0
))
adaptive_print <- capture.output(print(adaptive_summary))
stopifnot(any(grepl("Adaptive pilot ESS", adaptive_print, fixed = TRUE)),
          any(grepl("Adaptive moments: used", adaptive_print, fixed = TRUE)))

# Full shrinkage is a useful boundary check: the independent pilot is still
# recorded and discarded, while the retained proposal moments equal Laplace.
full_shrinkage <- suppressWarnings(sp_mcml(
  laplace_fit,
  control = sp_mcml_control(nsim = 100L, maxit = 10L, adapt = TRUE,
    pilot_nsim = 100L, adapt_shrinkage = 1, ess_warning = .001),
  seed = 47L
))
stopifnot(close(full_shrinkage$mcml$proposal_mean,
                full_shrinkage$mcml$original_proposal_mean, 0),
          close(full_shrinkage$mcml$proposal_base_covariance,
                full_shrinkage$mcml$original_proposal_covariance, 0))

# A pilot below the declared ESS threshold is retained for audit but does not
# displace the original proposal moments. Lowering the threshold is the
# explicit opt-in when a user nevertheless wants to trust those moments.
fallback <- suppressWarnings(sp_mcml(
  laplace_fit,
  control = sp_mcml_control(nsim = 100L, maxit = 10L, adapt = TRUE,
    pilot_nsim = 100L, adapt_shrinkage = .5, ess_warning = .999999),
  seed = 48L
))
stopifnot(!fallback$mcml$adaptation$used,
          identical(fallback$mcml$adaptation$reason,
                    "pilot ESS fraction below ess_warning"),
          is.null(fallback$mcml$adaptation$weighted_mean),
          close(fallback$mcml$proposal_mean,
                fallback$mcml$original_proposal_mean, 0),
          close(fallback$mcml$proposal_base_covariance,
                fallback$mcml$original_proposal_covariance, 0),
          identical(summary(fallback)$mcml_adaptation_used, FALSE))

# Cross-validation repeats the declared MCML correction inside each training
# fold instead of reverting silently to the preliminary Laplace objective.
small_mc <- sp_mcml(laplace_fit,
  control = sp_mcml_control(nsim = 100, maxit = 20), seed = 6)
mc_cv <- sp_cross_validate(small_mc, folds = rep(1:2, each = 4))
stopifnot(nrow(mc_cv$predictions) == nrow(dat),
          all(is.finite(mc_cv$predictions$predicted)))

# The same likelihood machinery handles Poisson-log observations and retains
# a usable latent state for downstream prediction.
poisson_fit <- sp_fit(count ~ x, dat, family = poisson(),
                      coords = ~ east + north, control = fit_control)
poisson_mc <- sp_mcml(poisson_fit,
  control = sp_mcml_control(nsim = 400, maxit = 50), seed = 12)
poisson_prediction <- predict(poisson_mc, dat, type = "response")
stopifnot(all(is.finite(poisson_mc$mcml$log_weights)),
          all(is.finite(poisson_prediction)),
          all(poisson_prediction > 0))

expect_error(sp_mcml(sp_fit(y ~ x, dat, coords = ~ east + north,
  control = fit_control)), "binomial or Poisson")
expect_error(sp_mcml_control(nsim = 99), "not less than 100")
expect_error(sp_mcml_control(ess_warning = 1), "strictly between")
expect_error(sp_mcml_control(proposal_scale = c(1, 1)), "distinct")
expect_error(sp_mcml_control(adapt = NA), "TRUE or FALSE")
expect_error(sp_mcml_control(pilot_nsim = 99), "not less than 100")
expect_error(sp_mcml_control(adapt_shrinkage = 0), "greater than zero")
expect_error(sp_mcml_control(adapt_shrinkage = 1.01), "not greater")

# Independent repeats expose Monte Carlo variation without changing the input
# fit or the caller's random-number state.
repeat_control <- sp_mcml_control(nsim = 100, maxit = 20,
                                  ess_warning = 0.01)
set.seed(914)
repeat_seed_before <- .Random.seed
repeats <- suppressWarnings(sp_mcml_replicates(
  laplace_fit, nrep = 2, control = repeat_control, seed = 80
))
stopifnot(identical(.Random.seed, repeat_seed_before),
          identical(dim(repeats$estimates), c(2L, 2L)),
          identical(repeats$diagnostics$seed, 80:81),
          all(c("term", "mean", "between_sd", "min", "max") %in%
                names(repeats$summary)),
          close(repeats$summary$mean, colMeans(repeats$estimates), 0),
          close(repeats$summary$between_sd,
                apply(repeats$estimates, 2, sd), 0),
          close(repeats$summary$mean_mc_std_error,
                apply(repeats$estimates, 2, sd) / sqrt(2), 0),
          close(repeats$mc_covariance, stats::cov(repeats$estimates), 0),
          close(repeats$summary$mc_std_error,
                sqrt(diag(repeats$mc_covariance)), 0),
          close(repeats$natural_estimates, repeats$estimates, 0),
          close(repeats$natural_mc_covariance,
                stats::cov(repeats$natural_estimates), 0),
          close(repeats$natural_summary$mc_std_error,
                sqrt(diag(repeats$natural_mc_covariance)), 0),
          all(is.na(repeats$diagnostics$adaptation_used)),
          all(is.na(repeats$diagnostics$pilot_relative_mcse)))

# Replicate and precision histories expose pilot quality and the ESS-gated
# decision directly, while retaining the complete fits as before.
adaptive_history_control <- sp_mcml_control(
  nsim = 100L, maxit = 5L, adapt = TRUE, pilot_nsim = 100L,
  ess_warning = .999999
)
adaptive_repeats <- suppressWarnings(sp_mcml_replicates(
  laplace_fit, nrep = 2L, control = adaptive_history_control, seed = 82L
))
stopifnot(identical(adaptive_repeats$diagnostics$adaptation_used,
                    c(FALSE, FALSE)),
          all(is.finite(
            adaptive_repeats$diagnostics$pilot_effective_sample_fraction
          )),
          all(is.finite(adaptive_repeats$diagnostics$pilot_maximum_weight)),
          all(is.finite(adaptive_repeats$diagnostics$pilot_relative_mcse)))

# Natural-scale propagation exponentiates fitted log covariance parameters;
# unusable Hessian curvature stays missing while empirical MC covariance does
# not disappear.
estimated_laplace <- sp_fit(
  y ~ x, dat, family = binomial(), coords = ~ east + north,
  control = sp_control(distance_units = "km", maxit = 80L)
)
estimated_repeats <- suppressWarnings(sp_mcml_replicates(
  estimated_laplace, nrep = 2L,
  control = sp_mcml_control(nsim = 100L, maxit = 10L,
                            ess_warning = .001),
  seed = 121L
))
expected_natural <- do.call(rbind, lapply(estimated_repeats$fits, function(fit) {
  c(stats::setNames(coef(fit), paste0("mean::", names(coef(fit)))),
    stats::setNames(fit$covariance,
                    paste0("covariance::", names(fit$covariance))))
}))
stopifnot(close(estimated_repeats$natural_estimates, expected_natural, 0),
          close(estimated_repeats$natural_mc_covariance,
                stats::cov(expected_natural), 0),
          all(is.na(estimated_repeats$natural_sampling_covariance)),
          all(is.na(estimated_repeats$natural_total_covariance)))

# A permissive precision target stops at the first sample size and records an
# auditable history. An impossible target reaches the declared cap exactly.
precision <- suppressWarnings(sp_mcml_precision(
  laplace_fit, control = repeat_control, target_relative_mcse = 10,
  min_ess_fraction = 1e-6, max_nsim = 200, seed = 90
))
stopifnot(precision$mcml_precision$achieved,
          precision$mcml_precision$history$nsim == 100L,
          precision$mcml_precision$history$targets_met,
          is.na(precision$mcml_precision$history$adaptation_used))
adaptive_precision <- suppressWarnings(sp_mcml_precision(
  laplace_fit, control = adaptive_history_control,
  target_relative_mcse = 10, min_ess_fraction = 1e-6,
  max_nsim = 100L, seed = 92L
))
stopifnot(!adaptive_precision$mcml_precision$history$adaptation_used,
          is.finite(adaptive_precision$mcml_precision$history$
            pilot_effective_sample_fraction),
          is.finite(adaptive_precision$mcml_precision$history$
            pilot_relative_mcse))
capped <- suppressWarnings(sp_mcml_precision(
  laplace_fit, control = repeat_control, target_relative_mcse = 1e-12,
  min_ess_fraction = 0.999999, max_nsim = 200, growth = 2,
  seed = 91, refresh = FALSE
))
stopifnot(!capped$mcml_precision$achieved,
          identical(capped$mcml_precision$history$nsim, c(100, 200)),
          identical(capped$mcml_precision$stopping_reason,
                    "maximum nsim reached"))
expect_error(sp_mcml_replicates(laplace_fit, nrep = 1), "not less than two")
expect_error(sp_mcml_precision(laplace_fit, control = repeat_control,
  max_nsim = 99), "not less than")
expect_error(sp_mcml_precision(laplace_fit, control = repeat_control,
  growth = 1), "greater than one")

cat("Monte Carlo maximum-likelihood tests passed.\n")
