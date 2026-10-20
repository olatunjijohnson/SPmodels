if (file.exists(file.path("R", "nongaussian-fit.R"))) {
  for (file in list.files("R", full.names = TRUE)) source(file)
} else {
  library(SPmodels)
}

close <- function(a, b, tolerance = 1e-7) {
  isTRUE(all.equal(a, b, tolerance = tolerance, check.attributes = FALSE))
}
expect_error <- function(expression, pattern) {
  error <- tryCatch({ force(expression); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}

dat <- data.frame(
  successes = c(0, 1, 1, 2, 2, 3, 1, 3),
  trials = c(2, 2, 3, 3, 3, 4, 2, 4),
  exposure = c(-1, -0.6, -0.2, 0.1, 0.4, 0.8, 1.1, 1.5),
  east = c(0, 1, 0, 1, 2, 2, 3, 3),
  north = c(0, 0, 1, 1, 0, 1, 0, 1)
)
control <- sp_control(spatial_variance = 0.7, scale = 1.4,
                      nugget_variance = 0.2, estimate = FALSE,
                      distance_units = "km")
fit_counts <- sp_fit(
  cbind(successes, trials - successes) ~ exposure, dat,
  family = binomial(), coords = ~ east + north, control = control
)
stopifnot(fit_counts$family$family == "binomial",
          identical(fit_counts$likelihood_method, "Laplace"),
          fit_counts$latent_mode$converged,
          all(fitted(fit_counts) > 0),
          all(fitted(fit_counts) < dat$trials))
ordinary_binomial <- glm(
  cbind(successes, trials - successes) ~ exposure, data = dat,
  family = binomial()
)
ordinary_binomial_variogram <- sp_variogram(
  ordinary_binomial, coords = ~ east + north, data = dat, bins = 3
)
stopifnot(identical(ordinary_binomial_variogram$stage, "pre_fit"),
          identical(ordinary_binomial_variogram$family, "binomial"),
          identical(ordinary_binomial_variogram$residual_type, "response"),
          all(is.na(ordinary_binomial_variogram$values$fitted)))

# Independently reconstruct the Laplace log likelihood at the fitted values.
distance <- as.matrix(dist(as.matrix(dat[c("east", "north")])))
sigma <- 0.7 * matern_cor(distance, 1.4, 0.5) + diag(0.2, nrow(dat))
precision <- solve(sigma)
s <- fit_counts$latent_mode$mode
eta <- as.vector(model.matrix(~ exposure, dat) %*% coef(fit_counts) + s)
probability <- plogis(eta)
hessian <- precision + diag(dat$trials * probability * (1 - probability))
laplace_expected <- sum(dbinom(dat$successes, dat$trials, probability,
                               log = TRUE)) -
  0.5 * drop(crossprod(s, precision %*% s)) -
  0.5 * as.numeric(determinant(sigma, logarithm = TRUE)$modulus) -
  0.5 * as.numeric(determinant(hessian, logarithm = TRUE)$modulus)
stopifnot(close(as.numeric(logLik(fit_counts)), laplace_expected, 1e-9))

# The supported binomial encodings define the same likelihood.
dat$proportion <- dat$successes / dat$trials
fit_proportion <- sp_fit(
  proportion ~ exposure, dat, family = binomial(), weights = trials,
  coords = ~ east + north, control = control
)
stopifnot(close(coef(fit_counts), coef(fit_proportion), 1e-7),
          close(as.numeric(logLik(fit_counts)),
                as.numeric(logLik(fit_proportion)), 1e-7))

# A simultaneous row permutation changes neither estimates nor likelihood.
permutation <- c(8, 3, 6, 1, 5, 2, 7, 4)
fit_permuted <- sp_fit(
  cbind(successes, trials - successes) ~ exposure, dat[permutation, ],
  family = binomial(), coords = ~ east + north, control = control
)
stopifnot(close(coef(fit_counts), coef(fit_permuted), 1e-7),
          close(as.numeric(logLik(fit_counts)),
                as.numeric(logLik(fit_permuted)), 1e-7))

# Poisson counts accept exposure through the ordinary formula offset.
counts <- transform(dat, cases = c(0, 1, 1, 2, 2, 4, 3, 5),
                    person_time = c(1, 2, 1, 3, 2, 4, 3, 5))
fit_poisson <- sp_fit(
  cases ~ exposure + offset(log(person_time)), counts,
  family = poisson(), coords = ~ east + north, control = control
)
stopifnot(fit_poisson$latent_mode$converged,
          all(is.finite(coef(fit_poisson))), all(fitted(fit_poisson) > 0))

# Prediction integrates over the approximate Gaussian latent distribution.
new <- data.frame(exposure = c(-0.4, 0.9), east = c(0.5, 2.5),
                  north = c(0.4, 0.6), person_time = c(2, 3))
poisson_link <- predict(fit_poisson, new, type = "link", se.fit = TRUE,
                        target = "observation")
poisson_response <- predict(fit_poisson, new, type = "response",
                            se.fit = TRUE, target = "observation")
poisson_mean_expected <- exp(poisson_link$fit + poisson_link$se.fit^2 / 2)
poisson_variance_expected <- poisson_mean_expected +
  poisson_mean_expected^2 * expm1(poisson_link$se.fit^2)
stopifnot(close(poisson_response$fit, poisson_mean_expected, 1e-10),
          close(poisson_response$se.fit,
                sqrt(poisson_variance_expected), 1e-10))
poisson_score_observed <- c(0, 2)
poisson_score_data <- new[1:2, ]
poisson_score_link <- predict(fit_poisson, poisson_score_data, type = "link",
                              se.fit = TRUE, target = "observation")
poisson_score_expected <- vapply(1:2, function(i) -log(integrate(function(z) {
  dpois(poisson_score_observed[i],
        exp(poisson_score_link$fit[i] + poisson_score_link$se.fit[i] * z)) *
    dnorm(z)
}, -Inf, Inf, rel.tol = 1e-8)$value), numeric(1))
stopifnot(close(sp_predictive_log_score(
  fit_poisson, poisson_score_data, poisson_score_observed
), poisson_score_expected, 1e-9))

binomial_link <- predict(fit_counts, new, type = "link", se.fit = TRUE,
                         target = "observation")
binomial_response <- predict(fit_counts, new, type = "response",
                             se.fit = TRUE, target = "observation")
first_mean <- integrate(function(z) {
  plogis(binomial_link$fit[1] + binomial_link$se.fit[1] * z) * dnorm(z)
}, -Inf, Inf, rel.tol = 1e-10)$value
stopifnot(close(binomial_response$fit[1], first_mean, 1e-8),
          close(binomial_response$se.fit,
                sqrt(binomial_response$fit * (1 - binomial_response$fit)),
                1e-10))

# Fitted non-Gaussian simulation follows the declared latent Gaussian then
# response-family transform, and a supplied seed remains local.
set.seed(501)
seed_before <- .Random.seed
binomial_simulation <- simulate(fit_counts, nsim = 3, seed = 38)
stopifnot(identical(.Random.seed, seed_before))
set.seed(38)
simulation_z <- matrix(rnorm(nrow(dat) * 3), nrow = nrow(dat))
simulation_latent <- crossprod(fit_counts$fitted_state$chol_sigma,
                               simulation_z)
simulation_eta <- sweep(simulation_latent, 1,
  as.vector(fit_counts$x %*% coef(fit_counts) + fit_counts$offset), "+")
simulation_probability <- plogis(simulation_eta)
simulation_expected <- matrix(rbinom(
  length(simulation_probability), rep(fit_counts$trials, 3),
  as.vector(simulation_probability)
), nrow = nrow(dat))
stopifnot(close(as.matrix(binomial_simulation), simulation_expected, 0),
          all(as.matrix(binomial_simulation) >= 0),
          all(as.matrix(binomial_simulation) <= fit_counts$trials))
poisson_simulation <- simulate(fit_poisson, nsim = 2, seed = 39)
stopifnot(all(as.matrix(poisson_simulation) >= 0),
          all(as.matrix(poisson_simulation) ==
                floor(as.matrix(poisson_simulation))))
count_recovery <- sp_simulation_recovery(fit_counts, nsim = 1, seed = 40)
stopifnot(is.na(count_recovery$failures[1]),
          nrow(count_recovery$estimates) == 1L)

# Conditional residual diagnostics reproduce family definitions at the latent
# mode without claiming Gaussian decorrelation.
set.seed(606)
diagnostic_seed_before <- .Random.seed
binomial_diagnostics <- sp_diagnostics(fit_counts, seed = 77)
stopifnot(identical(.Random.seed, diagnostic_seed_before))
diagnostic_mean <- fitted(fit_counts)
diagnostic_response <- fit_counts$y - diagnostic_mean
diagnostic_probability <- diagnostic_mean / fit_counts$trials
diagnostic_variance <- fit_counts$trials * diagnostic_probability *
  (1 - diagnostic_probability)
diagnostic_deviance <- sign(diagnostic_response) * sqrt(
  binomial()$dev.resids(fit_counts$y / fit_counts$trials,
                        diagnostic_probability, fit_counts$trials)
)
set.seed(77)
diagnostic_uniform <- runif(nrow(dat),
  pbinom(fit_counts$y - 1, fit_counts$trials, diagnostic_probability),
  pbinom(fit_counts$y, fit_counts$trials, diagnostic_probability))
diagnostic_uniform <- pmin(1 - .Machine$double.eps,
                           pmax(.Machine$double.eps, diagnostic_uniform))
stopifnot(close(binomial_diagnostics$response_residual,
                diagnostic_response),
          close(binomial_diagnostics$pearson_residual,
                diagnostic_response / sqrt(diagnostic_variance)),
          close(binomial_diagnostics$deviance_residual,
                diagnostic_deviance),
          close(binomial_diagnostics$randomized_probability,
                diagnostic_uniform, 0),
          close(binomial_diagnostics$randomized_quantile_residual,
                qnorm(diagnostic_uniform), 0))
expect_error(sp_diagnostics(fit_counts, seed = 1.5), "finite integer")
set.seed(607)
variogram_seed_before <- .Random.seed
binomial_variogram <- sp_variogram(
  fit_counts, bins = 3, residual_type = "randomized_quantile", seed = 77
)
variogram_distance <- as.matrix(dist(fit_counts$coordinates))
variogram_triangle <- lower.tri(variogram_distance)
variogram_breaks <- seq(0, max(variogram_distance[variogram_triangle]),
                        length.out = 4)
variogram_group <- cut(variogram_distance[variogram_triangle],
  variogram_breaks, include.lowest = TRUE, right = TRUE)
variogram_pair <- 0.5 * outer(
  binomial_diagnostics$randomized_quantile_residual,
  binomial_diagnostics$randomized_quantile_residual, "-"
)[variogram_triangle]^2
variogram_expected <- as.numeric(tapply(
  variogram_pair, variogram_group, mean
))
stopifnot(identical(.Random.seed, variogram_seed_before),
          close(binomial_variogram$values$semivariance,
                variogram_expected),
          identical(binomial_variogram$stage, "post_fit"),
          all(is.na(binomial_variogram$values$fitted)),
          identical(binomial_variogram$family, "binomial"))
variogram_plot <- tempfile(fileext = ".pdf")
grDevices::pdf(variogram_plot)
plot(binomial_variogram)
grDevices::dev.off()
stopifnot(file.info(variogram_plot)$size > 0)
unlink(variogram_plot)
expect_error(sp_variogram(fit_counts, residual_type = "conditional"),
             "residual_type")
binomial_envelope <- sp_variogram_envelope(
  fit_counts, nsim = 19, method = "permutation", bins = 3,
  residual_type = "randomized_quantile", seed = 78
)
envelope_residual <- sp_diagnostics(
  fit_counts, seed = 78
)$randomized_quantile_residual
manual_envelope_curve <- function(value) {
  pair <- 0.5 * outer(value, value, "-")[variogram_triangle]^2
  as.numeric(tapply(pair, variogram_group, mean))
}
set.seed(78)
manual_envelope <- replicate(
  19, manual_envelope_curve(sample(envelope_residual))
)
stopifnot(close(binomial_envelope$simulations, manual_envelope),
          identical(binomial_envelope$family, "binomial"),
          binomial_envelope$p_value >= 1 / 20,
          binomial_envelope$p_value <= 1)
envelope_plot <- tempfile(fileext = ".pdf")
grDevices::pdf(envelope_plot)
plot(binomial_envelope)
grDevices::dev.off()
stopifnot(file.info(envelope_plot)$size > 0)
unlink(envelope_plot)
parametric_binomial <- suppressWarnings(sp_variogram_envelope(
  fit_counts, nsim = 9, method = "parametric", bins = 3,
  residual_type = "randomized_quantile", seed = 79
))
parametric_binomial_repeat <- suppressWarnings(sp_variogram_envelope(
  fit_counts, nsim = 9, method = "parametric", bins = 3,
  residual_type = "randomized_quantile", seed = 79
))
# Replay the response simulation, first complete refit and randomized residual
# draw independently; its binned curve must be the first envelope column.
set.seed(79)
parametric_responses <- simulate(fit_counts, nsim = 9)
parametric_data <- dat
parametric_data$successes <- parametric_responses[[1L]]
parametric_refit <- sp_fit(
  cbind(successes, trials - successes) ~ exposure, parametric_data,
  family = binomial(), coords = ~ east + north, control = control
)
parametric_residual <- sp_diagnostics(
  parametric_refit
)$randomized_quantile_residual
stopifnot(identical(parametric_binomial$simulations,
                    parametric_binomial_repeat$simulations),
          close(parametric_binomial$simulations[, 1L],
                manual_envelope_curve(parametric_residual)),
          identical(parametric_binomial$method, "parametric"),
          length(parametric_binomial$refit_convergence) == 9L,
          length(parametric_binomial$refit_warnings) == 9L,
          is.finite(parametric_binomial$p_value))
parametric_poisson <- suppressWarnings(sp_variogram_envelope(
  fit_poisson, nsim = 7, method = "parametric", bins = 3,
  residual_type = "pearson", seed = 80
))
parametric_poisson_repeat <- suppressWarnings(sp_variogram_envelope(
  fit_poisson, nsim = 7, method = "parametric", bins = 3,
  residual_type = "pearson", seed = 80
))
stopifnot(identical(parametric_poisson$simulations,
                    parametric_poisson_repeat$simulations),
          identical(parametric_poisson$family, "poisson"),
          length(parametric_poisson$refit_convergence) == 7L,
          is.finite(parametric_poisson$p_value))
expect_error(sp_fit(proportion ~ exposure, dat, family = binomial(),
                    weights = rep(2.5, nrow(dat)), coords = ~ east + north,
                    control = control), "integer trial")
expect_error(sp_fit(successes ~ exposure, dat, family = poisson(),
                    weights = trials, coords = ~ east + north,
                    control = control), "Poisson.*weights")

cat("Non-Gaussian Laplace fitting tests passed.\n")
