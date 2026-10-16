if (file.exists(file.path("R", "spacetime.R"))) {
  for (file in list.files("R", full.names = TRUE)) source(file)
} else {
  library(SPmodels)
}

close <- function(a, b, tolerance = 1e-9) {
  isTRUE(all.equal(a, b, tolerance = tolerance, check.attributes = FALSE))
}
expect_error <- function(expression, pattern) {
  error <- tryCatch({ force(expression); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}

dat <- data.frame(
  y = c(1.2, 0.7, 2.1, 1.8, 3.0, 2.6, 2.2, 3.4),
  x = c(-1, 0, 0.5, 1, 1.5, 2, 0.8, 2.4),
  east = c(0, 1, 0, 2, 1, 3, 2.5, 3.5),
  north = c(0, 0, 1, 1, 2, 2, 0.2, 1.2),
  year = c(0, 0.5, 1.5, 2, 3, 4, 5, 6)
)
pars <- c(spatial_variance = 1.2, scale = 1.5, nugget_variance = 0.25)
temporal <- sp_temporal(~ year, scale = 2.5, time_units = "years")
control <- sp_control(
  spatial_variance = pars[["spatial_variance"]], scale = pars[["scale"]],
  nugget_variance = pars[["nugget_variance"]], smoothness = 0.5,
  estimate = FALSE, distance_units = "km"
)
fit <- sp_fit(y ~ x, dat, coords = ~ east + north, time = temporal,
              control = control)

# The fitted coefficients and likelihood reproduce direct space-time GLS.
x <- model.matrix(y ~ x, dat)
coordinates <- as.matrix(dat[c("east", "north")])
distance <- as.matrix(dist(coordinates))
time_lag <- abs(outer(dat$year, dat$year, "-"))
v <- pars[["spatial_variance"]] * exp(-distance / pars[["scale"]]) *
  exp(-time_lag / temporal$scale) +
  diag(pars[["nugget_variance"]], nrow(dat))
v_inverse <- solve(v)
beta <- solve(crossprod(x, v_inverse %*% x),
              crossprod(x, v_inverse %*% dat$y))
residual <- dat$y - as.vector(x %*% beta)
log_likelihood <- -0.5 * (nrow(dat) * log(2 * pi) +
  as.numeric(determinant(v, logarithm = TRUE)$modulus) +
  drop(crossprod(residual, v_inverse %*% residual)))
stopifnot(close(unname(coef(fit)), as.vector(beta)),
          close(as.numeric(logLik(fit)), log_likelihood),
          identical(fit$temporal_values, dat$year),
          identical(fit$temporal_scale, 2.5),
          identical(attr(logLik(fit), "df"), 2L))

# New-time point and joint predictions use both cross-time and new-time lags.
new <- data.frame(x = c(0.25, 1.75, 0.8), east = c(0.5, 2.5, 1.2),
                  north = c(0.25, 1.5, 1.8), year = c(0.25, 4.5, 8))
x_new <- model.matrix(~ x, new)
coordinates_new <- as.matrix(new[c("east", "north")])
cross_distance <- sqrt(outer(rowSums(coordinates^2),
  rowSums(coordinates_new^2), "+") - 2 * tcrossprod(coordinates,
                                                     coordinates_new))
cross_lag <- abs(outer(dat$year, new$year, "-"))
c_matrix <- pars[["spatial_variance"]] *
  exp(-cross_distance / pars[["scale"]]) * exp(-cross_lag / temporal$scale)
prediction_expected <- as.vector(x_new %*% beta +
  crossprod(c_matrix, v_inverse %*% residual))
q <- x_new - crossprod(c_matrix, v_inverse %*% x)
beta_vcov <- solve(crossprod(x, v_inverse %*% x))
new_distance <- as.matrix(dist(coordinates_new))
new_lag <- abs(outer(new$year, new$year, "-"))
k_new <- pars[["spatial_variance"]] *
  exp(-new_distance / pars[["scale"]]) * exp(-new_lag / temporal$scale)
joint_covariance <- k_new - crossprod(c_matrix, v_inverse %*% c_matrix) +
  q %*% beta_vcov %*% t(q)
prediction <- predict(fit, new, se.fit = TRUE)
draws <- sp_predict_draws(fit, new, nsim = 4, type = "link", seed = 17)
stopifnot(close(prediction$fit, prediction_expected),
          close(prediction$se.fit^2, diag(joint_covariance)),
          close(draws$predictive_covariance, joint_covariance))

# Simulation and the fitted variogram are exact consequences of the same V.
set.seed(19)
z <- matrix(rnorm(nrow(dat) * 2), nrow(dat))
simulation_expected <- sweep(crossprod(chol(v), z), 1L,
                             as.vector(x %*% beta), "+")
stopifnot(close(as.matrix(simulate(fit, nsim = 2, seed = 19)),
                simulation_expected))
variogram <- sp_variogram(fit, bins = 3)
triangle <- lower.tri(distance)
breaks <- seq(0, max(distance[triangle]), length.out = 4)
group <- cut(distance[triangle], breaks, include.lowest = TRUE, right = TRUE)
pair_fitted <- pars[["nugget_variance"]] + pars[["spatial_variance"]] *
  (1 - exp(-distance[triangle] / pars[["scale"]]) *
     exp(-time_lag[triangle] / temporal$scale))
stopifnot(close(variogram$values$fitted,
                as.numeric(tapply(pair_fitted, group, mean))))

# A fixed temporal scale is also carried through non-Gaussian Laplace fitting.
binary <- transform(dat, y = c(0, 1, 0, 1, 1, 0, 1, 1))
binary_fit <- sp_fit(y ~ x, binary, family = binomial(),
  coords = ~ east + north, time = temporal,
  control = sp_control(spatial_variance = 0.7, scale = 1.4,
    nugget_variance = 0.2, estimate = FALSE, distance_units = "km"))
binary_prediction <- predict(binary_fit, new, type = "response",
                             se.fit = TRUE)
binary_draws <- sp_predict_draws(binary_fit, new, nsim = 4, seed = 4)
stopifnot(all(is.finite(binary_prediction$fit)),
          all(binary_prediction$fit > 0 & binary_prediction$fit < 1),
          all(is.finite(binary_draws$predictive_covariance)))
binary_mcml <- sp_mcml(binary_fit,
  control = sp_mcml_control(nsim = 100, maxit = 20), seed = 23)
stopifnot(identical(binary_mcml$temporal_values, binary$year),
          identical(binary_mcml$temporal_scale, temporal$scale),
          all(is.finite(binary_mcml$mcml$log_weights)))

# Reporting and refitting retain the temporal contract.
temporal_row <- sp_model_table(fit)
temporal_row <- temporal_row[temporal_row$component == "temporal", ]
stopifnot(nrow(temporal_row) == 1L, temporal_row$term == "scale",
          temporal_row$status == "fixed", temporal_row$units == "years")
cross_validation <- sp_cross_validate(fit, folds = rep(1:2, 4))
stopifnot(all(is.finite(cross_validation$predictions$predicted)))

# Gaussian temporal-scale estimation is positive and counted in model df.
estimated <- sp_fit(y ~ x, dat, coords = ~ east + north,
  time = sp_temporal(~ year, scale = 2, time_units = "years", estimate = TRUE),
  control = sp_control(distance_units = "km", maxit = 200L))
stopifnot(is.finite(estimated$temporal_scale), estimated$temporal_scale > 0,
          isTRUE(estimated$temporal_scale_estimated),
          identical(attr(logLik(estimated), "df"), 6L),
          "temporal_scale" %in% rownames(estimated$variance_parameter_vcov))
estimated_recovery <- suppressWarnings(sp_simulation_recovery(
  estimated, nsim = 1, seed = 24
))
stopifnot("temporal_scale" %in% names(estimated_recovery$truth),
          "temporal_scale" %in% colnames(estimated_recovery$estimates),
          "temporal_scale" %in% colnames(
            estimated_recovery$standard_errors
          ),
          identical(unname(estimated_recovery$truth[["temporal_scale"]]),
                    unname(estimated$temporal_scale)))

# Binomial and Poisson Laplace fits can estimate the same positive temporal
# scale; MCML then re-estimates it on the log scale with the exact same
# space-time prior construction.
binary_temporal <- suppressWarnings(sp_fit(
  y ~ x, binary, family = binomial(), coords = ~ east + north,
  time = sp_temporal(~ year, 2, "years", estimate = TRUE),
  control = sp_control(distance_units = "km", maxit = 80L)
))
binary_temporal_table <- sp_model_table(binary_temporal)
binary_temporal_row <- binary_temporal_table$component == "temporal"
stopifnot(isTRUE(binary_temporal$temporal_scale_estimated),
          is.numeric(binary_temporal$temporal_scale),
          is.null(names(binary_temporal$temporal_scale)),
          binary_temporal$temporal_scale > 0,
          attr(logLik(binary_temporal), "df") == 6L,
          binary_temporal_table$status[binary_temporal_row] == "estimated")
expected_binary_sigma <- binary_temporal$covariance[["spatial_variance"]] *
  exp(-distance / binary_temporal$covariance[["scale"]]) *
  exp(-time_lag / binary_temporal$temporal_scale) +
  diag(binary_temporal$covariance[["nugget_variance"]], nrow(binary))
stopifnot(close(binary_temporal$fitted_state$sigma,
                expected_binary_sigma, tolerance = 1e-8))
binary_temporal_mcml <- suppressWarnings(sp_mcml(
  binary_temporal,
  control = sp_mcml_control(nsim = 100L, maxit = 2L,
                            ess_warning = .001), seed = 25L
))
binary_temporal_repeats <- suppressWarnings(sp_mcml_replicates(
  binary_temporal, nrep = 2L,
  control = sp_mcml_control(nsim = 100L, maxit = 1L,
                            ess_warning = .001), seed = 26L
))
stopifnot(binary_temporal_mcml$temporal_scale > 0,
          isTRUE(binary_temporal_mcml$temporal_scale_estimated),
          nrow(binary_temporal_mcml$parameter_vcov) == 6L,
          identical(binary_temporal_mcml$temporal$scale,
                    binary_temporal_mcml$temporal_scale),
          "temporal::scale" %in% binary_temporal_repeats$natural_summary$term,
          identical(dim(binary_temporal_repeats$natural_mc_covariance),
                    c(6L, 6L)))

# Invalid trial states receive a finite objective sentinel, so a difficult
# bounded line search returns an auditable fitted state/warning rather than
# failing with a non-finite-value error from optim().
robust_data <- expand.grid(z = c(-1, 0, 1, 2), group = seq_len(8L))
robust_data$east <- seq_len(nrow(robust_data))
robust_data$north <- rep(0:3, 8L)
robust_data$year <- rep(0:3, 8L) + rep(0:7, each = 4L) * .2
set.seed(42)
robust_data$y <- rpois(nrow(robust_data), exp(.2 + .15 * robust_data$z))
robust_fit <- suppressWarnings(sp_fit(
  y ~ z, robust_data, family = poisson(), coords = ~ east + north,
  time = sp_temporal(~ year, 1, "years", estimate = TRUE),
  control = sp_control(spatial_variance = .2, scale = 3,
    nugget_variance = .1, distance_units = "km", maxit = 60L)
))
stopifnot(is.finite(robust_fit$log_likelihood),
          robust_fit$temporal_scale > 0,
          !is.null(robust_fit$fitted_state$chol_sigma))

expect_error(sp_fit(y ~ x, dat, spatial = FALSE, time = temporal),
             "requires `spatial = TRUE`")
expect_error(predict(fit, transform(new, year = NA_real_)),
             "finite numeric")
constant_time <- transform(dat, year = 1)
expect_error(sp_fit(y ~ x, constant_time, coords = ~ east + north,
  time = sp_temporal(~ year, 2, "years", estimate = TRUE),
  control = sp_control(distance_units = "km")), "two distinct observed times")
constant_binary <- transform(binary, year = 1)
expect_error(sp_fit(y ~ x, constant_binary, family = binomial(),
  coords = ~ east + north,
  time = sp_temporal(~ year, 2, "years", estimate = TRUE),
  control = sp_control(distance_units = "km")), "two distinct observed times")

cat("Fitted space-time model tests passed.\n")
