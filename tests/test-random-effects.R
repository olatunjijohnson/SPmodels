if (file.exists(file.path("R", "random-effects.R"))) {
  for (file in list.files("R", full.names = TRUE)) source(file)
} else {
  library(SPmodels)
}

expect_error <- function(expr) {
  value <- try(force(expr), silent = TRUE)
  stopifnot(inherits(value, "try-error"))
}

dat <- data.frame(
  y = c(0.8, 1.4, 1.7, 2.5, 2.2, 3.1, 3.4, 4.0),
  x = c(-1, -.5, 0, .5, 1, 1.5, 2, 2.5),
  clinic = factor(rep(c("A", "B", "C", "D"), each = 2)),
  survey = factor(rep(c("early", "late"), 4)),
  east = c(0, .5, 1, 1.5, 2, 2.5, 3, 3.5),
  north = c(0, .2, .8, 1, 1.7, 1.9, 2.7, 2.9)
)
base_covariance <- c(spatial_variance = 0.7, scale = 1.3,
                     nugget_variance = 0.15)
random_variance <- c(clinic = 0.4, survey = 0.25)
control <- sp_control(
  smoothness = 0.5, spatial_variance = base_covariance[[1]],
  scale = base_covariance[[2]], nugget_variance = base_covariance[[3]],
  random_variance = random_variance, estimate = FALSE,
  distance_units = "km"
)
fit <- sp_fit(y ~ x, dat, coords = ~ east + north,
              random = ~ clinic + survey, control = control)

# Reconstruct the complete covariance, GLS fit and likelihood independently.
coordinates <- as.matrix(dat[c("east", "north")])
distance <- as.matrix(dist(coordinates))
z_clinic <- model.matrix(~ 0 + clinic, dat)
z_survey <- model.matrix(~ 0 + survey, dat)
v <- base_covariance[["spatial_variance"]] *
  exp(-distance / base_covariance[["scale"]]) +
  random_variance[["clinic"]] * tcrossprod(z_clinic) +
  random_variance[["survey"]] * tcrossprod(z_survey) +
  diag(base_covariance[["nugget_variance"]], nrow(dat))
x <- model.matrix(y ~ x, dat)
v_inverse <- solve(v)
beta <- solve(crossprod(x, v_inverse %*% x),
              crossprod(x, v_inverse %*% dat$y))
residual <- dat$y - as.vector(x %*% beta)
log_likelihood <- -0.5 * (nrow(dat) * log(2 * pi) +
  as.numeric(determinant(v, logarithm = TRUE)$modulus) +
  drop(crossprod(residual, v_inverse %*% residual)))
stopifnot(all.equal(unname(coef(fit)), as.vector(beta), tolerance = 1e-11))
stopifnot(all.equal(as.numeric(logLik(fit)), log_likelihood,
                    tolerance = 1e-11))
stopifnot(all.equal(unname(fit$fitted_state$v), unname(v),
                    tolerance = 1e-12))
stopifnot(identical(fit$random_variance, random_variance))
stopifnot(attr(logLik(fit), "df") == length(coef(fit)))

# Known-level prediction includes both random-intercept cross-covariances.
new <- data.frame(
  x = c(.25, 1.75), clinic = factor(c("A", "C"), levels = levels(dat$clinic)),
  survey = factor(c("late", "early"), levels = levels(dat$survey)),
  east = c(.75, 2.75), north = c(.4, 2.3),
  row.names = c("one", "two")
)
new_coordinates <- as.matrix(new[c("east", "north")])
cross_distance <- sqrt(outer(rowSums(coordinates^2),
                             rowSums(new_coordinates^2), "+") -
                         2 * tcrossprod(coordinates, new_coordinates))
z_clinic_new <- model.matrix(~ 0 + clinic, new)
z_survey_new <- model.matrix(~ 0 + survey, new)
cross <- base_covariance[["spatial_variance"]] *
  exp(-cross_distance / base_covariance[["scale"]]) +
  random_variance[["clinic"]] * z_clinic %*% t(z_clinic_new) +
  random_variance[["survey"]] * z_survey %*% t(z_survey_new)
x_new <- model.matrix(~ x, new)
expected_mean <- as.vector(x_new %*% beta +
  crossprod(cross, v_inverse %*% residual))
vinv_x <- v_inverse %*% x
mean_difference <- x_new - crossprod(cross, vinv_x)
expected_variance <- base_covariance[["spatial_variance"]] +
  sum(random_variance) - colSums(cross * (v_inverse %*% cross)) +
  rowSums((mean_difference %*% solve(crossprod(x, vinv_x))) *
            mean_difference)
prediction <- predict(fit, new, se.fit = TRUE)
stopifnot(all.equal(prediction$fit, expected_mean, tolerance = 1e-11))
stopifnot(all.equal(unname(prediction$se.fit^2), unname(expected_variance),
                    tolerance = 1e-11))

# Joint draws use the corresponding full new-location covariance.
joint <- sp_predict_draws(fit, new, nsim = 3, seed = 440,
                          target = "latent", type = "link")
new_distance <- as.matrix(dist(new_coordinates))
new_covariance <- base_covariance[["spatial_variance"]] *
  exp(-new_distance / base_covariance[["scale"]]) +
  random_variance[["clinic"]] * tcrossprod(z_clinic_new) +
  random_variance[["survey"]] * tcrossprod(z_survey_new) -
  crossprod(cross, v_inverse %*% cross) +
  mean_difference %*% solve(crossprod(x, vinv_x)) %*% t(mean_difference)
stopifnot(all.equal(joint$predictive_mean, expected_mean, tolerance = 1e-11))
stopifnot(all.equal(joint$predictive_covariance, new_covariance,
                    tolerance = 1e-11))

# Marginal prediction for a new level has zero training cross-covariance but
# retains its variance and sharing between prediction rows with that level.
new_unseen <- new
new_unseen$clinic <- c("E", "E")
cross_unseen <- base_covariance[["spatial_variance"]] *
  exp(-cross_distance / base_covariance[["scale"]]) +
  random_variance[["survey"]] * z_survey %*% t(z_survey_new)
unseen_mean <- as.vector(x_new %*% beta +
  crossprod(cross_unseen, v_inverse %*% residual))
unseen_difference <- x_new - crossprod(cross_unseen, vinv_x)
unseen_variance <- base_covariance[["spatial_variance"]] +
  sum(random_variance) -
  colSums(cross_unseen * (v_inverse %*% cross_unseen)) +
  rowSums((unseen_difference %*% solve(crossprod(x, vinv_x))) *
            unseen_difference)
unseen_prediction <- predict(
  fit, new_unseen, se.fit = TRUE, unseen = "marginal"
)
unseen_new_covariance <- base_covariance[["spatial_variance"]] *
  exp(-new_distance / base_covariance[["scale"]]) +
  random_variance[["clinic"]] * matrix(1, 2, 2) +
  random_variance[["survey"]] * tcrossprod(z_survey_new) -
  crossprod(cross_unseen, v_inverse %*% cross_unseen) +
  unseen_difference %*% solve(crossprod(x, vinv_x)) %*%
    t(unseen_difference)
unseen_joint <- sp_predict_draws(
  fit, new_unseen, nsim = 3, seed = 444, type = "link",
  unseen = "marginal"
)
stopifnot(all.equal(unseen_prediction$fit, unseen_mean, tolerance = 1e-11),
          all.equal(unname(unseen_prediction$se.fit^2),
                    unname(unseen_variance), tolerance = 1e-11),
          all.equal(unseen_joint$predictive_covariance,
                    unseen_new_covariance, tolerance = 1e-11))

# Simulation is exactly the fitted full-covariance Cholesky transform.
set.seed(441)
z <- matrix(rnorm(nrow(dat) * 2L), nrow(dat))
expected_simulation <- sweep(crossprod(chol(v), z), 1L,
                             as.vector(x %*% beta), "+")
simulation <- simulate(fit, nsim = 2, seed = 441)
stopifnot(all.equal(unname(as.matrix(simulation)), unname(expected_simulation),
                    tolerance = 1e-12))

model_table <- sp_model_table(fit)
stopifnot(identical(model_table$term[model_table$component == "random"],
                    c("clinic", "survey")))
stopifnot(all(model_table$status[model_table$component == "random"] == "fixed"))

# Fitted variogram adds the exact average group-difference contribution by bin.
variogram <- sp_variogram(fit, bins = 3)
triangle <- lower.tri(distance)
pair_distance <- distance[triangle]
breaks <- seq(0, max(pair_distance), length.out = 4L)
bin <- cut(pair_distance, breaks = breaks, include.lowest = TRUE, right = TRUE)
random_covariance <- random_variance[["clinic"]] * tcrossprod(z_clinic) +
  random_variance[["survey"]] * tcrossprod(z_survey)
random_gamma <- 0.5 * outer(diag(random_covariance),
                            diag(random_covariance), "+") - random_covariance
random_by_bin <- as.numeric(tapply(random_gamma[triangle], bin, mean))
base_fitted <- base_covariance[["nugget_variance"]] +
  base_covariance[["spatial_variance"]] *
  (1 - exp(-variogram$values$distance / base_covariance[["scale"]]))
stopifnot(all.equal(variogram$values$fitted,
                    base_fitted + random_by_bin, tolerance = 1e-12))

# Refit workflows retain the grouping specification when all held-out levels
# remain represented in their corresponding training folds.
cross_validation <- sp_cross_validate(fit, folds = c(1, 2, 2, 1, 1, 2, 2, 1))
stopifnot(all(is.finite(cross_validation$predictions$predicted)))
recovery <- sp_simulation_recovery(fit, nsim = 1, seed = 442)
stopifnot(is.na(recovery$failures[1L]))
expect_error(predict(fit, transform(new, clinic = factor(c("new", "C")))))
group_parameter_draws <- sp_predict_draws(
  fit, new, nsim = 5, parameter_uncertainty = TRUE, seed = 443,
  type = "link"
)
stopifnot(identical(dim(group_parameter_draws$parameter_draws), c(5L, 2L)),
          all(is.finite(group_parameter_draws$draws)))
expect_error(sp_fit(y ~ x, dat, coords = ~ east + north,
                    random = ~ clinic:survey, control = control))
expect_error(sp_fit(y ~ x, dat, coords = ~ east + north,
                    random = ~ clinic,
                    control = sp_control(
                      spatial_variance = .7, scale = 1.3,
                      nugget_variance = .15, estimate = FALSE,
                      distance_units = "km")))

# Estimated random variance is positive and represented in its Hessian block.
estimated_control <- sp_control(
  smoothness = .5, spatial_variance = .7, scale = 1.3,
  nugget_variance = .15, random_variance = c(clinic = .3),
  estimate = TRUE, distance_units = "km", maxit = 100
)
estimated <- suppressWarnings(sp_fit(
  y ~ x, dat, coords = ~ east + north, random = ~ clinic,
  control = estimated_control
))
stopifnot(length(estimated$random_variance) == 1L,
          estimated$random_variance > 0,
          identical(dim(estimated$random_variance_vcov), c(1L, 1L)),
          attr(logLik(estimated), "df") == length(coef(estimated)) + 4L)
estimated_recovery <- suppressWarnings(sp_simulation_recovery(
  estimated, nsim = 1, seed = 445
))
stopifnot("random::clinic" %in% names(estimated_recovery$truth),
          "random::clinic" %in% colnames(estimated_recovery$estimates),
          "random::clinic" %in% colnames(estimated_recovery$standard_errors),
          identical(unname(estimated_recovery$truth[["random::clinic"]]),
                    unname(estimated$random_variance[["clinic"]])))

cat("Gaussian grouped random-intercept tests passed.\n")
