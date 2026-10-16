if (file.exists(file.path("R", "random-effects.R"))) {
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
  y = c(0, 1, 0, 1, 1, 0, 1, 1),
  x = seq(-1, 1, length.out = 8),
  east = c(0, 1, 0, 1, 2, 3, 2, 3),
  north = c(0, 0, 1, 1, 0, 0, 1, 1),
  clinic = factor(rep(c("A", "B"), each = 4))
)
control <- sp_control(
  spatial_variance = 0.7, scale = 1.3, nugget_variance = 0.2,
  random_variance = c(clinic = 0.45), estimate = FALSE,
  distance_units = "km"
)
fit <- sp_fit(y ~ x, dat, family = binomial(), coords = ~ east + north,
              random = ~ clinic, control = control)

coordinates <- as.matrix(dat[c("east", "north")])
distance <- as.matrix(dist(coordinates))
z <- model.matrix(~ clinic - 1, dat)
sigma <- 0.7 * exp(-distance / 1.3) + 0.45 * tcrossprod(z) +
  diag(0.2, nrow(dat))
stopifnot(close(fit$fitted_state$sigma, sigma),
          identical(fit$random_variance, c(clinic = 0.45)),
          attr(logLik(fit), "df") == length(coef(fit)))

new <- data.frame(x = c(-0.2, 0.7), east = c(0.5, 2.5),
                  north = c(0.4, 0.6),
                  clinic = factor(c("A", "B"), levels = levels(dat$clinic)))
x_new <- model.matrix(~ x, new)
coordinates_new <- as.matrix(new[c("east", "north")])
cross_distance <- sqrt(outer(rowSums(coordinates^2), rowSums(coordinates_new^2),
  "+") - 2 * tcrossprod(coordinates, coordinates_new))
z_new <- model.matrix(~ clinic - 1, new)
cross_covariance <- 0.7 * exp(-cross_distance / 1.3) +
  0.45 * tcrossprod(z, z_new)
precision_cross <- solve(sigma, cross_covariance)
mean_expected <- as.vector(x_new %*% coef(fit) +
  crossprod(cross_covariance, solve(sigma, fit$latent_mode$mode)))
variance_expected <- 0.7 + 0.45 -
  colSums(cross_covariance * precision_cross) +
  colSums(precision_cross * (fit$latent_mode$covariance %*% precision_cross))
prediction <- predict(fit, new, type = "link", se.fit = TRUE)
stopifnot(close(prediction$fit, mean_expected),
          close(prediction$se.fit^2, variance_expected))

new_distance <- as.matrix(dist(coordinates_new))
k_new <- 0.7 * exp(-new_distance / 1.3) + 0.45 * tcrossprod(z_new)
joint_expected <- k_new - crossprod(cross_covariance, precision_cross) +
  crossprod(precision_cross,
            fit$latent_mode$covariance %*% precision_cross)
joint <- sp_predict_draws(fit, new, nsim = 4, type = "link", seed = 51)
stopifnot(close(joint$predictive_covariance, joint_expected))

# A shared unseen level is integrated out: it has no covariance with training
# random effects, but contributes its variance to both new rows and their
# joint covariance.
new_unseen <- new
new_unseen$clinic <- c("C", "C")
cross_unseen <- 0.7 * exp(-cross_distance / 1.3)
precision_unseen <- solve(sigma, cross_unseen)
mean_unseen <- as.vector(x_new %*% coef(fit) +
  crossprod(cross_unseen, solve(sigma, fit$latent_mode$mode)))
variance_unseen <- 0.7 + 0.45 -
  colSums(cross_unseen * precision_unseen) +
  colSums(precision_unseen *
    (fit$latent_mode$covariance %*% precision_unseen))
prediction_unseen <- predict(
  fit, new_unseen, type = "link", se.fit = TRUE, unseen = "marginal"
)
joint_unseen_expected <- 0.7 * exp(-new_distance / 1.3) +
  0.45 * matrix(1, 2, 2) -
  crossprod(cross_unseen, precision_unseen) +
  crossprod(precision_unseen,
            fit$latent_mode$covariance %*% precision_unseen)
joint_unseen <- sp_predict_draws(
  fit, new_unseen, nsim = 4, type = "link", seed = 56,
  unseen = "marginal"
)
stopifnot(close(prediction_unseen$fit, mean_unseen),
          close(prediction_unseen$se.fit^2, variance_unseen),
          close(joint_unseen$predictive_covariance,
                joint_unseen_expected))

# The fitted latent covariance automatically drives seeded response simulation.
set.seed(52)
simulation_z <- matrix(rnorm(nrow(dat) * 2), nrow(dat))
simulation_eta <- sweep(crossprod(chol(sigma), simulation_z), 1L,
                        as.vector(fit$x %*% coef(fit)), "+")
simulation_expected <- matrix(rbinom(length(simulation_eta), 1,
                                     plogis(simulation_eta)), nrow(dat))
stopifnot(close(as.matrix(simulate(fit, nsim = 2, seed = 52)),
                simulation_expected, 0))

random_row <- sp_model_table(fit)
random_row <- random_row[random_row$component == "random", ]
stopifnot(nrow(random_row) == 1L, random_row$term == "clinic",
          random_row$status == "fixed")
cv <- sp_cross_validate(fit, folds = c(1, 2, 2, 1, 1, 2, 2, 1))
stopifnot(all(is.finite(cv$predictions$predicted)))

# Joint estimation adds the random variance on the log scale.
estimated <- suppressWarnings(sp_fit(
  y ~ x, dat, family = binomial(), coords = ~ east + north,
  random = ~ clinic,
  control = sp_control(random_variance = c(clinic = 0.3), estimate = TRUE,
                       distance_units = "km", maxit = 100)
))
stopifnot(estimated$random_variance[["clinic"]] > 0,
          identical(dim(estimated$random_variance_vcov), c(1L, 1L)),
          attr(logLik(estimated), "df") == length(coef(estimated)) + 4L)
estimated_mcml <- suppressWarnings(sp_mcml(
  estimated, sp_mcml_control(nsim = 100, maxit = 25, ess_warning = .001),
  seed = 54
))
stopifnot(estimated_mcml$random_variance[["clinic"]] > 0,
          "log_random::clinic" %in% names(estimated_mcml$optimisation$par),
          identical(dim(estimated_mcml$random_variance_vcov), c(1L, 1L)),
          all(is.finite(estimated_mcml$mcml$log_weights)))

expect_error(predict(fit, transform(new, clinic = factor(c("C", "B")))),
             "unseen")
mcml_fit <- suppressWarnings(sp_mcml(
  fit, sp_mcml_control(nsim = 100, maxit = 25, ess_warning = .001), seed = 53
))
stopifnot(all(is.finite(mcml_fit$mcml$log_weights)),
          identical(mcml_fit$random_variance, fit$random_variance),
          close(mcml_fit$fitted_state$sigma, sigma))
parameter_draws <- sp_predict_draws(
  fit, new, nsim = 3, type = "link", parameter_uncertainty = TRUE, seed = 55
)
stopifnot(identical(dim(parameter_draws$parameter_draws), c(3L, 2L)),
          all(is.finite(parameter_draws$draws)))
expect_error(sp_fit(y ~ x, dat, family = binomial(),
  coords = ~ east + north, random = ~ clinic,
  control = sp_control(spatial_variance = .7, scale = 1.3,
    nugget_variance = .2, estimate = FALSE, distance_units = "km")),
  "random_variance")

cat("Non-Gaussian grouped random-intercept tests passed.\n")
