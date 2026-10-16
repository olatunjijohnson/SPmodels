if (file.exists(file.path("R", "predict-draws.R"))) {
  for (file in list.files("R", full.names = TRUE)) source(file)
} else {
  library(SPmodels)
}

close <- function(a, b, tolerance = 1e-10) {
  isTRUE(all.equal(a, b, tolerance = tolerance, check.attributes = FALSE))
}
expect_error <- function(expression, pattern) {
  error <- tryCatch({ force(expression); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}

dat <- data.frame(
  y = c(1.2, 0.7, 2.1, 1.8, 3.0, 2.6),
  x = c(-1, 0, 0.5, 1, 1.5, 2),
  east = c(0, 1, 0, 2, 1, 3), north = c(0, 0, 1, 1, 2, 2)
)
control <- sp_control(spatial_variance = 1.4, scale = 1.7,
                      nugget_variance = 0.3, smoothness = 0.5,
                      estimate = FALSE, distance_units = "km")
fit <- sp_fit(y ~ x, dat, coords = ~ east + north, control = control)
new <- data.frame(x = c(0.25, 1.75, 0.8), east = c(0.5, 2.5, 1.2),
                  north = c(0.25, 1.5, 1.8))

# Reconstruct the complete universal-kriging covariance and seeded transform.
x_new <- model.matrix(~ x, new)
coordinates <- as.matrix(dat[c("east", "north")])
coordinates_new <- as.matrix(new[c("east", "north")])
distance <- as.matrix(dist(coordinates))
cross_distance <- sqrt(outer(rowSums(coordinates^2),
  rowSums(coordinates_new^2), "+") - 2 * tcrossprod(coordinates,
                                                     coordinates_new))
new_distance <- as.matrix(dist(coordinates_new))
v <- 1.4 * exp(-distance / 1.7) + diag(0.3, nrow(dat))
c_matrix <- 1.4 * exp(-cross_distance / 1.7)
k_new <- 1.4 * exp(-new_distance / 1.7)
v_inverse <- solve(v)
x <- model.matrix(~ x, dat)
beta_vcov <- solve(crossprod(x, v_inverse %*% x))
q <- x_new - crossprod(c_matrix, v_inverse %*% x)
covariance_expected <- k_new - crossprod(c_matrix, v_inverse %*% c_matrix) +
  q %*% beta_vcov %*% t(q)
mean_expected <- predict(fit, new)

set.seed(707)
seed_before <- .Random.seed
draws <- sp_predict_draws(fit, new, nsim = 12, type = "link",
                          target = "latent", seed = 22)
stopifnot(identical(.Random.seed, seed_before))
set.seed(22)
z <- matrix(rnorm(12 * nrow(new)), nrow = 12)
draws_expected <- sweep(z %*% chol(covariance_expected), 2L,
                        mean_expected, "+")
stopifnot(close(draws$predictive_mean, mean_expected),
          close(draws$predictive_covariance, covariance_expected),
          close(draws$draws, draws_expected))

# Optional parameter draws are local, reproducible, recorded and alter the
# joint predictive draws. Fixed covariance leaves only beta uncertainty.
set.seed(909)
parameter_seed_before <- .Random.seed
parameter_draws_one <- sp_predict_draws(
  fit, new, nsim = 15, type = "link", target = "latent", seed = 31,
  parameter_uncertainty = TRUE
)
stopifnot(identical(.Random.seed, parameter_seed_before))
parameter_draws_two <- sp_predict_draws(
  fit, new, nsim = 15, type = "link", target = "latent", seed = 31,
  parameter_uncertainty = TRUE
)
stopifnot(identical(parameter_draws_one$draws, parameter_draws_two$draws),
          identical(parameter_draws_one$parameter_draws,
                    parameter_draws_two$parameter_draws),
          identical(dim(parameter_draws_one$parameter_draws), c(15L, 2L)),
          isTRUE(parameter_draws_one$parameter_uncertainty))

# Non-Gaussian joint link covariance has the same marginal variances as
# pointwise prediction, and response draws respect their support.
binomial_data <- transform(dat, y = c(0, 1, 0, 1, 1, 0))
binomial_fit <- sp_fit(y ~ x, binomial_data, family = binomial(),
  coords = ~ east + north,
  control = sp_control(spatial_variance = 0.7, scale = 1.4,
    nugget_variance = 0.2, estimate = FALSE, distance_units = "km"))
link_prediction <- predict(binomial_fit, new, type = "link", se.fit = TRUE,
                           target = "observation")
joint_link <- sp_predict_draws(binomial_fit, new, nsim = 20, type = "link",
                               target = "observation", seed = 8)
stopifnot(close(joint_link$predictive_mean, link_prediction$fit),
          close(diag(joint_link$predictive_covariance),
                link_prediction$se.fit^2))
probability_draws <- sp_predict_draws(
  binomial_fit, new, nsim = 20, type = "response", target = "latent", seed = 9
)$draws
outcome_draws <- sp_predict_draws(
  binomial_fit, new, nsim = 20, type = "response",
  target = "observation", seed = 9
)$draws
stopifnot(all(probability_draws > 0 & probability_draws < 1),
          all(outcome_draws %in% c(0, 1)))
binomial_parameter_draws <- sp_predict_draws(
  binomial_fit, new, nsim = 10, type = "response", target = "latent",
  seed = 12, parameter_uncertainty = TRUE
)
stopifnot(all(binomial_parameter_draws$draws > 0 &
              binomial_parameter_draws$draws < 1),
          ncol(binomial_parameter_draws$parameter_draws) == 2L)

poisson_data <- transform(dat, y = c(0, 1, 1, 2, 3, 4))
poisson_fit <- sp_fit(y ~ x, poisson_data, family = poisson(),
  coords = ~ east + north,
  control = sp_control(spatial_variance = 0.7, scale = 1.4,
    nugget_variance = 0.2, estimate = FALSE, distance_units = "km"))
count_draws <- sp_predict_draws(
  poisson_fit, new, nsim = 20, type = "response",
  target = "observation", seed = 10
)$draws
stopifnot(all(count_draws >= 0), all(count_draws == floor(count_draws)))

expect_error(sp_predict_draws(fit, new, nsim = 0), "positive integer")
cat("Joint predictive-draw tests passed.\n")
