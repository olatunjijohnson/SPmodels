if (file.exists(file.path("R", "control.R"))) {
  for (file in list.files("R", full.names = TRUE)) source(file)
} else {
  library(SPmodels)
}

expect_error <- function(expression, pattern) {
  error <- tryCatch({ force(expression); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}

# The control contract distinguishes the overall covariance optimisation from
# the nugget choice. A supplied zero fixes the nugget automatically, while a
# positive value remains a starting value unless explicitly fixed.
zero_control <- sp_control(nugget_variance = 0, distance_units = "unit")
stopifnot(isTRUE(zero_control$estimate), !zero_control$estimate_nugget)
positive_start <- sp_control(nugget_variance = 0.2, distance_units = "unit")
stopifnot(isTRUE(positive_start$estimate_nugget))
positive_fixed <- sp_control(
  nugget_variance = 0.2, estimate_nugget = FALSE, distance_units = "unit"
)
stopifnot(!positive_fixed$estimate_nugget)
expect_error(sp_control(estimate_nugget = FALSE), "supplied")
expect_error(sp_control(nugget_variance = 0, estimate_nugget = TRUE),
             "positive starting value")
expect_error(sp_control(nugget_variance = -0.1), "non-negative")

set.seed(1201)
n <- 14L
dat <- data.frame(
  east = runif(n), north = runif(n), exposure = rnorm(n)
)
distance <- as.matrix(dist(as.matrix(dat[c("east", "north")])))
truth <- c(spatial_variance = 0.8, scale = 0.45, nugget_variance = 0)
sigma <- truth[["spatial_variance"]] * exp(-distance / truth[["scale"]])
dat$response <- as.vector(
  0.4 - 0.3 * dat$exposure + t(chol(sigma)) %*% rnorm(n)
)

# Dense Gaussian ML estimates the two spatial parameters and keeps the nugget
# exactly at the requested boundary value. Fixed parameters have no invented
# standard error and do not count as fitted likelihood parameters.
dense <- suppressWarnings(sp_fit(
  response ~ exposure, dat, coords = ~ east + north,
  control = sp_control(
    spatial_variance = 0.7, scale = 0.5, nugget_variance = 0,
    distance_units = "unit", maxit = 150L
  )
))
stopifnot(
  identical(unname(dense$covariance[["nugget_variance"]]), 0),
  identical(dense$covariance_estimated_names,
            c("spatial_variance", "scale")),
  !"log_nugget_variance" %in% names(dense$optimisation$par),
  is.na(dense$covariance_vcov["nugget_variance", "nugget_variance"]),
  identical(attr(logLik(dense), "df"), dense$rank + 2L),
  all(is.finite(predict(dense, dat[1:3, ], se.fit = TRUE)$fit)),
  nrow(simulate(dense, nsim = 2L, seed = 9L)) == n
)
dense_table <- sp_model_table(dense)
dense_nugget <- dense_table[
  dense_table$component == "covariance" &
    dense_table$term == "nugget_variance", , drop = FALSE
]
stopifnot(identical(dense_nugget$status, "fixed"),
          identical(dense_nugget$estimate, 0),
          is.na(dense_nugget$std_error))
uncertain_draws <- sp_predict_draws(
  dense, dat[1:3, ], nsim = 4L, seed = 10L,
  parameter_uncertainty = TRUE
)
stopifnot(identical(dim(uncertain_draws$draws), c(4L, 3L)),
          all(is.finite(uncertain_draws$draws)))

# A positive nugget can likewise be held fixed while the other covariance
# parameters are estimated.
fixed_positive <- suppressWarnings(sp_fit(
  response ~ exposure, dat, coords = ~ east + north,
  control = sp_control(
    spatial_variance = 0.7, scale = 0.5, nugget_variance = 0.15,
    estimate_nugget = FALSE, distance_units = "unit", maxit = 150L
  )
))
stopifnot(
  identical(unname(fixed_positive$covariance[["nugget_variance"]]), 0.15),
  identical(fixed_positive$covariance_estimated_names,
            c("spatial_variance", "scale")),
  !"log_nugget_variance" %in% names(fixed_positive$optimisation$par)
)

# The nearest-predecessor engine obeys the same control contract.
neighbor <- suppressWarnings(sp_fit_neighbor(
  response ~ exposure, dat, coords = ~ east + north, neighbors = 4L,
  control = sp_control(
    spatial_variance = 0.7, scale = 0.5, nugget_variance = 0,
    distance_units = "unit", maxit = 100L
  )
))
stopifnot(
  identical(unname(neighbor$covariance[["nugget_variance"]]), 0),
  identical(neighbor$covariance_estimated_names,
            c("spatial_variance", "scale")),
  !"log_nugget_variance" %in% names(neighbor$optimisation$par),
  identical(attr(logLik(neighbor), "df"), neighbor$rank + 2L)
)

# Laplace binomial fitting also keeps zero out of the log-parameter vector.
binary <- dat
probability <- plogis(-0.25 + 0.6 * binary$exposure)
binary$response <- rbinom(n, 1L, probability)
binomial_fit <- suppressWarnings(sp_fit(
  response ~ exposure, binary, family = binomial(),
  coords = ~ east + north,
  control = sp_control(
    spatial_variance = 0.5, scale = 0.5, nugget_variance = 0,
    distance_units = "unit", maxit = 80L
  )
))
stopifnot(
  identical(unname(binomial_fit$covariance[["nugget_variance"]]), 0),
  identical(binomial_fit$covariance_estimated_names,
            c("spatial_variance", "scale")),
  !"log_nugget_variance" %in% names(binomial_fit$optimisation$par),
  identical(attr(logLik(binomial_fit), "df"), binomial_fit$rank + 2L)
)

# The MCML refinement uses the fitted parameter map rather than reintroducing
# a log nugget internally.
mcml_fit <- suppressWarnings(sp_mcml(
  binomial_fit,
  control = sp_mcml_control(
    nsim = 120L, maxit = 25L, ess_warning = 0.001
  ),
  seed = 19L
))
stopifnot(
  identical(unname(mcml_fit$covariance[["nugget_variance"]]), 0),
  identical(mcml_fit$covariance_estimated_names,
            c("spatial_variance", "scale")),
  !"log_nugget_variance" %in% names(mcml_fit$optimisation$par),
  identical(nrow(mcml_fit$parameter_vcov), mcml_fit$rank + 2L)
)

# Fully fixed zero remains supported and is algebraically checked against the
# exact Gaussian likelihood, rather than merely checking that the call runs.
fixed_zero <- sp_fit(
  response ~ exposure, dat, coords = ~ east + north,
  control = sp_control(
    spatial_variance = truth[["spatial_variance"]],
    scale = truth[["scale"]], nugget_variance = 0,
    estimate = FALSE, distance_units = "unit"
  )
)
x <- model.matrix(response ~ exposure, dat)
inverse <- solve(sigma)
beta <- solve(crossprod(x, inverse %*% x),
              crossprod(x, inverse %*% dat$response))
residual <- dat$response - as.vector(x %*% beta)
expected <- -0.5 * (
  n * log(2 * pi) + as.numeric(determinant(sigma, logarithm = TRUE)$modulus) +
    drop(crossprod(residual, inverse %*% residual))
)
stopifnot(
  isTRUE(all.equal(as.numeric(logLik(fixed_zero)), expected,
                   tolerance = 1e-10)),
  identical(attr(logLik(fixed_zero), "df"), fixed_zero$rank)
)
