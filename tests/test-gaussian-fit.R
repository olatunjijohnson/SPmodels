if (file.exists(file.path("R", "fit.R"))) {
  for (file in list.files("R", full.names = TRUE)) source(file)
} else {
  library(SPmodels)
}

dat <- data.frame(
  response = c(1.2, 0.7, 2.1, 1.8, 3.0, 2.6),
  exposure = c(-1, 0, 0.5, 1, 1.5, 2),
  east = c(0, 1, 0, 2, 1, 3),
  north = c(0, 0, 1, 1, 2, 2)
)
pars <- c(spatial_variance = 1.4, scale = 1.7, nugget_variance = 0.3)
ctl <- sp_control(
  smoothness = 0.5,
  spatial_variance = pars[["spatial_variance"]],
  scale = pars[["scale"]],
  nugget_variance = pars[["nugget_variance"]],
  estimate = FALSE,
  distance_units = "km"
)
fit <- sp_fit(response ~ exposure, data = dat, coords = ~ east + north,
              control = ctl)

# Direct generalized least-squares and likelihood calculation.
x <- model.matrix(response ~ exposure, dat)
coordinates <- as.matrix(dat[c("east", "north")])
distance <- as.matrix(dist(coordinates))
v <- pars[["spatial_variance"]] *
  exp(-distance / pars[["scale"]]) +
  diag(pars[["nugget_variance"]], nrow(dat))
v_inverse <- solve(v)
beta_expected <- solve(crossprod(x, v_inverse %*% x),
                       crossprod(x, v_inverse %*% dat$response))
stopifnot(all.equal(unname(coef(fit)), as.vector(beta_expected),
                    tolerance = 1e-11))
residual_expected <- dat$response - as.vector(x %*% beta_expected)
log_likelihood_expected <- -0.5 * (
  nrow(dat) * log(2 * pi) + as.numeric(determinant(v, logarithm = TRUE)$modulus) +
    drop(crossprod(residual_expected, v_inverse %*% residual_expected))
)
stopifnot(all.equal(as.numeric(logLik(fit)), log_likelihood_expected,
                    tolerance = 1e-11))
stopifnot(identical(attr(logLik(fit), "df"), 2L))

# Direct universal-kriging mean and variance calculation.
new <- data.frame(exposure = c(0.25, 1.75), east = c(0.5, 2.5),
                  north = c(0.25, 1.5))
x_new <- model.matrix(~ exposure, new)
cross_distance <- sqrt(outer(rowSums(coordinates^2),
                             rowSums(as.matrix(new[c("east", "north")])^2), "+") -
                         2 * tcrossprod(coordinates,
                                       as.matrix(new[c("east", "north")])))
cross_covariance <- pars[["spatial_variance"]] *
  exp(-cross_distance / pars[["scale"]])
prediction_expected <- as.vector(
  x_new %*% beta_expected +
    crossprod(cross_covariance, v_inverse %*% residual_expected)
)
mean_difference <- x_new - crossprod(cross_covariance, v_inverse %*% x)
beta_vcov <- solve(crossprod(x, v_inverse %*% x))
latent_variance_expected <- pars[["spatial_variance"]] -
  colSums(cross_covariance * (v_inverse %*% cross_covariance)) +
  rowSums((mean_difference %*% beta_vcov) * mean_difference)
prediction <- predict(fit, newdata = new, se.fit = TRUE)
stopifnot(all.equal(prediction$fit, prediction_expected, tolerance = 1e-11))
stopifnot(all.equal(prediction$se.fit, unname(sqrt(latent_variance_expected)),
                    tolerance = 1e-11))
observation <- predict(fit, newdata = new, se.fit = TRUE,
                       target = "observation")
stopifnot(all.equal(observation$se.fit,
                    unname(sqrt(latent_variance_expected + pars[["nugget_variance"]])),
                    tolerance = 1e-11))

# A row permutation changes no fitted quantity or prediction.
permutation <- c(6, 2, 4, 1, 5, 3)
fit_permuted <- sp_fit(response ~ exposure, data = dat[permutation, ],
                       coords = ~ east + north, control = ctl)
stopifnot(all.equal(unname(coef(fit_permuted)), unname(coef(fit)),
                    tolerance = 1e-11))
stopifnot(all.equal(predict(fit_permuted, newdata = new),
                    predict(fit, newdata = new), tolerance = 1e-11))

# Estimated covariance parameters remain positive and optimisation completes.
estimated <- sp_fit(response ~ exposure, data = dat, coords = ~ east + north,
                    control = sp_control(distance_units = "km", maxit = 200L))
stopifnot(estimated$optimisation$convergence == 0L)
stopifnot(all(is.finite(estimated$covariance)), all(estimated$covariance > 0))
stopifnot(identical(attr(logLik(estimated), "df"), 5L))
stopifnot(identical(dim(estimated$covariance_vcov), c(3L, 3L)))
stopifnot(isSymmetric(estimated$covariance_vcov))
stopifnot(identical(colnames(summary(estimated)$covariance_table),
                    c("Estimate", "Std. Error")))

# Contract failures are explicit.
expect_error <- function(expression, pattern) {
  error <- tryCatch({ force(expression); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}
expect_error(sp_fit(response ~ exposure, dat, coords = ~ east + north),
             "explicit coordinate units")
expect_error(sp_fit(response ~ exposure, dat, coords = ~ east,
                    control = ctl), "exactly two")
expect_error(sp_fit(response ~ exposure, dat, coords = ~ east + north,
                    family = binomial(), control = ctl), "binomial response")
expect_error(sp_fit(response ~ exposure, transform(dat, east = 1, north = 1),
                    coords = ~ east + north, control = ctl), "distinct")

# Simulation is the direct Gaussian transform and a supplied seed is local.
set.seed(101)
seed_before <- .Random.seed
simulated <- simulate(fit, nsim = 3, seed = 44)
stopifnot(identical(.Random.seed, seed_before))
set.seed(44)
z <- matrix(rnorm(nrow(dat) * 3), nrow = nrow(dat))
simulation_expected <- sweep(
  crossprod(fit$fitted_state$chol_v, z), 1L,
  as.vector(fit$x %*% coef(fit) + fit$offset), "+"
)
stopifnot(all.equal(as.matrix(simulated), simulation_expected,
                    check.attributes = FALSE, tolerance = 1e-12))

# Variogram bins and diagnostics reproduce their defining calculations.
variogram <- sp_variogram(fit, bins = 3)
stopifnot(inherits(variogram, "spmodels_variogram"))
stopifnot(sum(variogram$values$pairs) == choose(nrow(dat), 2))
stopifnot(all(variogram$values$fitted >= 0))
diagnostics <- sp_diagnostics(fit)
decorrelated_expected <- forwardsolve(t(fit$fitted_state$chol_v),
                                      fit$residual_mean)
stopifnot(all.equal(diagnostics$decorrelated_residual,
                    as.vector(decorrelated_expected), tolerance = 1e-12))

# Projected point sf data supply coordinates, CRS and units without duplication.
if (requireNamespace("sf", quietly = TRUE)) {
  dat_sf <- sf::st_as_sf(dat, coords = c("east", "north"), crs = 3857)
  fit_sf <- sp_fit(response ~ exposure, data = dat_sf,
                   control = sp_control(
                     smoothness = 0.5,
                     spatial_variance = pars[["spatial_variance"]],
                     scale = pars[["scale"]],
                     nugget_variance = pars[["nugget_variance"]],
                     estimate = FALSE
                   ))
  stopifnot(identical(fit_sf$coordinate_source, "sf"))
  stopifnot(identical(fit_sf$distance_units, "metre"))
  stopifnot(length(predict(fit_sf, dat_sf[1:2, ])) == 2L)

  # British National Grid coordinates are metres; requesting kilometres
  # converts both training and prediction coordinates before covariance use.
  projected <- sf::st_as_sf(
    data.frame(response = dat$response, exposure = dat$exposure,
      east = 371000 + dat$east * 1000,
      north = 408000 + dat$north * 1000),
    coords = c("east", "north"), crs = 27700
  )
  projected_km <- sp_fit(
    response ~ exposure, projected,
    control = sp_control(smoothness = 0.5, spatial_variance = 1.4,
      scale = 1.7, nugget_variance = 0.3, estimate = FALSE,
      distance_units = "km")
  )
  projected_m <- sp_fit(
    response ~ exposure, projected,
    control = sp_control(smoothness = 0.5, spatial_variance = 1.4,
      scale = 1700, nugget_variance = 0.3, estimate = FALSE,
      distance_units = "m")
  )
  stopifnot(identical(projected_km$distance_units, "km"),
    isTRUE(all.equal(projected_km$coordinates,
      sf::st_coordinates(projected) / 1000)),
    all.equal(as.numeric(logLik(projected_km)),
      as.numeric(logLik(projected_m)), tolerance = 1e-10),
    all.equal(predict(projected_km, projected[1:2, ]),
      predict(projected_m, projected[1:2, ]), tolerance = 1e-10))
  expect_error(sp_fit(response ~ exposure, dat_sf, coords = ~ exposure,
                      control = ctl), "Do not supply")
  longlat <- sf::st_transform(dat_sf, 4326)
  expect_error(sp_fit(response ~ exposure, longlat,
                      control = sp_control()), "must be projected")
  wrong_crs <- sf::st_transform(dat_sf[1:2, ], 32630)
  expect_error(predict(fit_sf, wrong_crs), "same coordinate reference")
}

# The independent limiting model agrees with Gaussian linear-model algebra.
independent <- sp_fit(response ~ exposure, dat, spatial = FALSE)
lm_fit <- lm(response ~ exposure, dat)
stopifnot(all.equal(unname(coef(independent)), unname(coef(lm_fit)),
                    tolerance = 1e-12))
independent_sigma2 <- sum(residuals(lm_fit)^2) / nrow(dat)
stopifnot(all.equal(independent$covariance[["nugget_variance"]],
                    independent_sigma2, tolerance = 1e-12))
stopifnot(identical(attr(logLik(independent), "df"), 3L))
stopifnot(all.equal(as.numeric(logLik(independent)),
                    as.numeric(stats::logLik(lm_fit)), tolerance = 1e-12))
independent_with_coords <- sp_fit(
  response ~ exposure, dat, spatial = FALSE, coords = ~ east + north,
  control = sp_control(distance_units = "km")
)
prefit_variogram <- sp_variogram(
  lm_fit, coords = ~ east + north, data = dat, bins = 3
)
nonspatial_variogram <- sp_variogram(independent_with_coords, bins = 3)
direct_variogram <- sp_variogram(
  coords = as.matrix(dat[c("east", "north")]),
  values = residuals(lm_fit), bins = 3, distance_units = "km"
)
stopifnot(identical(prefit_variogram$stage, "pre_fit"),
          identical(prefit_variogram$residual_type, "response"),
          identical(prefit_variogram$values$semivariance,
                    direct_variogram$values$semivariance),
          isTRUE(all.equal(nonspatial_variogram$values$semivariance,
                           prefit_variogram$values$semivariance,
                           tolerance = 1e-12)),
          all(is.na(prefit_variogram$values$fitted)))
conditional_variogram <- sp_variogram(fit, residual_type = "conditional",
                                      bins = 3)
stopifnot(identical(conditional_variogram$stage, "post_fit"))

beta_standard_error_expected <- sqrt(diag(
  solve(crossprod(x, v_inverse %*% x))
))
stopifnot(all.equal(unname(sqrt(diag(vcov(fit)))),
                    unname(beta_standard_error_expected), tolerance = 1e-12))
independent_prediction <- predict(independent, newdata = new, se.fit = TRUE,
                                  target = "observation")
stopifnot(all(is.finite(independent_prediction$fit)),
          all(independent_prediction$se.fit > 0))
expect_error(sp_variogram(independent), "Supply two coordinate columns")

cat("Gaussian fitting and prediction tests passed.\n")
