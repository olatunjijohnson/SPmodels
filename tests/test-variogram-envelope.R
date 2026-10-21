if (file.exists(file.path("R", "assessment.R"))) {
  for (file in list.files("R", full.names = TRUE)) source(file)
} else {
  library(SPmodels)
}

close <- function(a, b, tolerance = 1e-12) {
  isTRUE(all.equal(a, b, tolerance = tolerance, check.attributes = FALSE))
}
expect_error <- function(expression, pattern) {
  error <- tryCatch({ force(expression); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}

dat <- data.frame(
  y = c(1.0, 1.5, 0.8, 2.1, 1.7, 2.5, 2.2, 3.0, 2.7, 3.4, 3.1, 3.8),
  x = seq(-1, 1.2, length.out = 12),
  east = rep(0:3, 3), north = rep(0:2, each = 4)
)
fit <- sp_fit(
  y ~ x, dat, coords = ~ east + north,
  control = sp_control(
    spatial_variance = 0.8, scale = 1.3, nugget_variance = 0.25,
    estimate = FALSE, distance_units = "km"
  )
)

# Reconstruct every seeded permutation curve and the global statistic directly.
set.seed(450)
seed_before <- .Random.seed
envelope <- sp_variogram_envelope(
  fit, nsim = 39, method = "permutation", level = 0.9,
  bins = 4, residual_type = "mean", seed = 73
)
stopifnot(identical(.Random.seed, seed_before),
          identical(dim(envelope$simulations), c(4L, 39L)),
          envelope$p_value >= 1 / 40, envelope$p_value <= 1,
          all(envelope$values$lower <= envelope$values$upper,
              na.rm = TRUE))

distance <- as.matrix(dist(fit$coordinates))
triangle <- lower.tri(distance)
pair_distance <- distance[triangle]
breaks <- seq(0, envelope$max_distance, length.out = 5)
group <- cut(pair_distance, breaks = breaks, include.lowest = TRUE,
             right = TRUE, labels = FALSE)
manual_curve <- function(values) {
  semivariance <- 0.5 * outer(values, values, "-")[triangle]^2
  grouped <- tapply(semivariance, group, mean)
  answer <- rep(NA_real_, 4)
  answer[as.integer(names(grouped))] <- grouped
  answer
}
set.seed(73)
manual_simulations <- replicate(39, manual_curve(sample(fit$residual_mean)))
stopifnot(close(envelope$simulations, manual_simulations))
manual_centre <- rowMeans(manual_simulations)
manual_spread <- apply(manual_simulations, 1, sd)
manual_observed <- envelope$values$semivariance
usable <- is.finite(manual_observed) & manual_spread > 0
manual_observed_statistic <- sum(
  ((manual_observed[usable] - manual_centre[usable]) /
     manual_spread[usable])^2
)
manual_simulated_statistic <- apply(manual_simulations[usable, , drop = FALSE],
  2, function(curve) sum(
    ((curve - manual_centre[usable]) / manual_spread[usable])^2
  ))
manual_p <- (1 + sum(manual_simulated_statistic >=
                       manual_observed_statistic)) / 40
stopifnot(close(envelope$observed_statistic, manual_observed_statistic),
          close(envelope$simulated_statistic,
                manual_simulated_statistic),
          close(envelope$p_value, manual_p))

# Parametric envelopes are locally seeded and support both residual choices.
parametric_one <- sp_variogram_envelope(
  fit, nsim = 29, method = "parametric", level = 0.8,
  bins = 4, residual_type = "conditional", seed = 74
)
parametric_two <- sp_variogram_envelope(
  fit, nsim = 29, method = "parametric", level = 0.8,
  bins = 4, residual_type = "conditional", seed = 74
)
stopifnot(identical(parametric_one$simulations,
                    parametric_two$simulations),
          identical(parametric_one$method, "parametric"),
          is.finite(parametric_one$p_value))

independent <- sp_fit(y ~ x, dat, spatial = FALSE)
independent_with_coords <- sp_fit(
  y ~ x, dat, spatial = FALSE, coords = ~ east + north,
  control = sp_control(distance_units = "km")
)
ordinary <- lm(y ~ x, dat)
ordinary_envelope <- sp_variogram_envelope(
  ordinary, coords = ~ east + north, data = dat, nsim = 39,
  method = "permutation", bins = 4, distance_units = "km", seed = 101
)
nonspatial_envelope <- sp_variogram_envelope(
  independent_with_coords, nsim = 39, method = "permutation", bins = 4,
  seed = 101
)
raw_envelope <- sp_variogram_envelope(
  coords = as.matrix(dat[c("east", "north")]),
  values = residuals(ordinary), nsim = 39, method = "permutation",
  bins = 4, distance_units = "km", seed = 101
)
stopifnot(identical(ordinary_envelope$stage, "pre_fit"),
          identical(ordinary_envelope$method, "permutation"),
          close(ordinary_envelope$simulations,
                nonspatial_envelope$simulations),
          close(ordinary_envelope$simulations, raw_envelope$simulations),
          identical(ordinary_envelope$distance_units, "km"))

binary <- transform(dat, present = as.integer(y > median(y)))
ordinary_glm <- glm(present ~ x, data = binary, family = binomial())
glm_envelope <- sp_variogram_envelope(
  ordinary_glm, coords = ~ east + north, data = binary, nsim = 19,
  method = "permutation", bins = 4, seed = 102
)
stopifnot(identical(glm_envelope$stage, "pre_fit"),
          identical(glm_envelope$family, "binomial"))

if (requireNamespace("sf", quietly = TRUE)) {
  projected <- sf::st_as_sf(
    transform(dat, east = 371000 + east * 1000,
              north = 408000 + north * 1000),
    coords = c("east", "north"), crs = 27700
  )
  ordinary_projected <- lm(y ~ x, data = projected)
  variogram_km <- sp_variogram(
    ordinary_projected, data = projected, distance_units = "km", bins = 4
  )
  variogram_m <- sp_variogram(
    ordinary_projected, data = projected, distance_units = "m", bins = 4
  )
  stopifnot(identical(variogram_km$distance_units, "km"),
            identical(variogram_m$distance_units, "m"),
            close(variogram_m$values$distance,
                  variogram_km$values$distance * 1000))
  spatial_variogram_km <- sp_variogram(fit, bins = 4)
  spatial_variogram_m <- sp_variogram(fit, distance_units = "m", bins = 4)
  stopifnot(close(spatial_variogram_m$values$distance,
                  spatial_variogram_km$values$distance * 1000))
}

expect_error(sp_variogram_envelope(independent, nsim = 9),
             "Supply two coordinate columns")
expect_error(sp_variogram_envelope(independent_with_coords,
                                   method = "parametric", nsim = 9),
             "requires a fitted spatial model")
expect_error(sp_variogram_envelope(fit, nsim = 0), "positive integer")
expect_error(sp_variogram_envelope(fit, nsim = 9, level = 1),
             "strictly between")

cat("Variogram-envelope tests passed.\n")
