if (file.exists(file.path("R", "neighbor-prediction.R"))) {
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
  y = c(1.1, 1.8, 1.3, 2.4, 2.1, 3.0),
  z = c(-1, -.5, 0, .5, 1, 1.5),
  east = c(0, 1, 0, 2, 1, 3), north = c(0, 0, 1, 1, 2, 2),
  village = rep(c("a", "b", "c"), each = 2)
)
control <- sp_control(
  spatial_variance = .8, scale = 1.4, nugget_variance = .25,
  random_variance = .3, estimate = FALSE, distance_units = "km"
)
fit <- sp_fit(y ~ z, dat, coords = ~ east + north, random = ~ village,
              control = control)

# Retaining every predecessor is the ordinary covariance factorization written
# sequentially. With the same standard normals it reproduces exact joint draws.
exact <- sp_predict_draws(fit, dat, nsim = 40L, target = "observation",
                          type = "link", seed = 612L)
sequential_exact <- sp_predict_draws_neighbor(
  fit, dat, nsim = 40L, neighbors = nrow(dat) - 1L,
  target = "observation", type = "link", seed = 612L
)
stopifnot(close(sequential_exact$draws, exact$draws, 2e-8),
          close(sequential_exact$predictive_mean, exact$predictive_mean,
                1e-10),
          is.null(sequential_exact$predictive_covariance),
          identical(sequential_exact$prediction_order, seq_len(nrow(dat))),
          identical(sequential_exact$conditional_neighbors[[6L]], 1:5))
exact_covariance_diagnostic <- sp_neighbor_covariance(
  fit, dat, neighbors = nrow(dat), target = "observation", reference = TRUE
)
stopifnot(close(exact_covariance_diagnostic$covariance,
                exact$predictive_covariance, 2e-10),
          exact_covariance_diagnostic$errors[["relative_frobenius"]] < 1e-10,
          any(grepl("predictive covariance diagnostic",
                    capture.output(print(exact_covariance_diagnostic)),
                    fixed = TRUE)))

# The optional fitted-parameter layer uses the same parameter draws and random
# ordering as the ordinary exact route. Its all-predecessor boundary therefore
# agrees draw for draw while still omitting a stored full covariance matrix.
parameter_exact <- sp_predict_draws(
  fit, dat, nsim = 20L, target = "latent", type = "link", seed = 611L,
  parameter_uncertainty = TRUE
)
parameter_neighbor <- sp_predict_draws_neighbor(
  fit, dat, nsim = 20L, neighbors = nrow(dat), target = "latent",
  type = "link", seed = 611L, parameter_uncertainty = TRUE
)
stopifnot(close(parameter_neighbor$parameter_draws,
                parameter_exact$parameter_draws, 0),
          close(parameter_neighbor$draws, parameter_exact$draws, 2e-8),
          parameter_neighbor$parameter_uncertainty,
          all(is.finite(parameter_neighbor$conditional_variance)))

# A reordered full-predecessor factorization still defines the exact joint
# Gaussian distribution and returns columns in original row order.
reverse_order <- nrow(dat):1L
reordered <- sp_predict_draws_neighbor(
  fit, dat, nsim = 6000L, neighbors = nrow(dat),
  prediction_order = reverse_order, target = "latent", type = "link",
  seed = 613L
)
direct <- sp_predict_draws(fit, dat, nsim = 6000L, target = "latent",
                           type = "link", seed = 614L)
stopifnot(max(abs(colMeans(reordered$draws) - exact$predictive_mean)) < .06,
          max(abs(stats::cov(reordered$draws) -
                  direct$predictive_covariance)) < .07,
          identical(reordered$prediction_order, reverse_order))

# A small predecessor budget remains reproducible, records the local graph,
# avoids a full covariance result and preserves the caller's RNG state.
set.seed(614L)
before <- .Random.seed
local_one <- sp_predict_draws_neighbor(
  fit, dat, nsim = 50L, neighbors = 2L, seed = 99L, type = "link"
)
stopifnot(identical(.Random.seed, before))
local_two <- sp_predict_draws_neighbor(
  fit, dat, nsim = 50L, neighbors = 2L, seed = 99L, type = "link"
)
stopifnot(close(local_one$draws, local_two$draws, 0),
          all(lengths(local_one$conditional_neighbors) <= 2L),
          identical(local_one$approximation,
                    "nearest-predecessor Gaussian"),
          any(grepl("nearest-predecessor", capture.output(print(local_one)),
                    fixed = TRUE)))
local_covariance_diagnostic <- sp_neighbor_covariance(
  fit, dat, neighbors = 2L, reference = TRUE
)
stopifnot(all(is.finite(local_covariance_diagnostic$errors)),
          min(eigen(local_covariance_diagnostic$covariance, symmetric = TRUE,
                    only.values = TRUE)$values) > -1e-10,
          all(lengths(local_covariance_diagnostic$conditional_neighbors) <= 2L))
local_target <- sp_predict_target(
  local_one, weights = rep(1, ncol(local_one$draws)), threshold = 2
)
stopifnot(close(local_target$draws[, 1L], rowMeans(local_one$draws), 1e-14))

# Direct target streaming reproduces post-hoc aggregation for Gaussian draws
# while retaining only graph-active point vectors.
target_weights <- c(1, 2, 1, 3, 2, 1)
target_groups <- rep(c("west", "east"), each = 3L)
posthoc_target <- sp_predict_target(
  local_one, weights = target_weights, groups = target_groups,
  threshold = 2
)
set.seed(619L)
target_seed_before <- .Random.seed
streamed_target <- sp_predict_target_neighbor(
  fit, dat, weights = target_weights, groups = target_groups,
  nsim = 50L, neighbors = 2L, seed = 99L, threshold = 2
)
streamed_table <- sp_prediction_table(streamed_target)
stopifnot(identical(.Random.seed, target_seed_before),
          close(streamed_target$draws, posthoc_target$draws, 1e-13),
          close(streamed_target$summary$mean,
                posthoc_target$summary$mean, 1e-13),
          streamed_target$peak_active_locations <= nrow(dat),
          ncol(streamed_target$draws) == 2L,
          nrow(streamed_target$weights) == nrow(dat),
          close(streamed_table$mean, streamed_target$summary$mean, 0),
          any(grepl("Peak active", capture.output(print(streamed_target)),
                    fixed = TRUE)))

# Deterministic maximin ordering starts at row one, then repeatedly chooses the
# farthest row from the selected set without constructing a full distance
# matrix. An explicit permutation takes precedence and is labelled as user.
maximin <- sp_predict_draws_neighbor(
  fit, dat, nsim = 20L, neighbors = 2L, ordering = "maximin", seed = 98L,
  type = "link"
)
stopifnot(identical(maximin$prediction_order, c(1L, 6L, 5L, 4L, 2L, 3L)),
          identical(maximin$ordering, "maximin"),
          identical(reordered$ordering, "user"))

# Attributes identifying a shared unseen grouping level survive local design
# subsetting, so the full-predecessor boundary also matches there.
new_group <- dat[1:3, ]
new_group$village <- c("new", "new", "a")
group_exact <- sp_predict_draws(
  fit, new_group, nsim = 30L, target = "latent", type = "link",
  unseen = "marginal", seed = 617L
)
group_neighbor <- sp_predict_draws_neighbor(
  fit, new_group, nsim = 30L, neighbors = 3L,
  target = "latent", type = "link", unseen = "marginal", seed = 617L
)
stopifnot(close(group_neighbor$draws, group_exact$draws, 2e-8))

# Coincident latent prediction rows define a valid singular Gaussian
# distribution. The sequential route uses a checked PSD solve, giving the
# duplicate rows identical draws instead of requiring artificial jitter.
duplicated_new <- dat[c(1L, 1L, 2L), ]
row.names(duplicated_new) <- c("first", "duplicate", "other")
duplicated_draws <- sp_predict_draws_neighbor(
  fit, duplicated_new, nsim = 30L, neighbors = 3L,
  target = "latent", type = "link", seed = 618L
)
stopifnot(close(duplicated_draws$draws[, "first"],
                duplicated_draws$draws[, "duplicate"], 1e-10),
          duplicated_draws$conditional_variance[["duplicate"]] < 1e-12)
duplicated_covariance <- sp_neighbor_covariance(
  fit, duplicated_new, neighbors = 3L, target = "latent", reference = TRUE
)
stopifnot(close(duplicated_covariance$covariance["first", ],
                duplicated_covariance$covariance["duplicate", ], 1e-10))

# Non-Gaussian link draws use the same local Gaussian predictive approximation;
# latent response transformation respects support.
binary <- transform(dat, infected = c(0, 1, 0, 1, 1, 1))
bin_fit <- sp_fit(infected ~ z, binary, family = binomial(),
  coords = ~ east + north,
  control = sp_control(spatial_variance = .7, scale = 1.2,
    nugget_variance = .2, estimate = FALSE, distance_units = "km"))
bin_exact <- sp_predict_draws(bin_fit, binary, nsim = 30L,
  target = "latent", type = "link", seed = 615L)
bin_neighbor <- sp_predict_draws_neighbor(bin_fit, binary, nsim = 30L,
  neighbors = nrow(binary), target = "latent", type = "link", seed = 615L)
bin_response <- sp_predict_draws_neighbor(bin_fit, binary, nsim = 30L,
  neighbors = 2L, target = "latent", type = "response", seed = 616L)
stopifnot(close(bin_neighbor$draws, bin_exact$draws, 2e-8),
          all(bin_response$draws >= 0 & bin_response$draws <= 1))
bin_covariance <- sp_neighbor_covariance(
  bin_fit, binary, neighbors = nrow(binary), reference = TRUE
)
stopifnot(bin_covariance$errors[["relative_frobenius"]] < 1e-10)

# Repeated polygon allocation indices (as produced by fractional-area
# overlays) stream into several targets without a cell-by-target dense matrix.
polygon_allocation <- structure(list(
  index = c(1L, 1L, 2L, 3L, 4L, 5L, 6L),
  weights = c(.4, .6, 1, 1, 1, 1, 1),
  groups = c("A", "B", "A", "A", "B", "B", "B"),
  point_count = nrow(binary)
), class = "spmodels_polygon_weights")
polygon_posthoc <- sp_predict_target(bin_response, polygon_allocation)
polygon_stream <- sp_predict_target_neighbor(
  bin_fit, binary, weights = polygon_allocation, nsim = 30L,
  neighbors = 2L, target = "latent", type = "response", seed = 616L
)
stopifnot(close(polygon_stream$draws, polygon_posthoc$draws, 1e-13),
          nrow(polygon_stream$weights) == length(polygon_allocation$index))

expect_error(sp_predict_draws_neighbor(lm(y ~ z, dat), dat), "spatial")
expect_error(sp_predict_draws_neighbor(fit, dat, neighbors = 0), "positive")
expect_error(sp_predict_draws_neighbor(fit, dat, ordering = "random"),
             "arg")
expect_error(sp_predict_draws_neighbor(
  fit, dat, parameter_uncertainty = NA
), "TRUE or FALSE")
expect_error(sp_predict_draws_neighbor(
  fit, dat, prediction_order = rep(1, nrow(dat))
), "permutation")
expect_error(sp_predict_target_neighbor(
  fit, dat, weights = rep(0, nrow(dat))
), "positive total")
expect_error(sp_predict_target_neighbor(
  fit, dat, weights = 1:2
), "matching prediction")
expect_error(sp_neighbor_covariance(fit, dat, reference = NA),
             "TRUE or FALSE")

cat("Nearest-predecessor predictive-draw tests passed.\n")

