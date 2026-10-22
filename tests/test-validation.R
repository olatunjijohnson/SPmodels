if (file.exists(file.path("R", "validation.R"))) {
  for (file in list.files("R", full.names = TRUE)) source(file)
} else {
  library(SPmodels)
}

close <- function(a, b, tolerance = 1e-9) {
  isTRUE(all.equal(a, b, tolerance = tolerance, check.attributes = FALSE))
}

dat <- data.frame(
  y = c(1.0, 1.5, 0.8, 2.1, 1.7, 2.5, 2.2, 3.0, 2.7, 3.4, 3.1, 3.8),
  x = seq(-1, 1.2, length.out = 12),
  east = rep(0:3, 3), north = rep(0:2, each = 4)
)
control <- sp_control(spatial_variance = 0.8, scale = 1.3,
                      nugget_variance = 0.25, estimate = FALSE,
                      distance_units = "km")
fit <- sp_fit(y ~ x, dat, coords = ~ east + north, control = control)
folds <- rep(1:3, each = 4)
cv <- sp_cross_validate(fit, folds)
stopifnot(identical(cv$folds, folds), nrow(cv$predictions) == nrow(dat),
          cv$summary$n == nrow(dat), is.finite(cv$summary$rmse),
          is.finite(cv$summary$mean_log_score))

# One fold agrees with an independently performed refit and prediction.
manual_fit <- sp_fit(y ~ x, dat[folds != 1, ], coords = ~ east + north,
                     control = control)
manual_prediction <- predict(manual_fit, dat[folds == 1, ], se.fit = TRUE,
                             target = "observation")
stopifnot(close(cv$predictions$predicted[folds == 1], manual_prediction$fit),
          close(cv$predictions$se.fit[folds == 1], manual_prediction$se.fit))
gaussian_score <- sp_predictive_log_score(
  manual_fit, dat[folds == 1, ], dat$y[folds == 1]
)
stopifnot(close(gaussian_score,
  -dnorm(dat$y[folds == 1], manual_prediction$fit,
         manual_prediction$se.fit, log = TRUE)))

# Generated balanced folds are seed-local and reproducible.
set.seed(88)
seed_before <- .Random.seed
random_cv_one <- sp_cross_validate(fit, folds = 3, seed = 19)
stopifnot(identical(.Random.seed, seed_before))
random_cv_two <- sp_cross_validate(fit, folds = 3, seed = 19)
stopifnot(identical(random_cv_one$folds, random_cv_two$folds),
          all(tabulate(random_cv_one$folds) == 4L))

# Geometry-driven folds are seed-local, use canonical labels, and never split
# repeated coordinates across training and held-out sets.
duplicate_coordinates <- rbind(fit$coordinates, fit$coordinates[c(1, 4), ])
fold_object <- structure(list(
  spatial = TRUE, coordinates = duplicate_coordinates
), class = "spmodels_fit")
set.seed(992)
spatial_seed_before <- .Random.seed
spatial_folds_one <- sp_spatial_folds(fold_object, folds = 3, seed = 27)
stopifnot(identical(.Random.seed, spatial_seed_before))
spatial_folds_two <- sp_spatial_folds(fold_object, folds = 3, seed = 27)
stopifnot(identical(as.integer(spatial_folds_one),
                    as.integer(spatial_folds_two)),
          spatial_folds_one[1] == spatial_folds_one[13],
          spatial_folds_one[4] == spatial_folds_one[14],
          identical(sort(unique(as.integer(spatial_folds_one))), 1:3),
          identical(dim(attr(spatial_folds_one, "centres")), c(3L, 2L)))
fit_spatial_folds <- sp_spatial_folds(fit, folds = 3, seed = 28)
spatial_cv <- sp_cross_validate(fit, fit_spatial_folds)
stopifnot(all(is.finite(spatial_cv$predictions$predicted)))

# Regular blocks retain whole occupied cells and expose an independently
# reconstructable canonical allocation.
block_folds <- sp_block_folds(
  fit, block_size = c(2, 2), folds = 2, assignment = "systematic"
)
block_index <- floor(fit$coordinates /
  rep(c(2, 2), each = nrow(fit$coordinates)))
block_key <- paste(block_index[, 1], block_index[, 2], sep = "\r")
unique_block <- !duplicated(block_key)
block_coordinates <- block_index[unique_block, , drop = FALSE]
canonical_block <- order(-block_coordinates[, 2], block_coordinates[, 1])
block_assignment <- integer(nrow(block_coordinates))
block_assignment[canonical_block] <- rep(1:2,
  length.out = nrow(block_coordinates))
manual_block_folds <- block_assignment[match(
  block_key, block_key[unique_block]
)]
set.seed(281)
block_seed_before <- .Random.seed
random_blocks_one <- sp_block_folds(
  fit, 2, folds = 2, assignment = "random", seed = 29
)
random_blocks_two <- sp_block_folds(
  fit, 2, folds = 2, assignment = "random", seed = 29
)
block_cv <- sp_cross_validate(fit, block_folds)
stopifnot(inherits(block_folds, "spmodels_block_folds"),
          identical(as.integer(block_folds), as.integer(manual_block_folds)),
          identical(random_blocks_one, random_blocks_two),
          identical(.Random.seed, block_seed_before),
          all(vapply(split(as.integer(block_folds), block_key),
                     function(value) length(unique(value)) == 1L,
                     logical(1))),
          all(is.finite(block_cv$predictions$predicted)))
stopifnot(inherits(try(sp_block_folds(fit, 100, folds = 2), silent = TRUE),
                   "try-error"),
          inherits(try(sp_block_folds(fit, 0, folds = 2), silent = TRUE),
                   "try-error"))

# Environmental folds cluster unique locations in standardized covariate
# space and leave the caller's RNG state unchanged.
set.seed(282)
environment_seed_before <- .Random.seed
environment_folds_one <- sp_environmental_folds(
  fit, ~ x, folds = 3, seed = 30, nstart = 10
)
environment_folds_two <- sp_environmental_folds(
  fit, "x", folds = 3, seed = 30, nstart = 10
)
environment_cv <- sp_cross_validate(fit, environment_folds_one)
stopifnot(inherits(environment_folds_one, "spmodels_environmental_folds"),
          identical(as.integer(environment_folds_one),
                    as.integer(environment_folds_two)),
          identical(attr(environment_folds_one, "centres"),
                    attr(environment_folds_two, "centres")),
          identical(.Random.seed, environment_seed_before),
          identical(attr(environment_folds_one, "variables"), "x"),
          identical(dim(attr(environment_folds_one, "centres")), c(3L, 1L)),
          all(is.finite(environment_cv$predictions$predicted)),
          inherits(try(sp_environmental_folds(
            fit, "missing", folds = 2
          ), silent = TRUE), "try-error"))

# Group-constrained folds keep every declared unit wholly on one side and
# balance row counts by a directly reconstructable least-loaded allocation.
group_fit <- fit
group_fit$data_used$region <- rep(c("a", "b", "c", "d"), c(4, 3, 3, 2))
group_folds <- sp_group_folds(group_fit, "region", folds = 3)
group_metadata <- attr(group_folds, "groups")
set.seed(283)
group_seed_before <- .Random.seed
random_group_one <- sp_group_folds(
  group_fit, "region", folds = 3, assignment = "random", seed = 31
)
random_group_two <- sp_group_folds(
  group_fit, "region", folds = 3, assignment = "random", seed = 31
)
group_cv <- sp_cross_validate(group_fit, group_folds)
stopifnot(inherits(group_folds, "spmodels_group_folds"),
          identical(group_metadata$rows, c(4L, 3L, 3L, 2L)),
          identical(group_metadata$fold, c(1L, 2L, 3L, 2L)),
          all(vapply(split(as.integer(group_folds), group_fit$data_used$region),
                     function(value) length(unique(value)) == 1L,
                     logical(1))),
          identical(random_group_one, random_group_two),
          identical(.Random.seed, group_seed_before),
          all(is.finite(group_cv$predictions$predicted)))
stopifnot(inherits(try(sp_group_folds(
  group_fit, "region", folds = 5
), silent = TRUE), "try-error"))

# Buffered folds remove every near-test training row independently per fold.
buffered_folds <- sp_buffered_folds(fit, folds = 3, buffer = 1.01, seed = 28)
buffered_cv <- sp_cross_validate(fit, buffered_folds)
buffer_test <- as.integer(buffered_folds) == 1L
buffer_candidates <- which(!buffer_test)
buffer_candidate_coordinates <- fit$coordinates[buffer_candidates, , drop = FALSE]
buffer_test_coordinates <- fit$coordinates[buffer_test, , drop = FALSE]
buffer_distance <- sqrt(pmax(
  outer(rowSums(buffer_candidate_coordinates^2),
        rowSums(buffer_test_coordinates^2), "+") -
    2 * tcrossprod(buffer_candidate_coordinates, buffer_test_coordinates),
  0
))
buffer_train <- buffer_candidates[!apply(buffer_distance <= 1.01, 1L, any)]
buffer_manual_fit <- sp_fit(y ~ x, dat[buffer_train, ],
  coords = ~ east + north, control = control)
buffer_manual_prediction <- predict(
  buffer_manual_fit, dat[buffer_test, ], se.fit = TRUE,
  target = "observation"
)
stopifnot(inherits(buffered_folds, "spmodels_buffered_folds"),
          identical(attr(buffered_folds, "distance_units"), "km"),
          buffered_cv$training_size[["1"]] == length(buffer_train),
          length(buffer_train) < sum(!buffer_test),
          close(buffered_cv$predictions$predicted[buffer_test],
                buffer_manual_prediction$fit),
          close(buffered_cv$predictions$se.fit[buffer_test],
                buffer_manual_prediction$se.fit))

# The same refit contract supports non-Gaussian held-out probabilities.
binomial_data <- transform(dat, y = c(0,1,0,1,1,0,1,0,1,0,1,1))
binomial_fit <- sp_fit(y ~ x, binomial_data, family = binomial(),
  coords = ~ east + north,
  control = sp_control(spatial_variance = 0.7, scale = 1.2,
    nugget_variance = 0.2, estimate = FALSE, distance_units = "km"))
binomial_cv <- sp_cross_validate(binomial_fit, folds)
stopifnot(all(binomial_cv$predictions$predicted > 0 &
              binomial_cv$predictions$predicted < 1),
          all(is.finite(binomial_cv$predictions$log_score)))
score_data <- binomial_data[1:2, ]
score_link <- predict(binomial_fit, score_data, type = "link", se.fit = TRUE,
                      target = "observation")
score_expected <- vapply(1:2, function(i) -log(integrate(function(z) {
  dbinom(binomial_fit$y[i], binomial_fit$trials[i],
         plogis(score_link$fit[i] + score_link$se.fit[i] * z)) * dnorm(z)
}, -Inf, Inf, rel.tol = 1e-8)$value), numeric(1))
stopifnot(close(sp_predictive_log_score(
  binomial_fit, score_data, binomial_fit$y[1:2],
  trials = binomial_fit$trials[1:2]
), score_expected, 1e-9))

# Draw assessment statistics and empirical CRPS match direct calculations.
draw_matrix <- matrix(c(
  0, 1, 2, 3,
  1, 2, 3, 4,
  2, 3, 4, 5
), nrow = 4, ncol = 3, dimnames = list(NULL, c("a", "b", "c")))
draw_object <- structure(list(
  draws = draw_matrix, family = "gaussian", target = "latent",
  type = "response"
), class = "spmodels_prediction_draws")
truth <- c(1.5, 2.5, 3.5)
assessment <- sp_simulation_assess(draw_object, truth, interval = 0.5,
                                   threshold = 2.5)
brute_crps <- function(x, y) mean(abs(x - y)) -
  0.5 * mean(abs(outer(x, x, "-")))
stopifnot(close(assessment$locations$mean, colMeans(draw_matrix)),
          close(assessment$locations$crps,
                vapply(seq_along(truth), function(j)
                  brute_crps(draw_matrix[, j], truth[j]), numeric(1))),
          close(assessment$summary$rmse, 0),
          assessment$summary$coverage == 1,
          all(assessment$locations$exceedance_probability ==
                colMeans(draw_matrix > 2.5)))

# Repeated simulation/refitting reproduces an independently refitted first
# replicate and reports direct parameter-recovery summaries.
set.seed(174)
recovery_seed_before <- .Random.seed
recovery <- sp_simulation_recovery(fit, nsim = 4, seed = 31)
stopifnot(identical(.Random.seed, recovery_seed_before),
          all(is.na(recovery$failures)),
          nrow(recovery$estimates) == 4L,
          nrow(recovery$summary) == length(coef(fit)))
manual_recovery_data <- dat
manual_recovery_data$y <- recovery$simulations[[1L]]
manual_recovery_fit <- sp_fit(y ~ x, manual_recovery_data,
  coords = ~ east + north, control = control)
stopifnot(close(recovery$estimates[1L, ], coef(manual_recovery_fit)),
          close(recovery$summary$bias,
                colMeans(recovery$estimates) - coef(fit)),
          close(recovery$summary$rmse,
                sqrt(colMeans(sweep(recovery$estimates, 2L,
                                    coef(fit), "-")^2))))
set.seed(176)
parallel_recovery_seed_before <- .Random.seed
parallel_recovery <- if (.Platform$OS.type == "windows") recovery else
  sp_simulation_recovery(fit, nsim = 4, seed = 31, workers = 2L)
stopifnot(identical(parallel_recovery$simulations, recovery$simulations),
          identical(parallel_recovery$estimates, recovery$estimates),
          identical(parallel_recovery$standard_errors,
                    recovery$standard_errors),
          identical(parallel_recovery$convergence, recovery$convergence),
          identical(.Random.seed, parallel_recovery_seed_before),
          parallel_recovery$workers == if (.Platform$OS.type == "windows")
            1L else 2L)

# Cartesian scenario construction changes only declared generative truths,
# records their exact order, rebuilds covariance factors, and composes with
# the existing multi-scenario recovery runner.
grid_template <- fit
grid_template$covariance_estimated <- TRUE
grid_template$covariance_estimated_names <- names(grid_template$covariance)
grid_template$control$estimate <- TRUE
simulation_grid <- sp_simulation_grid(
  grid_template, list(`(Intercept)` = c(1.2, 1.8), scale = c(0.8, 1.6)),
  prefix = "design"
)
grid_values <- attr(simulation_grid, "grid")
grid_variance <- simulation_grid[[4L]]$fitted_state$v
manual_variance <- 0.8 * matern_cor(
  sqrt(pmax(outer(rowSums(fit$coordinates^2),
                  rowSums(fit$coordinates^2), "+") -
              2 * tcrossprod(fit$coordinates), 0)),
  scale = 1.6, smoothness = fit$smoothness
)
diag(manual_variance) <- diag(manual_variance) + 0.25
grid_run <- sp_simulation_scenarios(simulation_grid[1:2], nsim = 1, seed = 35)
stopifnot(inherits(simulation_grid, "spmodels_simulation_grid"),
          identical(names(simulation_grid), paste0("design_", 1:4)),
          identical(grid_values$`(Intercept)`, c(1.2, 1.8, 1.2, 1.8)),
          identical(grid_values$scale, c(0.8, 0.8, 1.6, 1.6)),
          close(simulation_grid[[4L]]$coefficients[["(Intercept)"]], 1.8),
          close(simulation_grid[[4L]]$covariance[["scale"]], 1.6),
          close(grid_variance, manual_variance),
          nrow(grid_run$scenarios) == 2L)

# Precision-directed recovery stops on a passing first attempt or evaluates
# the exact declared cap without weakening its success criterion.
set.seed(175)
precision_seed_before <- .Random.seed
precision_pass <- sp_simulation_precision(
  fit, initial_nsim = 3, max_nsim = 6, growth = 2,
  target_mcse = 1e6, min_success = 2, seed = 33
)
precision_cap <- sp_simulation_precision(
  fit, initial_nsim = 2, max_nsim = 3, growth = 2,
  target_mcse = 1e6, min_success = 4, seed = 34
)
precision_parallel <- if (.Platform$OS.type == "windows") precision_pass else
  sp_simulation_precision(
    fit, initial_nsim = 3, max_nsim = 6, growth = 2,
    target_mcse = 1e6, min_success = 2, seed = 33, workers = 2L
  )
manual_precision_mcse <- apply(precision_pass$result$estimates, 2, sd) /
  sqrt(colSums(is.finite(precision_pass$result$estimates)))
stopifnot(identical(.Random.seed, precision_seed_before),
          isTRUE(precision_pass$passed),
          identical(precision_pass$history$replicates, 3L),
          close(precision_pass$mcse, manual_precision_mcse),
          identical(precision_parallel$result$estimates,
                    precision_pass$result$estimates),
          identical(precision_parallel$history, precision_pass$history),
          !precision_cap$passed,
          identical(precision_cap$history$replicates, c(2L, 3L)),
          nrow(precision_cap$result$estimates) == 3L)

# Two-column binomial response reconstruction keeps every simulated success
# count paired with its fitted trial total.
binomial_recovery <- sp_simulation_recovery(binomial_fit, nsim = 2, seed = 32)
stopifnot(all(is.na(binomial_recovery$failures)),
          all(binomial_recovery$simulations >= 0),
          all(binomial_recovery$simulations <= binomial_fit$trials))
binomial_grid <- sp_simulation_grid(
  binomial_fit, list(x = c(-0.25, 0.25)), prefix = "binomial"
)
binomial_grid_simulation <- simulate(binomial_grid[[2L]], nsim = 2, seed = 36)
stopifnot(close(binomial_grid[[2L]]$coefficients[["x"]], 0.25),
          close(crossprod(binomial_grid[[2L]]$fitted_state$chol_sigma),
                binomial_grid[[2L]]$fitted_state$sigma),
          all(binomial_grid_simulation >= 0 & binomial_grid_simulation <= 1))

# Multiple named fitted scenarios preserve each complete experiment, use
# deterministic scenario seeds and expose an unpooled combined summary.
scenario_fit <- sp_fit(y ~ x, transform(dat, y = y + 0.4 * x),
  coords = ~ east + north, control = control)
set.seed(813)
scenario_seed_before <- .Random.seed
scenarios <- sp_simulation_scenarios(
  list(baseline = fit, shifted = scenario_fit),
  nsim = c(shifted = 3L, baseline = 2L), seed = 41
)
stopifnot(identical(.Random.seed, scenario_seed_before),
          identical(names(scenarios$results), c("baseline", "shifted")),
          identical(scenarios$scenarios$replicates, c(2L, 3L)),
          identical(scenarios$scenarios$seed, c(41L, 42L)),
          nrow(scenarios$summary) == 2L * length(coef(fit)),
          identical(unique(scenarios$summary$scenario),
                    c("baseline", "shifted")))
baseline_direct <- sp_simulation_recovery(fit, nsim = 2, seed = 41)
shifted_direct <- sp_simulation_recovery(scenario_fit, nsim = 3, seed = 42)
parallel_scenarios <- if (.Platform$OS.type == "windows") scenarios else
  sp_simulation_scenarios(
    list(baseline = fit, shifted = scenario_fit),
    nsim = c(baseline = 2L, shifted = 3L), seed = 41, workers = 2L
  )
stopifnot(identical(scenarios$results$baseline$estimates,
                    baseline_direct$estimates),
          identical(scenarios$results$shifted$estimates,
                    shifted_direct$estimates),
          identical(parallel_scenarios$results$baseline$estimates,
                    scenarios$results$baseline$estimates),
          identical(parallel_scenarios$results$shifted$estimates,
                    scenarios$results$shifted$estimates),
          nrow(sp_validation_table(scenarios)) > 0L)
expect_error <- function(expression, pattern) {
  error <- tryCatch({ force(expression); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}
expect_error(sp_simulation_scenarios(list(fit, fit), nsim = 1),
             "scenario names")
expect_error(sp_simulation_scenarios(list(a = fit, b = fit),
  nsim = c(a = 1, c = 1)), "match")
unsupported_response_fit <- fit
unsupported_response_fit$formula <- log(y) ~ x
recorded_refit_failure <- sp_simulation_recovery(
  unsupported_response_fit, nsim = 2, seed = 43, on_error = "record"
)
stopifnot(all(grepl("named response", recorded_refit_failure$failures)))
expect_error(sp_simulation_recovery(
  unsupported_response_fit, nsim = 2, seed = 43, on_error = "stop"
), "named response")
expect_error(sp_simulation_grid(grid_template,
  list(nugget_variance = 0)), "positive")
expect_error(sp_simulation_grid(fit, list(scale = 1)), "not recoverable")
ambiguous_template <- grid_template
names(ambiguous_template$coefficients)[2L] <- "scale"
expect_error(sp_simulation_grid(ambiguous_template, list(scale = 1)),
             "ambiguous")
expect_error(sp_buffered_folds(fit, folds = 3, buffer = 0), "positive")
expect_error(sp_cross_validate(fit,
  sp_buffered_folds(fit, folds = 3, buffer = 100, seed = 28)),
  "no training rows")

cat("Cross-validation and simulation-assessment tests passed.\n")
