args <- commandArgs(trailingOnly = TRUE)
replicates <- if (length(args) >= 1L) as.integer(args[[1L]]) else 25L
sizes <- if (length(args) >= 2L) {
  as.integer(strsplit(args[[2L]], ",", fixed = TRUE)[[1L]])
} else c(100L, 225L, 400L)
archive_sha256 <- if (length(args) >= 3L) args[[3L]] else NA_character_
output_path <- if (length(args) >= 4L) args[[4L]] else
  "docs/benchmarks/data/neighbor_likelihood_study.csv"

if (is.na(replicates) || replicates < 1L || anyNA(sizes) ||
    any(sizes < 25L) || !length(sizes)) {
  stop("Supply positive replicates and sizes of at least 25.", call. = FALSE)
}
for (file in list.files("R", full.names = TRUE)) source(file)

smoothness <- 0.5
spatial_variance <- 0.8
nugget_variance <- 0.2
covariate_effect <- 0.4
intercept <- 0.3
range_fractions <- c(0.05, 0.20, 0.50)
budgets <- c(5L, 10L, 20L, 40L)
scenario_index <- -1L
results <- list()
result_index <- 0L

for (n in sizes) {
  for (range_fraction in range_fractions) {
    scenario_index <- scenario_index + 1L
    scenario_seed <- 460001L + scenario_index * 10000L
    for (replicate in seq_len(replicates)) {
      response_seed <- scenario_seed + replicate
      set.seed(response_seed)
      total_n <- n + 50L
      coordinates <- cbind(east = runif(total_n), north = runif(total_n))
      scale <- range_fraction * sqrt(2)
      distance <- .sp_dist(coordinates)
      covariance <- spatial_variance *
        .sp_matern_cor_unchecked(distance, scale, smoothness)
      diag(covariance) <- diag(covariance) + nugget_variance
      covariate <- rnorm(total_n)
      response <- intercept + covariate_effect * covariate +
        as.vector(t(chol(covariance)) %*% rnorm(total_n))
      data <- data.frame(
        y = response[seq_len(n)], z = covariate[seq_len(n)],
        east = coordinates[seq_len(n), 1L],
        north = coordinates[seq_len(n), 2L]
      )
      test <- data.frame(
        z = covariate[n + seq_len(50L)], east = coordinates[n +
          seq_len(50L), 1L], north = coordinates[n + seq_len(50L), 2L]
      )
      truth <- response[n + seq_len(50L)]
      control <- sp_control(
        smoothness = smoothness, spatial_variance = spatial_variance,
        scale = scale, nugget_variance = nugget_variance,
        estimate = FALSE, distance_units = "unit"
      )
      dense <- sp_fit(y ~ z, data, coords = ~ east + north,
                      control = control)
      dense_prediction <- predict(dense, test)
      dense_rmse <- sqrt(mean((dense_prediction - truth)^2))
      configuration <- expand.grid(
        ordering = c("input", "maximin"),
        budget = unique(c(budgets[budgets < n], n - 1L)),
        stringsAsFactors = FALSE
      )
      for (row in seq_len(nrow(configuration))) {
        ordering <- configuration$ordering[[row]]
        budget <- configuration$budget[[row]]
        fit_elapsed <- system.time(fit <- sp_fit_neighbor(
          y ~ z, data, coords = ~ east + north,
          neighbors = budget, ordering = ordering, control = control
        ))[["elapsed"]]
        prediction <- predict(fit, test)
        result_index <- result_index + 1L
        results[[result_index]] <- data.frame(
          archive_sha256 = archive_sha256,
          r_version = as.character(getRversion()),
          platform = R.version$platform,
          n = n, range_fraction = range_fraction,
          range_scale = scale, replicate = replicate,
          response_seed = response_seed,
          ordering = ordering, neighbors = budget,
          exact_boundary = budget >= n - 1L,
          fit_elapsed_seconds = unname(fit_elapsed),
          dense_log_likelihood = dense$log_likelihood,
          neighbor_log_likelihood = as.numeric(logLik(fit)),
          log_likelihood_error = as.numeric(logLik(fit)) - dense$log_likelihood,
          coefficient_error_intercept = unname(
            coef(fit)[["(Intercept)"]] - dense$coefficients[["(Intercept)"]]
          ),
          coefficient_error_slope = unname(
            coef(fit)[["z"]] - dense$coefficients[["z"]]
          ),
          coefficient_error_in_dense_se = max(abs(
            (coef(fit) - dense$coefficients) /
              sqrt(diag(dense$coefficient_vcov))
          )),
          dense_prediction_rmse = dense_rmse,
          neighbor_prediction_rmse = sqrt(mean((prediction - truth)^2)),
          prediction_rmse_ratio = sqrt(mean((prediction - truth)^2)) /
            dense_rmse,
          retained_object_bytes = as.numeric(object.size(fit)),
          stringsAsFactors = FALSE
        )
      }
      cat("completed n=", n, " range_fraction=", range_fraction,
          " replicate=", replicate, "\n", sep = "")
    }
  }
}

answer <- do.call(rbind, results)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(answer, output_path, row.names = FALSE, na = "")
cat("wrote ", nrow(answer), " rows to ", output_path, "\n", sep = "")
