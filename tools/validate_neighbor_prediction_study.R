args <- commandArgs(trailingOnly = TRUE)
input_paths <- if (length(args) >= 1L) args else
  "docs/validation/data/neighbor_prediction_gaussian.csv"
data <- do.call(rbind, lapply(input_paths, function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, na.strings = "")
}))
required <- c("archive_sha256", "case", "grid_size", "replicate",
  "response_seed", "nsim", "ordering", "neighbors", "exact_boundary",
  "exact_boundary_same_seed_equal", "target_id", "error",
  "relative_frobenius_covariance_error",
  "maximum_absolute_marginal_variance_error", "exact_target_sd",
  "approximate_target_sd", "threshold_probability_difference")
stopifnot(nrow(data) > 0L, all(required %in% names(data)))
stopifnot(length(unique(data$archive_sha256)) == 1L,
          grepl("^[[:xdigit:]]{64}$", unique(data$archive_sha256)))
stopifnot(all(data$nsim == 20000L),
          all(data$grid_size %in% c(25L, 100L, 225L)),
          all(data$ordering %in% c("input", "maximin")))
stopifnot(!anyDuplicated(data[c("case", "grid_size", "replicate",
  "ordering", "neighbors", "target_id")]))

target_ids <- c("overall", "west_south", "west_north", "east_south",
                "east_north")
for (case in unique(data$case)) {
  rows <- data[data$case == case, , drop = FALSE]
  stopifnot(nrow(rows) == 600L,
            setequal(unique(rows$replicate), 1:5),
            all(rows$response_seed == 430000L + rows$replicate))
  for (grid_size in c(25L, 100L, 225L)) {
    grid <- rows[rows$grid_size == grid_size, , drop = FALSE]
    expected_budgets <- unique(pmin(c(5L, 10L, 20L, 40L), grid_size - 1L))
    stopifnot(setequal(unique(grid$neighbors), expected_budgets),
      nrow(grid) == 200L,
      setequal(unique(grid$target_id), target_ids),
      all(table(grid$replicate, grid$ordering, grid$neighbors,
                grid$target_id) == 1L))
  }
}

boundary <- data$exact_boundary
input_boundary <- boundary & data$ordering == "input"
stopifnot(all(data$exact_boundary_same_seed_equal[input_boundary]),
          all(is.na(data$exact_boundary_same_seed_equal[
            boundary & data$ordering == "maximin"])),
          all(data$relative_frobenius_covariance_error[input_boundary] <
                2e-10),
          all(data$maximum_absolute_marginal_variance_error[
            input_boundary] < 2e-10))
successful <- is.na(data$error) | !nzchar(data$error)
stopifnot(all(is.finite(data$relative_frobenius_covariance_error[successful])),
          all(is.finite(data$exact_target_sd[successful])),
          all(is.finite(data$approximate_target_sd[successful])),
          all(is.finite(data$threshold_probability_difference[successful])))
cat("NNP raw design, seeds, retained rows and exact-boundary checks passed; ",
    sum(!successful), " failed configurations retained.\n", sep = "")
