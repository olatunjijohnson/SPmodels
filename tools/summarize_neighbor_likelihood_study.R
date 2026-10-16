args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1L) args[[1L]] else
  "docs/benchmarks/data/neighbor_likelihood_study.csv"
output_path <- if (length(args) >= 2L) args[[2L]] else
  "docs/benchmarks/data/neighbor_likelihood_summary.csv"

data <- utils::read.csv(input_path, stringsAsFactors = FALSE)
required <- c(
  "n", "range_fraction", "ordering", "neighbors", "exact_boundary",
  "coefficient_error_in_dense_se", "prediction_rmse_ratio",
  "log_likelihood_error", "retained_object_bytes"
)
if (!all(required %in% names(data))) {
  stop("The input file does not contain the required study fields.",
       call. = FALSE)
}
data$passes <- data$coefficient_error_in_dense_se < 0.05 &
  abs(log(data$prediction_rmse_ratio)) <= log(1.02)
groups <- split(data, interaction(data$n, data$range_fraction,
                                  data$ordering, data$neighbors, drop = TRUE))
summary <- do.call(rbind, lapply(groups, function(rows) {
  result <- data.frame(
    n = rows$n[[1L]], range_fraction = rows$range_fraction[[1L]],
    ordering = rows$ordering[[1L]], neighbors = rows$neighbors[[1L]],
    exact_boundary = all(rows$exact_boundary),
    pass_fraction = mean(rows$passes),
    median_coefficient_error_in_dense_se =
      median(rows$coefficient_error_in_dense_se),
    median_absolute_log_likelihood_error =
      median(abs(rows$log_likelihood_error)),
    median_absolute_log_prediction_rmse_ratio =
      median(abs(log(rows$prediction_rmse_ratio))),
    median_retained_object_kb = median(rows$retained_object_bytes) / 1024
  )
  if ("fit_elapsed_seconds" %in% names(rows)) {
    result$median_fit_elapsed_seconds <- median(rows$fit_elapsed_seconds)
  }
  result
}))
summary <- summary[order(summary$n, summary$range_fraction,
                         summary$ordering, summary$neighbors), ]
utils::write.csv(summary, output_path, row.names = FALSE)

scenarios <- split(summary, interaction(summary$n, summary$range_fraction,
                                        drop = TRUE))
selected <- do.call(rbind, lapply(scenarios, function(rows) {
  eligible <- rows[rows$pass_fraction >= 0.90, ]
  selected_budget <- NA_integer_
  for (budget in sort(unique(eligible$neighbors))) {
    both_orders <- eligible[eligible$neighbors == budget, ]
    if (all(c("input", "maximin") %in% both_orders$ordering)) {
      selected_budget <- as.integer(budget)
      break
    }
  }
  data.frame(
    n = rows$n[[1L]], range_fraction = rows$range_fraction[[1L]],
    selected_budget = selected_budget,
    status = if (is.na(selected_budget)) "inconclusive" else
      if (selected_budget == rows$n[[1L]] - 1L) "dense boundary required" else
        "scenario-specific budget passed"
  )
}))
cat("Rows:", nrow(data), "\n")
cat("Unique response datasets:", length(unique(data$response_seed)), "\n")
cat("Archive checksums:", paste(unique(data$archive_sha256), collapse = ", "), "\n")
print(selected, row.names = FALSE)
