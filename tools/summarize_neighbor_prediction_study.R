args <- commandArgs(trailingOnly = TRUE)
prefix_arg <- grep("^--output-prefix=", args, value = TRUE)
output_prefix <- if (length(prefix_arg)) {
  sub("^--output-prefix=", "", prefix_arg[[1L]])
} else NULL
input_paths <- args[!startsWith(args, "--output-prefix=")]
if (!length(input_paths)) {
  input_paths <- "docs/validation/data/neighbor_prediction_gaussian.csv"
}
if (is.null(output_prefix)) {
  output_prefix <- sub("\\.csv$", "", input_paths[[1L]])
}
target_output <- paste0(output_prefix, "_target_summary.csv")
recommendation_output <- paste0(output_prefix, "_recommendations.csv")
order_output <- paste0(output_prefix, "_order_sensitivity.csv")

data <- do.call(rbind, lapply(input_paths, function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, na.strings = "")
}))
if (!nrow(data)) stop("No NNP results were supplied.", call. = FALSE)
required <- c("case", "grid_size", "replicate", "response_seed",
  "ordering", "neighbors", "target_id", "exact_target_sd",
  "approximate_target_sd", "threshold_probability_difference",
  "error")
if (!all(required %in% names(data))) {
  stop("NNP input is missing required result columns.", call. = FALSE)
}
if (anyDuplicated(data[c("case", "grid_size", "replicate", "ordering",
                         "neighbors", "target_id")])) {
  stop("NNP inputs contain duplicate configuration/target rows.",
       call. = FALSE)
}

keys <- c("case", "grid_size", "ordering", "neighbors", "target_id")
groups <- split(data, do.call(interaction,
  c(unname(as.list(data[keys])), list(drop = TRUE, lex.order = TRUE))))
summaries <- lapply(groups, function(rows) {
  usable <- (is.na(rows$error) | !nzchar(rows$error)) &
    is.finite(rows$exact_target_sd) & rows$exact_target_sd > 0 &
    is.finite(rows$approximate_target_sd) &
    is.finite(rows$threshold_probability_difference)
  differences <- rows$threshold_probability_difference[usable]
  n <- length(differences)
  mean_difference <- if (n) mean(differences) else NA_real_
  mcse <- if (n > 1L) stats::sd(differences) / sqrt(n) else NA_real_
  critical <- if (n > 1L) stats::qt(0.975, df = n - 1L) else NA_real_
  lower <- mean_difference - critical * mcse
  upper <- mean_difference + critical * mcse
  interval_status <- if (n < 5L || !is.finite(lower) || !is.finite(upper)) {
    "inconclusive"
  } else if (lower >= -0.02 && upper <= 0.02) {
    "pass"
  } else if (lower > 0.02 || upper < -0.02) {
    "fail"
  } else "inconclusive"
  relative_sd_error <- if (n) {
    abs(mean(rows$approximate_target_sd[usable]) /
          mean(rows$exact_target_sd[usable]) - 1)
  } else NA_real_
  sd_status <- if (n < 5L || !is.finite(relative_sd_error)) {
    "inconclusive"
  } else if (relative_sd_error <= 0.05) "pass" else "fail"
  status <- if (sd_status == "fail" || interval_status == "fail") "fail" else
    if (sd_status == "pass" && interval_status == "pass") "pass" else
      "inconclusive"
  data.frame(
    case = rows$case[[1L]], grid_size = rows$grid_size[[1L]],
    ordering = rows$ordering[[1L]], neighbors = rows$neighbors[[1L]],
    target_id = rows$target_id[[1L]], attempted_replicates = nrow(rows),
    usable_replicates = n, failed_replicates = sum(!usable),
    relative_sd_error = relative_sd_error,
    sd_status = sd_status,
    threshold_probability_difference = mean_difference,
    threshold_difference_mcse = mcse,
    threshold_difference_ci_lower = lower,
    threshold_difference_ci_upper = upper,
    threshold_status = interval_status, status = status,
    median_relative_frobenius_covariance_error =
      if (any(usable)) median(rows$relative_frobenius_covariance_error[usable])
      else NA_real_,
    median_maximum_marginal_variance_error =
      if (any(usable)) median(rows$maximum_absolute_marginal_variance_error[usable])
      else NA_real_,
    median_target_elapsed_seconds =
      if (any(usable)) median(rows$target_elapsed_seconds[usable]) else NA_real_,
    median_peak_active_point_vectors =
      if (any(usable)) median(rows$peak_active_point_vectors[usable]) else
        NA_real_,
    stringsAsFactors = FALSE
  )
})
target_summary <- do.call(rbind, summaries)
utils::write.csv(target_summary, target_output, row.names = FALSE, na = "")

configuration_groups <- split(target_summary, do.call(interaction,
  c(unname(as.list(target_summary[
    c("case", "grid_size", "ordering", "neighbors")])),
    list(drop = TRUE, lex.order = TRUE))))
configuration <- do.call(rbind, lapply(configuration_groups, function(rows) {
  complete_targets <- setequal(rows$target_id,
    c("overall", "west_south", "west_north", "east_south", "east_north"))
  all_pass <- complete_targets && nrow(rows) == 5L && all(rows$status == "pass")
  data.frame(case = rows$case[[1L]], grid_size = rows$grid_size[[1L]],
    ordering = rows$ordering[[1L]], neighbors = rows$neighbors[[1L]],
    all_targets_pass = all_pass,
    target_status = paste(rows$target_id, rows$status, sep = ":", collapse = ";"),
    stringsAsFactors = FALSE)
}))
recommendation_groups <- split(configuration, do.call(interaction,
  c(unname(as.list(configuration[c("case", "grid_size")])),
    list(drop = TRUE, lex.order = TRUE))))
recommendations <- do.call(rbind, lapply(recommendation_groups, function(rows) {
  common <- sort(unique(rows$neighbors))
  selected <- NA_integer_
  for (budget in common) {
    both <- rows[rows$neighbors == budget, , drop = FALSE]
    if (all(c("input", "maximin") %in% both$ordering) &&
        all(both$all_targets_pass)) {
      selected <- as.integer(budget)
      break
    }
  }
  data.frame(case = rows$case[[1L]], grid_size = rows$grid_size[[1L]],
    recommended_neighbors = selected,
    status = if (is.na(selected)) "no budget passed both orders" else "pass",
    stringsAsFactors = FALSE)
}))
utils::write.csv(recommendations, recommendation_output,
                 row.names = FALSE, na = "")

order_keys <- c("case", "grid_size", "replicate", "response_seed",
                "neighbors", "target_id")
order_metrics <- c("target_mean_difference", "target_sd_difference",
  "threshold_probability_difference", "relative_frobenius_covariance_error",
  "maximum_absolute_marginal_variance_error", "target_elapsed_seconds",
  "peak_active_point_vectors")
input_order <- data[data$ordering == "input", c(order_keys, order_metrics),
                    drop = FALSE]
maximin_order <- data[data$ordering == "maximin",
                      c(order_keys, order_metrics), drop = FALSE]
paired <- merge(input_order, maximin_order, by = order_keys,
                suffixes = c("_input", "_maximin"), all = TRUE)
for (metric in order_metrics) {
  paired[[paste0(metric, "_maximin_minus_input")]] <-
    paired[[paste0(metric, "_maximin")]] -
    paired[[paste0(metric, "_input")]]
}
utils::write.csv(paired, order_output, row.names = FALSE, na = "")
cat("Target summaries: ", nrow(target_summary), " -> ", target_output,
    "\nRecommendations: ", nrow(recommendations), " -> ",
    recommendation_output, "\nOrder sensitivity: ", nrow(paired), " -> ",
    order_output, "\n", sep = "")
print(recommendations, row.names = FALSE)
