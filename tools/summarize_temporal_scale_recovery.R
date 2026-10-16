args <- commandArgs(trailingOnly = TRUE)
fit_path <- if (length(args) >= 1L) args[[1L]] else
  "docs/validation/data/temporal_scale_recovery.csv"
mcml_path <- if (length(args) >= 2L) args[[2L]] else
  "docs/validation/data/temporal_scale_mcml_repeats.csv"
summary_path <- if (length(args) >= 3L) args[[3L]] else
  "docs/validation/data/temporal_scale_recovery_summary.csv"
mcml_summary_path <- if (length(args) >= 4L) args[[4L]] else
  "docs/validation/data/temporal_scale_mcml_summary.csv"

fits <- utils::read.csv(fit_path, stringsAsFactors = FALSE,
                        check.names = FALSE)
scenario_truth <- list(
  "TSC-01" = c(temporal_scale = 0.75,
               log_temporal_scale = log(0.75), spatial_variance = 0.45,
               spatial_scale = 2, nugget_variance = 0.20,
               `(Intercept)` = -0.5, z = 0.35),
  "TSC-02" = c(temporal_scale = 2.00,
               log_temporal_scale = log(2.00), spatial_variance = 0.45,
               spatial_scale = 2, nugget_variance = 0.20,
               `(Intercept)` = -0.5, z = 0.35),
  "TSC-03" = c(temporal_scale = 0.75,
               log_temporal_scale = log(0.75), spatial_variance = 0.30,
               spatial_scale = 2, nugget_variance = 0.15,
               `(Intercept)` = 0.2, z = 0.20),
  "TSC-04" = c(temporal_scale = 2.00,
               log_temporal_scale = log(2.00), spatial_variance = 0.30,
               spatial_scale = 2, nugget_variance = 0.15,
               `(Intercept)` = 0.2, z = 0.20)
)

fit_summaries <- lapply(names(scenario_truth), function(scenario) {
  rows <- fits[fits$scenario == scenario, , drop = FALSE]
  usable <- rows[rows$usable, , drop = FALSE]
  truth <- scenario_truth[[scenario]]
  target <- max(0.05, 0.05 * truth[["temporal_scale"]])
  estimates <- if (nrow(usable)) {
    vapply(names(truth), function(name) {
      usable[[paste0("estimate__", name)]]
    }, numeric(nrow(usable)))
  } else matrix(numeric(), nrow = 0L, ncol = length(truth),
                dimnames = list(NULL, names(truth)))
  mean_estimate <- if (nrow(estimates)) colMeans(estimates) else
    stats::setNames(rep(NA_real_, length(truth)), names(truth))
  bias_mcse <- if (nrow(estimates) > 1L) apply(estimates, 2L, stats::sd) /
    sqrt(nrow(estimates)) else stats::setNames(rep(NA_real_, length(truth)),
                                                names(truth))
  temporal_mcse <- bias_mcse[["temporal_scale"]]
  usable_fraction <- nrow(usable) / nrow(rows)
  status <- if (usable_fraction < 0.90) {
    "inconclusive: fewer than 90% usable"
  } else if (is.finite(temporal_mcse) && temporal_mcse <= target) {
    "temporal bias-MCSE target met (precision only)"
  } else "inconclusive: temporal bias-MCSE target not met"
  answer <- data.frame(
    scenario = scenario, attempted = nrow(rows), usable = nrow(usable),
    usable_fraction = usable_fraction,
    median_fit_elapsed_seconds = median(rows$fit_elapsed_seconds),
    finite_covariance_uncertainty = sum(rows$covariance_uncertainty_finite),
    status = status, stringsAsFactors = FALSE
  )
  for (name in names(truth)) {
    coverage <- rows[[paste0("covered__", name)]]
    answer[[paste0("truth__", name)]] <- truth[[name]]
    answer[[paste0("mean__", name)]] <- mean_estimate[[name]]
    answer[[paste0("bias__", name)]] <-
      mean_estimate[[name]] - truth[[name]]
    answer[[paste0("bias_mcse__", name)]] <- bias_mcse[[name]]
    answer[[paste0("wald_n__", name)]] <- sum(!is.na(coverage))
    answer[[paste0("local_wald_coverage__", name)]] <-
      if (any(!is.na(coverage))) mean(coverage, na.rm = TRUE) else NA_real_
  }
  answer$temporal_bias_mcse_target <- target
  answer$log_temporal_bias_mcse <- bias_mcse[["log_temporal_scale"]]
  answer
})
fit_summary <- do.call(rbind, fit_summaries)
utils::write.csv(fit_summary, summary_path, row.names = FALSE, na = "")
print(fit_summary[, c("scenario", "attempted", "usable",
                      "temporal_bias_mcse_target",
                      "bias_mcse__temporal_scale", "status")],
      row.names = FALSE)

if (file.exists(mcml_path)) {
  mcml <- utils::read.csv(mcml_path, stringsAsFactors = FALSE,
                          check.names = FALSE)
  response_groups <- split(mcml, interaction(mcml$scenario,
    mcml$response_seed, drop = TRUE))
  response_summary <- do.call(rbind, lapply(response_groups, function(rows) {
    data.frame(
      scenario = rows$scenario[[1L]], response_seed = rows$response_seed[[1L]],
      repeats = nrow(rows), all_three_repeat_gate = nrow(rows) == 3L &&
        all(rows$repeat_gate_pass),
      median_ess_fraction = median(rows$effective_sample_fraction,
                                   na.rm = TRUE),
      maximum_relative_mcse = max(rows$relative_mcse, na.rm = TRUE),
      mean_log_temporal_movement = mean(
        rows$log_temporal_scale_movement, na.rm = TRUE
      ),
      between_repeat_temporal_covariance =
        rows$between_repeat_temporal_covariance[[1L]],
      warnings = paste(unique(stats::na.omit(rows$warnings)), collapse = " | "),
      errors = paste(unique(stats::na.omit(rows$error_message)), collapse = " | "),
      stringsAsFactors = FALSE
    )
  }))
  mcml_summary <- do.call(rbind, lapply(names(scenario_truth), function(scenario) {
    rows <- response_summary[response_summary$scenario == scenario, ,
                             drop = FALSE]
    data.frame(
      scenario = scenario, responses = nrow(rows),
      all_repeats_passed = sum(rows$all_three_repeat_gate),
      fraction_all_repeats_passed = if (nrow(rows))
        mean(rows$all_three_repeat_gate) else NA_real_,
      median_ess_fraction = if (nrow(rows))
        median(rows$median_ess_fraction, na.rm = TRUE) else NA_real_,
      median_maximum_relative_mcse = if (nrow(rows))
        median(rows$maximum_relative_mcse, na.rm = TRUE) else NA_real_,
      median_log_temporal_movement = if (nrow(rows))
        median(rows$mean_log_temporal_movement, na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  utils::write.csv(mcml_summary, mcml_summary_path,
                   row.names = FALSE, na = "")
  print(mcml_summary, row.names = FALSE)
}
