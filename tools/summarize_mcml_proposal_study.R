args <- commandArgs(trailingOnly = TRUE)
fit_path <- if (length(args) >= 1L) args[[1L]] else
  "docs/validation/data/mcml_proposal_study.csv"
summary_path <- if (length(args) >= 2L) args[[2L]] else
  "docs/validation/data/mcml_proposal_summary.csv"
operational_path <- if (length(args) >= 3L) args[[3L]] else
  "docs/validation/data/mcml_proposal_operational.csv"

fits <- utils::read.csv(fit_path, stringsAsFactors = FALSE,
                        check.names = FALSE)
fits$attempt_gate <- fits$effective_sample_fraction >= 0.10 &
  fits$relative_mcse <= 0.05
groups <- split(fits, interaction(fits$case, fits$policy, fits$nsim,
                                  drop = TRUE))
summary_rows <- lapply(groups, function(rows) {
  data.frame(
    case = rows$case[[1L]], family = rows$family[[1L]], n = rows$n[[1L]],
    policy = rows$policy[[1L]], proposal_scales = rows$proposal_scales[[1L]],
    adapt_requested = rows$adapt_requested[[1L]], nsim = rows$nsim[[1L]],
    attempts = nrow(rows), gate_passes = sum(rows$attempt_gate, na.rm = TRUE),
    optimizer_converged = sum(rows$optimizer_convergence == 0L,
                              na.rm = TRUE),
    warning_attempts = sum(nzchar(rows$warnings)),
    error_attempts = sum(!is.na(rows$error_message) &
                           nzchar(rows$error_message)),
    adaptation_used = sum(rows$adapt_used %in% TRUE, na.rm = TRUE),
    median_pilot_ess_fraction = if (all(is.na(rows$pilot_ess_fraction)))
      NA_real_ else median(rows$pilot_ess_fraction, na.rm = TRUE),
    median_ess_fraction = median(rows$effective_sample_fraction,
                                 na.rm = TRUE),
    median_maximum_weight = median(rows$maximum_weight, na.rm = TRUE),
    median_relative_mcse = median(rows$relative_mcse, na.rm = TRUE),
    shared_objective_mean = mean(rows$shared_parameter_log_likelihood,
                                 na.rm = TRUE),
    shared_objective_sd = stats::sd(rows$shared_parameter_log_likelihood,
                                    na.rm = TRUE),
    optimized_objective_mean = mean(rows$optimized_log_likelihood,
                                    na.rm = TRUE),
    median_elapsed_seconds = median(rows$elapsed_seconds, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})
summary <- do.call(rbind, summary_rows)
utils::write.csv(summary, summary_path, row.names = FALSE, na = "")

at_2000 <- fits[fits$nsim == 2000L, , drop = FALSE]
operational <- do.call(rbind, lapply(split(at_2000,
  interaction(at_2000$case, at_2000$policy, drop = TRUE)), function(rows) {
    passed <- sum(rows$attempt_gate, na.rm = TRUE)
    data.frame(
      case = rows$case[[1L]], policy = rows$policy[[1L]],
      attempts = nrow(rows), gate_passes = passed,
      required_passes = 4L, operational = nrow(rows) == 5L && passed >= 4L,
      stringsAsFactors = FALSE
    )
  }))
utils::write.csv(operational, operational_path, row.names = FALSE,
                 na = "")
print(summary[, c("case", "policy", "nsim", "gate_passes",
                  "optimizer_converged", "warning_attempts",
                  "median_ess_fraction", "median_relative_mcse")],
      row.names = FALSE)
print(operational, row.names = FALSE)
