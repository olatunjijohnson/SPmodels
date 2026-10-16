args <- commandArgs(trailingOnly = TRUE)
fit_path <- if (length(args) >= 1L) args[[1L]] else
  "docs/validation/data/temporal_scale_recovery.csv"
mcml_path <- if (length(args) >= 2L) args[[2L]] else
  "docs/validation/data/temporal_scale_mcml_repeats.csv"

fits <- utils::read.csv(fit_path, stringsAsFactors = FALSE,
                        check.names = FALSE)
mcml <- utils::read.csv(mcml_path, stringsAsFactors = FALSE,
                        check.names = FALSE)
scenario_ids <- c("TSC-01", "TSC-02", "TSC-03", "TSC-04")
stopifnot(nrow(fits) == 4L * 400L)
stopifnot(identical(unique(fits$scenario), scenario_ids))
stopifnot(length(unique(fits$archive_sha256)) == 1L)
stopifnot(all(is.finite(fits$fit_elapsed_seconds) &
                fits$fit_elapsed_seconds > 0))

for (i in seq_along(scenario_ids)) {
  rows <- fits[fits$scenario == scenario_ids[[i]], , drop = FALSE]
  stopifnot(nrow(rows) == 400L)
  stopifnot(identical(rows$replicate, seq_len(400L)))
  stopifnot(identical(rows$response_seed,
                      410001L + 10000L * i + seq_len(400L)))
  finite_estimates <- apply(rows[grep("^estimate__", names(rows))], 1L,
                            function(x) all(is.finite(x)))
  expected_usable <- finite_estimates &
    rows$optimizer_convergence == 0L & rows$latent_mode_converged &
    rows$covariance_positive_definite & !nzchar(rows$fit_warnings) &
    is.na(rows$fit_error)
  stopifnot(all(rows$usable == expected_usable))

  repeats <- mcml[mcml$scenario == scenario_ids[[i]], , drop = FALSE]
  successful_replicates <- rows$replicate[rows$usable]
  expected_replicates <- head(successful_replicates, 25L)
  response_replicates <- sort(unique(repeats$replicate))
  stopifnot(all(response_replicates == expected_replicates))
  for (replicate in expected_replicates) {
    one <- repeats[repeats$replicate == replicate, , drop = FALSE]
    stopifnot(nrow(one) == 3L)
    stopifnot(all(one$repeat_index == 1:3))
    response_seed <- 410001L + 10000L * i + replicate
    stopifnot(all(one$mcml_seed == response_seed + 5000L + 0:2))
  }
}

expected_repeat_gate <- mcml$optimizer_convergence == 0L &
  !nzchar(mcml$warnings) &
  is.finite(mcml$effective_sample_fraction) &
  mcml$effective_sample_fraction >= 0.10 &
  is.finite(mcml$relative_mcse) & mcml$relative_mcse <= 0.05 &
  (is.na(mcml$error_message) | !nzchar(mcml$error_message))
stopifnot(all(mcml$repeat_gate_pass == expected_repeat_gate))
stopifnot(nrow(mcml) == 4L * 25L * 3L)

cat("TSC raw-row, seed, usability and MCML repeat/gate checks passed.\n")
