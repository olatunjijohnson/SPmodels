args <- commandArgs(trailingOnly = TRUE)
fit_path <- if (length(args) >= 1L) args[[1L]] else
  "docs/validation/data/mcml_proposal_study.csv"
fixture_path <- if (length(args) >= 2L) args[[2L]] else
  "docs/validation/data/mcml_proposal_fixtures.csv"
base_path <- if (length(args) >= 3L) args[[3L]] else
  "docs/validation/data/mcml_proposal_study_base.csv"
loa_path <- if (length(args) >= 4L) args[[4L]] else
  "docs/validation/data/mcml_proposal_loa.csv"

fits <- utils::read.csv(fit_path, stringsAsFactors = FALSE,
                        check.names = FALSE)
fixtures <- utils::read.csv(fixture_path, stringsAsFactors = FALSE,
                           check.names = FALSE)
base <- utils::read.csv(base_path, stringsAsFactors = FALSE,
                        check.names = FALSE)
stopifnot(nrow(fits) == 2L * 4L * 3L * 5L)
stopifnot(nrow(fixtures) == 24L + 48L, nrow(base) == 2L)
stopifnot(length(unique(fits$archive_sha256)) == 1L)
stopifnot(all(base$base_convergence == 0L),
          all(is.na(base$base_warnings) | !nzchar(base$base_warnings)),
          all(is.na(base$base_error) | !nzchar(base$base_error)))
stopifnot(all(is.finite(fits$shared_parameter_log_likelihood)))
stopifnot(all(fits$usable_base_fit))

policy_names <- c("single_fixed", "mixture_fixed",
                  "single_adaptive", "mixture_adaptive")
for (case_index in 1:2) for (policy_index in 1:4) {
  policy <- policy_names[[policy_index]]
  rows <- fits[fits$case_index == case_index &
    fits$policy_index == policy_index, , drop = FALSE]
  stopifnot(nrow(rows) == 15L)
  for (size_index in 1:3) {
    one <- rows[rows$size_index == size_index, , drop = FALSE]
    stopifnot(nrow(one) == 5L)
    stopifnot(all(one$attempt == 1:5))
    expected_seed <- 520001L + 100000L * case_index +
      10000L * policy_index + 100L * size_index + 1:5
    stopifnot(all(one$seed == expected_seed))
  }
  expected_adapt <- policy %in% c("single_adaptive", "mixture_adaptive")
  stopifnot(all(rows$adapt_requested == expected_adapt))
}

expected_gate <- fits$effective_sample_fraction >= 0.10 &
  fits$relative_mcse <= 0.05
stopifnot(length(expected_gate) == nrow(fits), !anyNA(expected_gate))
cat("MCP fixture, design, seed and attempt-gate checks passed.\n")

if (file.exists(loa_path)) {
  loa <- utils::read.csv(loa_path, stringsAsFactors = FALSE,
                         check.names = FALSE)
  stopifnot(nrow(loa) == 20L,
            all(loa$lsda_version == "0.2.0"),
            all(loa$input_sha256 ==
                  "0acc9426c9842f93c264914a86e62641c160bdcdf134391cc57be7facb3c4811"),
            all(loa$optimizer_convergence == 0L),
            all(nzchar(loa$warnings)),
            all(is.na(loa$error_message) | !nzchar(loa$error_message)))
  for (policy_index in 1:4) {
    rows <- loa[loa$policy_index == policy_index, , drop = FALSE]
    stopifnot(nrow(rows) == 5L,
              all(rows$attempt == 1:5),
              all(rows$seed == 730001L + 10000L * policy_index + 1:5))
  }
  cat("Loa input identity, policy/seed and diagnostics checks passed.\n")
}
