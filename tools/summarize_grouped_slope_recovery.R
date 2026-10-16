args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1L) args[[1L]] else
  "docs/validation/data/grouped_slope_recovery.csv"
output_path <- if (length(args) >= 2L) args[[2L]] else
  "docs/validation/data/grouped_slope_recovery_summary.csv"

data <- utils::read.csv(input_path, stringsAsFactors = FALSE,
                        check.names = FALSE)
if (!"mathematical_checks_pass" %in% names(data)) {
  stop("Run validate_grouped_slope_recovery.R before summarizing.",
       call. = FALSE)
}
data$usable <- data$usable & data$mathematical_checks_pass
truths <- list(
  "GRS-01" = c("group:z::var((Intercept))" = 0.50,
                "group:z::cov(z,(Intercept))" = 0,
                "group:z::var(z)" = 0.20,
                spatial_variance = 0.30, scale = 3,
                nugget_variance = 0.15, "(Intercept)" = 0.3, z = 0.4),
  "GRS-02" = c("group:z::var((Intercept))" = 0.70,
                "group:z::cov(z,(Intercept))" = 0.18,
                "group:z::var(z)" = 0.25,
                spatial_variance = 0.30, scale = 3,
                nugget_variance = 0.15, "(Intercept)" = 0.3, z = 0.4),
  "GRS-03" = c("group:z::var((Intercept))" = 0.70,
                "group:z::cov(z,(Intercept))" = 0.18,
                "group:z::var(z)" = 0.25,
                spatial_variance = 0.30, scale = 3,
                nugget_variance = 0.15, "(Intercept)" = 0.3, z = 0.4),
  "GRS-04" = c("group:z::var((Intercept))" = 0.30,
                "group:z::cov(z,(Intercept))" = -0.12,
                "group:z::var(z)" = 0.25,
                spatial_variance = 0.30, scale = 3,
                nugget_variance = 0.15, "(Intercept)" = 0.3, z = 0.4)
)
thresholds <- c(0.03, 0.02, 0.03, 0.02, 0.15, 0.02, 0.03, 0.03)
names(thresholds) <- names(truths[[1L]])

groups <- split(data, data$scenario)
summary_rows <- lapply(names(truths), function(scenario) {
  rows <- groups[[scenario]]
  truth <- truths[[scenario]]
  usable <- rows[rows$usable, , drop = FALSE]
  estimates <- if (nrow(usable)) {
    vapply(names(truth), function(name) {
      usable[[paste0("estimate__", name)]]
    }, numeric(nrow(usable)))
  } else matrix(numeric(), nrow = 0L, ncol = length(truth),
                dimnames = list(NULL, names(truth)))
  means <- if (nrow(estimates)) colMeans(estimates) else
    stats::setNames(rep(NA_real_, length(truth)), names(truth))
  bias <- means - truth
  bias_mcse <- if (nrow(estimates) > 1L) apply(estimates, 2L, stats::sd) /
    sqrt(nrow(estimates)) else stats::setNames(rep(NA_real_, length(truth)),
                                                names(truth))
  usable_fraction <- nrow(usable) / nrow(rows)
  precision_pass <- all(is.finite(bias_mcse)) &&
    all(bias_mcse <= thresholds[names(truth)])
  status <- if (usable_fraction < 0.90) {
    "inconclusive: fewer than 90% usable"
  } else if (precision_pass) {
    "bias-MCSE targets passed (precision only)"
  } else {
    "inconclusive: bias-MCSE targets not met"
  }
  answer <- data.frame(
    scenario = scenario, groups = rows$groups[[1L]],
    attempted = nrow(rows), usable = nrow(usable),
    usable_fraction = usable_fraction,
    mathematical_checks_passed = sum(rows$mathematical_checks_pass),
    finite_hessian_covariance = sum(rows$finite_covariance_uncertainty),
    hessian_covariance_fraction = mean(rows$finite_covariance_uncertainty),
    median_fit_elapsed_seconds = median(rows$fit_elapsed_seconds),
    status = status, stringsAsFactors = FALSE
  )
  for (name in names(truth)) {
    estimate <- rows[[paste0("estimate__", name)]]
    standard_error <- rows[[paste0("se__", name)]]
    coverage <- ifelse(
      is.finite(estimate) & is.finite(standard_error),
      abs(estimate - truth[[name]]) <= 1.96 * standard_error, NA
    )
    answer[[paste0("truth__", name)]] <- truth[[name]]
    answer[[paste0("mean__", name)]] <- means[[name]]
    answer[[paste0("bias__", name)]] <- bias[[name]]
    answer[[paste0("bias_mcse__", name)]] <- bias_mcse[[name]]
    answer[[paste0("bias_mcse_target__", name)]] <- thresholds[[name]]
    answer[[paste0("wald_n__", name)]] <- sum(!is.na(coverage))
    answer[[paste0("local_wald_coverage__", name)]] <-
      if (any(!is.na(coverage))) mean(coverage, na.rm = TRUE) else NA_real_
  }
  answer
})
summary <- do.call(rbind, summary_rows)
utils::write.csv(summary, output_path, row.names = FALSE, na = "")
print(summary[, c("scenario", "attempted", "usable", "usable_fraction",
                  "median_fit_elapsed_seconds", "status")], row.names = FALSE)
