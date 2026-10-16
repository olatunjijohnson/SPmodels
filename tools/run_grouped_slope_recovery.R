args <- commandArgs(trailingOnly = TRUE)
archive_sha256 <- if (length(args) >= 1L) args[[1L]] else NA_character_
output_path <- if (length(args) >= 2L) args[[2L]] else
  "docs/validation/data/grouped_slope_recovery.csv"
start_scenario <- if (length(args) >= 3L) as.integer(args[[3L]]) else 1L

for (file in list.files("R", full.names = TRUE)) source(file)

scenarios <- list(
  list(id = "GRS-01", groups = 8L,
       covariance = diag(c(0.50, 0.20))),
  list(id = "GRS-02", groups = 8L,
       covariance = matrix(c(0.70, 0.18, 0.18, 0.25), 2L)),
  list(id = "GRS-03", groups = 20L,
       covariance = matrix(c(0.70, 0.18, 0.18, 0.25), 2L)),
  list(id = "GRS-04", groups = 20L,
       covariance = matrix(c(0.30, -0.12, -0.12, 0.25), 2L))
)
targets <- c(100L, 200L, 400L)
mean_truth <- c("(Intercept)" = 0.3, z = 0.4)
parameter_names <- c(
  "group:z::var((Intercept))", "group:z::cov(z,(Intercept))",
  "group:z::var(z)", "spatial_variance", "scale", "nugget_variance",
  "(Intercept)", "z"
)
thresholds <- c(0.03, 0.02, 0.03, 0.02, 0.15, 0.02, 0.03, 0.03)
names(thresholds) <- parameter_names
truth_for <- function(covariance) c(
  "group:z::var((Intercept))" = covariance[1L, 1L],
  "group:z::cov(z,(Intercept))" = covariance[2L, 1L],
  "group:z::var(z)" = covariance[2L, 2L],
  "spatial_variance" = 0.30, scale = 3, nugget_variance = 0.15,
  mean_truth
)

make_data <- function(scenario, seed) {
  set.seed(seed)
  groups <- scenario$groups
  east_count <- ceiling(sqrt(groups))
  north_count <- ceiling(groups / east_count)
  centres <- expand.grid(
    east = seq(0, by = 2, length.out = east_count),
    north = seq(0, by = 2, length.out = north_count)
  )[seq_len(groups), , drop = FALSE]
  group <- rep(seq_len(groups), each = 3L)
  z <- rep(c(-1, 0, 1), times = groups)
  coordinates <- as.matrix(centres[group, , drop = FALSE]) +
    matrix(runif(2L * length(group), -0.1, 0.1), ncol = 2L)
  colnames(coordinates) <- c("east", "north")
  group_effects <- matrix(rnorm(2L * groups), groups, 2L) %*%
    chol(scenario$covariance)
  distance <- .sp_dist(coordinates)
  covariance <- 0.30 * exp(-distance / 3)
  diag(covariance) <- diag(covariance) + 0.15
  response <- 0.3 + 0.4 * z +
    rowSums(cbind(1, z) * group_effects[group, , drop = FALSE]) +
    as.vector(t(chol(covariance)) %*% rnorm(length(group)))
  data.frame(y = response, z = z, group = factor(group),
             east = coordinates[, 1L], north = coordinates[, 2L])
}

fit_one <- function(scenario, scenario_index, replicate) {
  response_seed <- 410001L + 10000L * scenario_index + replicate
  data <- make_data(scenario, response_seed)
  warnings <- character()
  error_message <- NA_character_
  fit <- NULL
  elapsed <- system.time({
    fit <- tryCatch(withCallingHandlers(
      sp_fit(
        y ~ z, data, coords = ~ east + north,
        random_slope = sp_random_slope(
          "group", "z", diag(c(0.4, 0.2)), estimate = TRUE
        ),
        control = sp_control(
          smoothness = 0.5, spatial_variance = 0.4, scale = 2.5,
          nugget_variance = 0.2, estimate = TRUE, distance_units = "unit"
        )
      ),
      warning = function(warning) {
        warnings <<- c(warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }
    ), error = function(error) {
      error_message <<- conditionMessage(error)
      NULL
    })
  })[["elapsed"]]

  estimates <- stats::setNames(rep(NA_real_, length(parameter_names)),
                               parameter_names)
  standard_errors <- estimates
  convergence <- NA_integer_
  optimizer_message <- NA_character_
  covariance_pd <- FALSE
  covariance_finite <- FALSE
  cholesky_reconstruction_error <- NA_real_
  jacobian_max_abs_error <- NA_real_
  if (!is.null(fit)) {
    estimates[names(.sp_random_slope_parameters(
      fit$random_slope_specification
    ))] <- .sp_random_slope_parameters(fit$random_slope_specification)
    estimates[names(fit$covariance)] <- fit$covariance
    estimates[names(fit$coefficients)] <- fit$coefficients
    standard_errors[names(estimates)[1:3]] <- if (
      all(is.finite(diag(fit$random_slope_vcov)) &
            diag(fit$random_slope_vcov) >= 0)
    ) sqrt(diag(fit$random_slope_vcov)) else NA_real_
    standard_errors[c("spatial_variance", "scale", "nugget_variance")] <-
      if (all(is.finite(diag(fit$covariance_vcov)) &
              diag(fit$covariance_vcov) >= 0)) {
        sqrt(diag(fit$covariance_vcov))
      } else NA_real_
    standard_errors[c("(Intercept)", "z")] <- if (
      all(is.finite(diag(fit$coefficient_vcov)) &
            diag(fit$coefficient_vcov) >= 0)
    ) sqrt(diag(fit$coefficient_vcov)) else NA_real_
    convergence <- fit$optimisation$convergence
    optimizer_message <- if (is.null(fit$optimisation$message)) "" else
      fit$optimisation$message
    group_covariance <- fit$random_slope[[1L]]$covariance
    covariance_pd <- min(eigen(group_covariance, symmetric = TRUE,
                               only.values = TRUE)$values) > 0
    covariance_finite <- all(is.finite(fit$random_slope_vcov)) &&
      all(is.finite(fit$covariance_vcov)) &&
      all(is.finite(fit$coefficient_vcov))
    group_chol <- chol(group_covariance)
    cholesky_reconstruction_error <- max(abs(
      t(group_chol) %*% group_chol - group_covariance
    ))
    parameterization <- .sp_random_slope_parameterization(
      fit$random_slope_specification
    )
    analytic_jacobian <- .sp_random_slope_natural_jacobian(
      fit$random_slope_specification
    )
    numeric_jacobian <- analytic_jacobian
    for (column in seq_along(parameterization$values)) {
      step <- 1e-6 * max(1, abs(parameterization$values[[column]]))
      plus <- parameterization$values
      minus <- parameterization$values
      plus[[column]] <- plus[[column]] + step
      minus[[column]] <- minus[[column]] - step
      natural_plus <- .sp_random_slope_parameters(.sp_random_slope_update(
        fit$random_slope_specification, plus, parameterization
      ))
      natural_minus <- .sp_random_slope_parameters(.sp_random_slope_update(
        fit$random_slope_specification, minus, parameterization
      ))
      numeric_jacobian[, column] <- (natural_plus - natural_minus) /
        (2 * step)
    }
    jacobian_max_abs_error <- max(abs(analytic_jacobian - numeric_jacobian))
  }
  mathematical_checks_pass <- is.finite(cholesky_reconstruction_error) &&
    cholesky_reconstruction_error <= 1e-12 &&
    is.finite(jacobian_max_abs_error) && jacobian_max_abs_error <= 1e-6
  usable <- !is.null(fit) && !length(warnings) &&
    identical(convergence, 0L) && covariance_pd &&
    mathematical_checks_pass && all(is.finite(estimates))
  answer <- data.frame(
    archive_sha256 = archive_sha256,
    r_version = as.character(getRversion()), platform = R.version$platform,
    operating_system = R.version$os,
    blas = unname(extSoftVersion()["BLAS"]),
    scenario = scenario$id, scenario_index = scenario_index,
    groups = scenario$groups, replicate = replicate,
    response_seed = response_seed, fit_elapsed_seconds = unname(elapsed),
    usable = usable, convergence_code = convergence,
    optimizer_message = optimizer_message,
    warnings = paste(warnings, collapse = " | "),
    error_message = error_message, group_covariance_pd = covariance_pd,
    finite_covariance_uncertainty = covariance_finite,
    cholesky_reconstruction_error = cholesky_reconstruction_error,
    jacobian_max_abs_error = jacobian_max_abs_error,
    mathematical_checks_pass = mathematical_checks_pass,
    stringsAsFactors = FALSE
  )
  for (name in parameter_names) answer[[paste0("estimate__", name)]] <-
    unname(estimates[[name]])
  for (name in parameter_names) answer[[paste0("se__", name)]] <-
    unname(standard_errors[[name]])
  true_values <- truth_for(scenario$covariance)
  for (name in parameter_names) {
    answer[[paste0("covered__", name)]] <- if (
      is.finite(estimates[[name]]) && is.finite(standard_errors[[name]])
    ) {
      abs(estimates[[name]] - true_values[[name]]) <=
        1.96 * standard_errors[[name]]
    } else NA
  }
  answer
}

if (is.na(start_scenario) || start_scenario < 1L ||
    start_scenario > length(scenarios)) {
  stop("start_scenario must be between 1 and 4.", call. = FALSE)
}
if (start_scenario > 1L) {
  if (!file.exists(output_path)) {
    stop("Resume requires an existing raw-results CSV.", call. = FALSE)
  }
  existing <- utils::read.csv(output_path, stringsAsFactors = FALSE,
                              check.names = FALSE)
  all_results <- lapply(seq_len(nrow(existing)), function(i) {
    existing[i, , drop = FALSE]
  })
} else {
  if (file.exists(output_path)) {
    stop("Refusing to overwrite existing results; choose a new output path or resume.",
         call. = FALSE)
  }
  all_results <- list()
}
result_index <- length(all_results)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
for (scenario_index in seq.int(start_scenario, length(scenarios))) {
  scenario <- scenarios[[scenario_index]]
  scenario_results <- list()
  success_count <- 0L
  attempted <- 0L
  final_status <- "inconclusive"
  for (target_successes in targets) {
    while (success_count < target_successes && attempted < 400L) {
      attempted <- attempted + 1L
      row <- fit_one(scenario, scenario_index, attempted)
      result_index <- result_index + 1L
      scenario_results[[length(scenario_results) + 1L]] <- row
      all_results[[result_index]] <- row
      success_count <- success_count + as.integer(row$usable)
      current <- do.call(rbind, all_results)
      utils::write.csv(current, output_path, row.names = FALSE, na = "")
      cat(scenario$id, "attempt", attempted, "usable", row$usable,
          "successful_total", success_count, "target", target_successes,
          "\n")
    }
    rows <- do.call(rbind, scenario_results)
    usable_rows <- rows[rows$usable, , drop = FALSE]
    usable_fraction <- nrow(usable_rows) / nrow(rows)
    if (usable_fraction < 0.90 || !nrow(usable_rows)) {
      final_status <- "inconclusive: fewer than 90% usable"
      break
    }
    truth <- truth_for(scenario$covariance)
    estimates <- vapply(names(truth), function(name) {
      usable_rows[[paste0("estimate__", name)]]
    }, numeric(nrow(usable_rows)))
    bias_mcse <- apply(estimates, 2L, stats::sd) / sqrt(nrow(estimates))
    if (all(is.finite(bias_mcse)) &&
        all(bias_mcse <= thresholds[names(truth)])) {
      final_status <- "bias-MCSE targets passed"
      break
    }
    if (target_successes == 400L || attempted >= 400L) {
      final_status <- "inconclusive: bias-MCSE targets not met by 400 attempts"
      break
    }
  }
  rows <- do.call(rbind, scenario_results)
  truth <- truth_for(scenario$covariance)
  usable_rows <- rows[rows$usable, , drop = FALSE]
  estimates <- vapply(names(truth), function(name) {
    usable_rows[[paste0("estimate__", name)]]
  }, numeric(nrow(usable_rows)))
  bias_mcse <- if (nrow(estimates) > 1L) apply(estimates, 2L, stats::sd) /
    sqrt(nrow(estimates)) else stats::setNames(rep(NA_real_, length(truth)),
                                                names(truth))
  cat("completed", scenario$id, "status:", final_status,
      "usable", nrow(usable_rows), "of", nrow(rows), "\n")
}
cat("wrote", nrow(do.call(rbind, all_results)), "attempts to", output_path,
    "\n")
