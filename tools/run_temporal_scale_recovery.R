args <- commandArgs(trailingOnly = TRUE)
archive_sha256 <- if (length(args) >= 1L) args[[1L]] else NA_character_
fit_path <- if (length(args) >= 2L) args[[2L]] else
  "docs/validation/data/temporal_scale_recovery.csv"
mcml_path <- if (length(args) >= 3L) args[[3L]] else
  "docs/validation/data/temporal_scale_mcml_repeats.csv"
for (file in list.files("R", full.names = TRUE)) source(file)

scenarios <- list(
  list(id = "TSC-01", family = "binomial", temporal_scale = 0.75,
       spatial_variance = 0.45, spatial_scale = 2, nugget_variance = 0.20,
       intercept = -0.5, slope = 0.35),
  list(id = "TSC-02", family = "binomial", temporal_scale = 2.00,
       spatial_variance = 0.45, spatial_scale = 2, nugget_variance = 0.20,
       intercept = -0.5, slope = 0.35),
  list(id = "TSC-03", family = "poisson", temporal_scale = 0.75,
       spatial_variance = 0.30, spatial_scale = 2, nugget_variance = 0.15,
       intercept = 0.2, slope = 0.20),
  list(id = "TSC-04", family = "poisson", temporal_scale = 2.00,
       spatial_variance = 0.30, spatial_scale = 2, nugget_variance = 0.15,
       intercept = 0.2, slope = 0.20)
)
targets <- c(100L, 200L, 400L)
design <- expand.grid(east = 0:3, north = 0:3,
                      year = c(0, 1, 2, 4))
design$z <- as.numeric(scale(design$east + 0.5 * design$north +
                              0.25 * design$year))

simulate_response <- function(scenario, seed) {
  set.seed(seed)
  distance <- .sp_dist(as.matrix(design[c("east", "north")]))
  time_lag <- abs(outer(design$year, design$year, "-"))
  covariance <- scenario$spatial_variance *
    exp(-distance / scenario$spatial_scale) *
    exp(-time_lag / scenario$temporal_scale)
  diag(covariance) <- diag(covariance) + scenario$nugget_variance
  latent <- as.vector(t(chol(covariance)) %*% rnorm(nrow(design)))
  eta <- scenario$intercept + scenario$slope * design$z + latent
  response <- if (scenario$family == "binomial") {
    stats::rbinom(nrow(design), size = 10L, prob = stats::plogis(eta))
  } else {
    stats::rpois(nrow(design), lambda = exp(eta))
  }
  data <- design
  data$y <- response
  data
}

fit_one <- function(scenario, scenario_index, replicate,
                    run_mcml = FALSE) {
  response_seed <- 410001L + 10000L * scenario_index + replicate
  data <- simulate_response(scenario, response_seed)
  formula <- if (scenario$family == "binomial") {
    stats::as.formula("cbind(y, 10 - y) ~ z")
  } else y ~ z
  family <- if (scenario$family == "binomial") stats::binomial() else
    stats::poisson()
  fit_warnings <- character()
  fit_error <- NA_character_
  fit <- NULL
  fit_elapsed <- system.time({
    fit <- tryCatch(withCallingHandlers(
      sp_fit(
        formula, data, family = family, coords = ~ east + north,
        time = sp_temporal(
          ~ year, scale = 1.25 * scenario$temporal_scale,
          time_units = "years", estimate = TRUE
        ),
        control = sp_control(
          smoothness = 0.5, spatial_variance = 0.4, scale = 1.5,
          nugget_variance = 0.2, estimate = TRUE,
          distance_units = "grid-unit"
        )
      ),
      warning = function(warning) {
        fit_warnings <<- c(fit_warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }
    ), error = function(error) {
      fit_error <<- conditionMessage(error)
      NULL
    })
  })[["elapsed"]]

  estimates <- c(temporal_scale = NA_real_, log_temporal_scale = NA_real_,
                 spatial_variance = NA_real_, spatial_scale = NA_real_,
                 nugget_variance = NA_real_,
                 `(Intercept)` = NA_real_, z = NA_real_)
  standard_errors <- stats::setNames(rep(NA_real_, length(estimates)),
                                     names(estimates))
  convergence <- NA_integer_
  latent_mode_converged <- FALSE
  covariance_positive_definite <- FALSE
  temporal_se <- NA_real_
  temporal_coverage <- NA
  covariance_uncertainty_finite <- FALSE
  fit_log_likelihood <- NA_real_
  if (!is.null(fit)) {
    estimates[c("temporal_scale", "spatial_variance", "spatial_scale",
                "nugget_variance", "(Intercept)", "z")] <- c(
      fit$temporal_scale,
      fit$covariance[c("spatial_variance", "scale", "nugget_variance")],
      fit$coefficients[c("(Intercept)", "z")]
    )
    estimates[["log_temporal_scale"]] <- log(fit$temporal_scale)
    convergence <- fit$optimisation$convergence
    latent_mode_converged <- isTRUE(fit$latent_mode$converged)
    covariance_positive_definite <- min(eigen(
      fit$fitted_state$sigma, symmetric = TRUE, only.values = TRUE
    )$values) > 0
    temporal_se <- fit$temporal_scale_se
    if (is.finite(temporal_se)) {
      temporal_coverage <- abs(fit$temporal_scale -
        scenario$temporal_scale) <= 1.96 * temporal_se
    }
    standard_errors[["temporal_scale"]] <- temporal_se
    if (all(is.finite(diag(fit$covariance_vcov))) &&
        all(diag(fit$covariance_vcov) >= 0)) {
      standard_errors[c("spatial_variance", "spatial_scale",
                        "nugget_variance")] <- sqrt(diag(
        fit$covariance_vcov[c("spatial_variance", "scale",
                              "nugget_variance"),
                            c("spatial_variance", "scale",
                              "nugget_variance"), drop = FALSE]
      ))
    }
    if (all(is.finite(diag(fit$coefficient_vcov))) &&
        all(diag(fit$coefficient_vcov) >= 0)) {
      standard_errors[c("(Intercept)", "z")] <- sqrt(diag(
        fit$coefficient_vcov[c("(Intercept)", "z"),
                             c("(Intercept)", "z"), drop = FALSE]
      ))
    }
    covariance_uncertainty_finite <-
      all(is.finite(fit$covariance_vcov)) &&
      all(is.finite(fit$parameter_vcov)) &&
      all(is.finite(fit$coefficient_vcov))
    fit_log_likelihood <- as.numeric(logLik(fit))
  }
  usable <- !is.null(fit) && !length(fit_warnings) &&
    identical(convergence, 0L) && latent_mode_converged &&
    covariance_positive_definite && all(is.finite(estimates))

  mcml_rows <- NULL
  mcml_error <- NA_character_
  mcml_warnings <- character()
  mcml_elapsed <- NA_real_
  if (!is.null(fit) && usable && run_mcml) {
    repeats <- NULL
    mcml_elapsed <- system.time({
      repeats <- tryCatch(withCallingHandlers(
        sp_mcml_replicates(
          fit, nrep = 3L,
          control = sp_mcml_control(nsim = 1000L),
          seed = response_seed + 5000L
        ),
        warning = function(warning) {
          mcml_warnings <<- c(mcml_warnings, conditionMessage(warning))
          invokeRestart("muffleWarning")
        }
      ), error = function(error) {
        mcml_error <<- conditionMessage(error)
        NULL
      })
    })[["elapsed"]]
    if (!is.null(repeats)) {
      temporal_name <- "temporal::scale"
      mcml_rows <- do.call(rbind, lapply(seq_along(repeats$fits), function(i) {
        mcml_fit <- repeats$fits[[i]]
        diagnostics <- repeats$diagnostics[i, , drop = FALSE]
        mcml_temporal <- mcml_fit$temporal_scale
        temporal_covariance <- repeats$natural_mc_covariance[
          temporal_name, temporal_name
        ]
        data.frame(
          archive_sha256 = archive_sha256, scenario = scenario$id,
          scenario_index = scenario_index, response_seed = response_seed,
          replicate = replicate, repeat_index = i,
          mcml_seed = diagnostics$seed,
          nsim = repeats$control$nsim,
          elapsed_seconds = mcml_elapsed,
          optimizer_convergence = mcml_fit$optimisation$convergence,
          optimizer_message = if (is.null(mcml_fit$optimisation$message)) "" else
            mcml_fit$optimisation$message,
          temporal_scale = mcml_temporal,
          log_temporal_scale_movement = log(mcml_temporal /
            fit$temporal_scale),
          log_likelihood = mcml_fit$log_likelihood,
          effective_sample_size = diagnostics$effective_sample_size,
          effective_sample_fraction = diagnostics$effective_sample_fraction,
          maximum_weight = diagnostics$maximum_weight,
          relative_mcse = diagnostics$relative_mcse,
          between_repeat_temporal_covariance = temporal_covariance,
          between_repeat_temporal_mcse = sqrt(temporal_covariance),
          repeat_gate_pass = isTRUE(mcml_fit$optimisation$convergence == 0L) &&
            !length(mcml_warnings) &&
            is.finite(diagnostics$effective_sample_fraction) &&
            diagnostics$effective_sample_fraction >= 0.10 &&
            is.finite(diagnostics$relative_mcse) &&
            diagnostics$relative_mcse <= 0.05,
          warnings = paste(mcml_warnings, collapse = " | "),
          error_message = mcml_error,
          stringsAsFactors = FALSE
        )
      }))
    } else {
      mcml_rows <- data.frame(
        archive_sha256 = archive_sha256, scenario = scenario$id,
        scenario_index = scenario_index, response_seed = response_seed,
        replicate = replicate, repeat_index = 1:3,
        mcml_seed = response_seed + 5000L + 0:2,
        nsim = 1000L, elapsed_seconds = mcml_elapsed,
        optimizer_convergence = NA_integer_, optimizer_message = NA_character_,
        temporal_scale = NA_real_, log_temporal_scale_movement = NA_real_,
        log_likelihood = NA_real_, effective_sample_size = NA_real_,
        effective_sample_fraction = NA_real_, maximum_weight = NA_real_,
        relative_mcse = NA_real_,
        between_repeat_temporal_covariance = NA_real_,
        between_repeat_temporal_mcse = NA_real_, repeat_gate_pass = FALSE,
        warnings = paste(mcml_warnings, collapse = " | "),
        error_message = mcml_error, stringsAsFactors = FALSE
      )
    }
  }
  row <- data.frame(
    archive_sha256 = archive_sha256,
    r_version = as.character(getRversion()), platform = R.version$platform,
    operating_system = R.version$os,
    blas = unname(extSoftVersion()["BLAS"]),
    scenario = scenario$id, scenario_index = scenario_index,
    family = scenario$family, replicate = replicate,
    response_seed = response_seed, fit_elapsed_seconds = unname(fit_elapsed),
    mcml_elapsed_seconds = unname(mcml_elapsed), usable = usable,
    optimizer_convergence = convergence,
    latent_mode_converged = latent_mode_converged,
    covariance_positive_definite = covariance_positive_definite,
    covariance_uncertainty_finite = covariance_uncertainty_finite,
    temporal_standard_error = temporal_se,
    temporal_local_wald_covered = temporal_coverage,
    fit_warnings = paste(fit_warnings, collapse = " | "),
    fit_error = fit_error, fit_log_likelihood = fit_log_likelihood,
    mcml_warnings = paste(mcml_warnings, collapse = " | "),
    mcml_error = mcml_error, stringsAsFactors = FALSE
  )
  for (name in names(estimates)) row[[paste0("estimate__", name)]] <-
    unname(estimates[[name]])
  for (name in names(standard_errors)) row[[paste0("se__", name)]] <-
    unname(standard_errors[[name]])
  truth <- c(
    temporal_scale = scenario$temporal_scale,
    log_temporal_scale = log(scenario$temporal_scale),
    spatial_variance = scenario$spatial_variance,
    spatial_scale = scenario$spatial_scale,
    nugget_variance = scenario$nugget_variance,
    `(Intercept)` = scenario$intercept, z = scenario$slope
  )
  for (name in names(estimates)) {
    row[[paste0("covered__", name)]] <- if (
      is.finite(estimates[[name]]) && is.finite(standard_errors[[name]])
    ) {
      abs(estimates[[name]] - truth[[name]]) <=
        1.96 * standard_errors[[name]]
    } else NA
  }
  list(fit = row, mcml = mcml_rows)
}

dir.create(dirname(fit_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(mcml_path), recursive = TRUE, showWarnings = FALSE)
fit_results <- list()
mcml_results <- list()
fit_count <- 0L
mcml_count <- 0L
if (file.exists(fit_path)) {
  existing_fits <- utils::read.csv(fit_path, stringsAsFactors = FALSE,
                                   check.names = FALSE)
  if (nrow(existing_fits)) {
    if (any(existing_fits$archive_sha256 != archive_sha256)) {
      stop("Existing fit rows use a different source archive SHA-256.",
           call. = FALSE)
    }
    fit_results <- lapply(seq_len(nrow(existing_fits)), function(i) {
      existing_fits[i, , drop = FALSE]
    })
    fit_count <- length(fit_results)
  }
}
if (file.exists(mcml_path)) {
  existing_mcml <- utils::read.csv(mcml_path, stringsAsFactors = FALSE,
                                   check.names = FALSE)
  if (nrow(existing_mcml)) {
    if (any(existing_mcml$archive_sha256 != archive_sha256)) {
      stop("Existing MCML rows use a different source archive SHA-256.",
           call. = FALSE)
    }
    mcml_results <- lapply(seq_len(nrow(existing_mcml)), function(i) {
      existing_mcml[i, , drop = FALSE]
    })
    mcml_count <- length(mcml_results)
  }
}
for (scenario_index in seq_along(scenarios)) {
  scenario <- scenarios[[scenario_index]]
  existing_scenario <- Filter(function(row) {
    identical(row$scenario[[1L]], scenario$id)
  }, fit_results)
  scenario_results <- existing_scenario
  usable_count <- sum(vapply(existing_scenario, function(row) {
    isTRUE(row$usable[[1L]])
  }, logical(1L)))
  attempted <- length(scenario_results)
  for (target_successes in targets) {
    while (usable_count < target_successes && attempted < 400L) {
      attempted <- attempted + 1L
      result <- fit_one(scenario, scenario_index, attempted,
                        run_mcml = usable_count < 25L)
      fit_count <- fit_count + 1L
      fit_results[[fit_count]] <- result$fit
      scenario_results[[length(scenario_results) + 1L]] <- result$fit
      usable_count <- usable_count + as.integer(result$fit$usable)
      if (!is.null(result$mcml)) {
        mcml_count <- mcml_count + 1L
        mcml_results[[mcml_count]] <- result$mcml
        utils::write.csv(do.call(rbind, mcml_results), mcml_path,
                         row.names = FALSE, na = "")
      }
      utils::write.csv(do.call(rbind, fit_results), fit_path,
                       row.names = FALSE, na = "")
      cat(scenario$id, "attempt", attempted, "usable",
          result$fit$usable, "successful_total", usable_count,
          "target", target_successes, "\n")
    }
    rows <- do.call(rbind, scenario_results)
    usable_rows <- rows[rows$usable, , drop = FALSE]
    scale_mcse <- if (nrow(usable_rows) > 1L) {
      stats::sd(usable_rows$`estimate__temporal_scale`) /
        sqrt(nrow(usable_rows))
    } else NA_real_
    target <- max(0.05, 0.05 * scenario$temporal_scale)
    if (is.finite(scale_mcse) && scale_mcse <= target) break
    if (target_successes == 400L || attempted >= 400L) break
  }
  rows <- do.call(rbind, scenario_results)
  usable_rows <- rows[rows$usable, , drop = FALSE]
  if (nrow(usable_rows) / nrow(rows) < 0.90) {
    status <- "inconclusive: fewer than 90% usable"
  } else {
    scale_mcse <- stats::sd(usable_rows$`estimate__temporal_scale`) /
      sqrt(nrow(usable_rows))
    target <- max(0.05, 0.05 * scenario$temporal_scale)
    status <- if (is.finite(scale_mcse) && scale_mcse <= target) {
      "temporal bias-MCSE target met (precision only)"
    } else "inconclusive: temporal bias-MCSE target not met"
  }
  cat("completed", scenario$id, "attempts", nrow(rows), "usable",
      nrow(usable_rows), "status:", status, "\n")
}
cat("wrote", fit_count, "Laplace fits and", mcml_count,
    "MCML replicate sets.\n")
