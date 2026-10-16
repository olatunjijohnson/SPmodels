args <- commandArgs(trailingOnly = TRUE)
archive_sha256 <- if (length(args) >= 1L) args[[1L]] else NA_character_
fit_path <- if (length(args) >= 2L) args[[2L]] else
  "docs/validation/data/mcml_proposal_study.csv"
fixture_path <- if (length(args) >= 3L) args[[3L]] else
  "docs/validation/data/mcml_proposal_fixtures.csv"
for (file in list.files("R", full.names = TRUE)) source(file)

cases <- list(
  list(id = "MCP-01", index = 1L, family = "binomial", width = 6L,
       height = 4L, seed = 510001L, trials = 8L,
       intercept = -0.4, slope = 0.3),
  list(id = "MCP-02", index = 2L, family = "poisson", width = 8L,
       height = 6L, seed = 510002L, trials = NA_integer_,
       intercept = 0.2, slope = 0.2)
)
policies <- list(
  list(id = "single_fixed", scales = 1, adapt = FALSE),
  list(id = "mixture_fixed", scales = c(0.75, 1.5), adapt = FALSE),
  list(id = "single_adaptive", scales = 1, adapt = TRUE),
  list(id = "mixture_adaptive", scales = c(0.75, 1.5), adapt = TRUE)
)
main_sizes <- c(500L, 1000L, 2000L)
attempts <- 5L

simulate_fixture <- function(case) {
  design <- expand.grid(east = seq_len(case$width) - 1L,
                       north = seq_len(case$height) - 1L)
  design$z <- as.numeric(scale(design$east + 0.5 * design$north))
  coords <- as.matrix(design[c("east", "north")])
  distance <- .sp_dist(coords)
  covariance <- 0.40 * exp(-distance / 2)
  diag(covariance) <- diag(covariance) + 0.20
  set.seed(case$seed)
  latent <- as.vector(t(chol(covariance)) %*% stats::rnorm(nrow(design)))
  eta <- case$intercept + case$slope * design$z + latent
  response <- if (case$family == "binomial") {
    stats::rbinom(nrow(design), case$trials, stats::plogis(eta))
  } else stats::rpois(nrow(design), exp(eta))
  design$y <- response
  design$trials <- if (case$family == "binomial") case$trials else NA_integer_
  design
}

fit_base <- function(case, data) {
  warnings <- character()
  error <- NA_character_
  model <- NULL
  elapsed <- system.time({
    model <- tryCatch(withCallingHandlers({
      formula <- if (case$family == "binomial") {
        stats::as.formula("cbind(y, trials - y) ~ z")
      } else y ~ z
      family <- if (case$family == "binomial") stats::binomial() else
        stats::poisson()
      sp_fit(formula, data, family = family, coords = ~ east + north,
        control = sp_control(
          smoothness = 0.5, spatial_variance = 0.4, scale = 1.5,
          nugget_variance = 0.2, estimate = TRUE,
          distance_units = "grid-unit", maxit = 200L
        ))
    }, warning = function(warning) {
      warnings <<- c(warnings, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }), error = function(condition) {
      error <<- conditionMessage(condition)
      NULL
    })
  })[["elapsed"]]
  list(model = model, elapsed = unname(elapsed), warnings = warnings,
       error = error)
}

fit_rows <- list()
fixture_rows <- list()
fixture_data <- list()
row_index <- 0L
for (case in cases) {
  data <- simulate_fixture(case)
  data$case <- case$id
  fixture_data[[case$index]] <- data
  base <- fit_base(case, data)
  fixture_rows[[case$index]] <- data.frame(
    archive_sha256 = archive_sha256, case = case$id,
    family = case$family, n = nrow(data), response_seed = case$seed,
    trials = case$trials, intercept_truth = case$intercept,
    slope_truth = case$slope,
    base_convergence = if (is.null(base$model)) NA_integer_ else
      base$model$optimisation$convergence,
    base_log_likelihood = if (is.null(base$model)) NA_real_ else
      as.numeric(logLik(base$model)), base_elapsed_seconds = base$elapsed,
    base_warnings = paste(base$warnings, collapse = " | "),
    base_error = base$error, stringsAsFactors = FALSE
  )
  base_usable <- !is.null(base$model) &&
    identical(base$model$optimisation$convergence, 0L) &&
    !length(base$warnings) && isTRUE(base$model$latent_mode$converged)
  if (is.null(base$model) || !base_usable) {
    for (policy_index in seq_along(policies)) for (size_index in seq_along(main_sizes)) {
      for (attempt in seq_len(attempts)) {
        row_index <- row_index + 1L
        seed <- 520001L + 100000L * case$index +
          10000L * policy_index + 100L * size_index + attempt
        fit_rows[[row_index]] <- data.frame(
          archive_sha256 = archive_sha256, case = case$id,
          case_index = case$index, family = case$family, n = nrow(data),
          policy = policy$id, policy_index = policy_index,
          proposal_scales = paste(policy$scales, collapse = ";"),
          adapt_requested = policy$adapt, adapt_used = NA,
          pilot_ess_fraction = NA_real_, pilot_relative_mcse = NA_real_,
          nsim = main_sizes[[size_index]], size_index = size_index,
          attempt = attempt, seed = seed, elapsed_seconds = NA_real_,
          usable_base_fit = FALSE, optimizer_convergence = NA_integer_,
          shared_parameter_log_likelihood = NA_real_,
          optimized_log_likelihood = NA_real_,
          effective_sample_size = NA_real_,
          effective_sample_fraction = NA_real_, maximum_weight = NA_real_,
          relative_mcse = NA_real_, estimate__intercept = NA_real_,
          estimate__z = NA_real_, estimate__spatial_variance = NA_real_,
          estimate__scale = NA_real_, estimate__nugget_variance = NA_real_,
          warnings = paste(base$warnings, collapse = " | "),
          error_message = if (is.null(base$model)) base$error else
            "Base fit failed the predeclared convergence/warning gate.",
          stringsAsFactors = FALSE
        )
      }
    }
    next
  }
  distance <- .sp_dist(base$model$coordinates)
  for (policy_index in seq_along(policies)) {
    policy <- policies[[policy_index]]
    for (size_index in seq_along(main_sizes)) {
      nsim <- main_sizes[[size_index]]
      for (attempt in seq_len(attempts)) {
        row_index <- row_index + 1L
        seed <- 520001L + 100000L * case$index +
          10000L * policy_index + 100L * size_index + attempt
        warnings <- character()
        error <- NA_character_
        model <- NULL
        elapsed <- system.time({
          model <- tryCatch(withCallingHandlers(
            sp_mcml(base$model,
              control = sp_mcml_control(
                nsim = nsim, proposal_scale = policy$scales,
                ess_warning = 0.10, adapt = policy$adapt,
                pilot_nsim = 500L, adapt_shrinkage = 0.5
              ), seed = seed
            ), warning = function(warning) {
              warnings <<- c(warnings, conditionMessage(warning))
              invokeRestart("muffleWarning")
            }
          ), error = function(condition) {
            error <<- conditionMessage(condition)
            NULL
          })
        })[["elapsed"]]
        diagnostics <- if (is.null(model)) NULL else model$mcml
        shared_objective <- if (is.null(model)) NA_real_ else tryCatch({
          weights <- .sp_mcml_log_weights(
            diagnostics$reference_parameters, base$model,
            diagnostics$latent_samples, diagnostics$log_proposal, distance
          )
          .sp_log_mean_exp(weights)
        }, error = function(condition) NA_real_)
        fit_rows[[row_index]] <- data.frame(
          archive_sha256 = archive_sha256, case = case$id,
          case_index = case$index, family = case$family, n = nrow(data),
          policy = policy$id, policy_index = policy_index,
          proposal_scales = paste(policy$scales, collapse = ";"),
          adapt_requested = policy$adapt, adapt_used = if (is.null(model))
            NA else isTRUE(model$mcml$adaptation$used),
          pilot_ess_fraction = if (is.null(model) ||
            is.null(model$mcml$adaptation)) NA_real_ else
              model$mcml$adaptation$pilot_diagnostics$effective_sample_fraction,
          pilot_relative_mcse = if (is.null(model) ||
            is.null(model$mcml$adaptation)) NA_real_ else
              model$mcml$adaptation$pilot_diagnostics$relative_mcse,
          nsim = nsim, size_index = size_index, attempt = attempt, seed = seed,
          elapsed_seconds = unname(elapsed), usable_base_fit = base_usable,
          optimizer_convergence = if (is.null(model)) NA_integer_ else
            model$optimisation$convergence,
          shared_parameter_log_likelihood = shared_objective,
          optimized_log_likelihood = if (is.null(model)) NA_real_ else
            model$log_likelihood,
          effective_sample_size = if (is.null(diagnostics)) NA_real_ else
            diagnostics$effective_sample_size,
          effective_sample_fraction = if (is.null(diagnostics)) NA_real_ else
            diagnostics$effective_sample_fraction,
          maximum_weight = if (is.null(diagnostics)) NA_real_ else
            diagnostics$maximum_weight,
          relative_mcse = if (is.null(diagnostics)) NA_real_ else
            diagnostics$relative_mcse,
          estimate__intercept = if (is.null(model)) NA_real_ else
            unname(model$coefficients[[1L]]),
          estimate__z = if (is.null(model)) NA_real_ else
            unname(model$coefficients[["z"]]),
          estimate__spatial_variance = if (is.null(model)) NA_real_ else
            unname(model$covariance[["spatial_variance"]]),
          estimate__scale = if (is.null(model)) NA_real_ else
            unname(model$covariance[["scale"]]),
          estimate__nugget_variance = if (is.null(model)) NA_real_ else
            unname(model$covariance[["nugget_variance"]]),
          warnings = paste(warnings, collapse = " | "),
          error_message = error, stringsAsFactors = FALSE
        )
        utils::write.csv(do.call(rbind, fit_rows), fit_path,
                         row.names = FALSE, na = "")
        cat(case$id, policy$id, "nsim", nsim, "attempt", attempt,
            "ESS", if (is.null(diagnostics)) NA else
              round(diagnostics$effective_sample_fraction, 4), "\n")
      }
    }
  }
}
utils::write.csv(do.call(rbind, fixture_data), fixture_path,
                 row.names = FALSE, na = "")
utils::write.csv(do.call(rbind, fit_rows), fit_path, row.names = FALSE,
                 na = "")
utils::write.csv(do.call(rbind, fixture_rows),
                 sub("\\.csv$", "_base.csv", fit_path),
                 row.names = FALSE, na = "")
cat("wrote", length(fit_rows), "MCML attempt rows\n")
