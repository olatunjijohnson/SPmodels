args <- commandArgs(trailingOnly = TRUE)
data_path <- if (length(args) >= 1L) args[[1L]] else stop(
  "Supply the authorised lsda 0.2.0 loaloa.rda path.", call. = FALSE
)
archive_sha256 <- if (length(args) >= 2L) args[[2L]] else NA_character_
output_path <- if (length(args) >= 3L) args[[3L]] else
  "docs/validation/data/mcml_proposal_loa.csv"
for (file in list.files("R", full.names = TRUE)) source(file)

loaded <- new.env(parent = globalenv())
load(data_path, envir = loaded)
if (!exists("loaloa", envir = loaded, inherits = FALSE)) {
  stop("The supplied lsda file does not contain `loaloa`.", call. = FALSE)
}
loaloa <- get("loaloa", envir = loaded)
loaloa$elevation_100m <- loaloa$elevation / 100
loaloa$mean_ndvi <- loaloa$mean9901
loaloa_sf <- sf::st_as_sf(
  loaloa, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE
)
loaloa_laea <- sf::st_transform(
  loaloa_sf,
  "+proj=laea +lat_0=7 +lon_0=12 +datum=WGS84 +units=km +no_defs"
)

base_warnings <- character()
base_error <- NA_character_
base_fit <- tryCatch(withCallingHandlers(sp_fit(
  cbind(no_inf, no_exam - no_inf) ~ elevation_100m + mean_ndvi,
  data = loaloa_laea, family = stats::binomial(),
  control = sp_control(
    smoothness = 0.5, spatial_variance = 0.5, scale = 150,
    nugget_variance = 0.5, distance_units = "km", maxit = 120L
  )
), warning = function(warning) {
  base_warnings <<- c(base_warnings, conditionMessage(warning))
  invokeRestart("muffleWarning")
}), error = function(condition) {
  base_error <<- conditionMessage(condition)
  NULL
})
base_usable <- !is.null(base_fit) &&
  identical(base_fit$optimisation$convergence, 0L) &&
  !length(base_warnings) && isTRUE(base_fit$latent_mode$converged)
if (!base_usable) {
  stop("The checked teaching fit failed its convergence/warning gate; no MCML attempt was run.",
       call. = FALSE)
}

policies <- list(
  list(id = "single_fixed", scales = 1, adapt = FALSE),
  list(id = "mixture_fixed", scales = c(0.75, 1.5), adapt = FALSE),
  list(id = "single_adaptive", scales = 1, adapt = TRUE),
  list(id = "mixture_adaptive", scales = c(0.75, 1.5), adapt = TRUE)
)
rows <- list()
index <- 0L
distance <- .sp_dist(base_fit$coordinates)
for (policy_index in seq_along(policies)) {
  policy <- policies[[policy_index]]
  for (attempt in 1:5) {
    index <- index + 1L
    seed <- 730001L + 10000L * policy_index + attempt
    warnings <- character()
    error <- NA_character_
    fit <- NULL
    elapsed <- system.time({
      fit <- tryCatch(withCallingHandlers(sp_mcml(
        base_fit,
        control = sp_mcml_control(
          nsim = 2000L, proposal_scale = policy$scales,
          ess_warning = 0.10, adapt = policy$adapt,
          pilot_nsim = 500L, adapt_shrinkage = 0.5
        ), seed = seed
      ), warning = function(warning) {
        warnings <<- c(warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }), error = function(condition) {
        error <<- conditionMessage(condition)
        NULL
      })
    })[["elapsed"]]
    mcml <- if (is.null(fit)) NULL else fit$mcml
    shared_objective <- if (is.null(fit)) NA_real_ else tryCatch(
      .sp_log_mean_exp(.sp_mcml_log_weights(
        mcml$reference_parameters, base_fit, mcml$latent_samples,
        mcml$log_proposal, distance
      )), error = function(condition) NA_real_
    )
    rows[[index]] <- data.frame(
      archive_sha256 = archive_sha256,
      input_sha256 = "0acc9426c9842f93c264914a86e62641c160bdcdf134391cc57be7facb3c4811",
      lsda_version = "0.2.0", case = "MCP-03-Loa-loa",
      policy = policy$id, policy_index = policy_index,
      proposal_scales = paste(policy$scales, collapse = ";"),
      adapt_requested = policy$adapt,
      adapt_used = if (is.null(fit)) NA else
        isTRUE(fit$mcml$adaptation$used),
      pilot_ess_fraction = if (is.null(fit) ||
        is.null(fit$mcml$adaptation)) NA_real_ else
          fit$mcml$adaptation$pilot_diagnostics$effective_sample_fraction,
      nsim = 2000L, attempt = attempt, seed = seed,
      base_log_likelihood = as.numeric(logLik(base_fit)),
      shared_parameter_log_likelihood = shared_objective,
      optimized_log_likelihood = if (is.null(fit)) NA_real_ else
        fit$log_likelihood,
      optimizer_convergence = if (is.null(fit)) NA_integer_ else
        fit$optimisation$convergence,
      effective_sample_fraction = if (is.null(mcml)) NA_real_ else
        mcml$effective_sample_fraction,
      maximum_weight = if (is.null(mcml)) NA_real_ else mcml$maximum_weight,
      relative_mcse = if (is.null(mcml)) NA_real_ else mcml$relative_mcse,
      elapsed_seconds = unname(elapsed),
      warnings = paste(warnings, collapse = " | "),
      error_message = error, stringsAsFactors = FALSE
    )
    utils::write.csv(do.call(rbind, rows), output_path,
                     row.names = FALSE, na = "")
    cat(policy$id, "attempt", attempt, "ESS",
        if (is.null(mcml)) NA else round(mcml$effective_sample_fraction, 4),
        "\n")
  }
}
cat("wrote", length(rows), "descriptive Loa loa MCML attempts\n")
