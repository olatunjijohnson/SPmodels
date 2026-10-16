args <- commandArgs(trailingOnly = TRUE)
replicates <- if (length(args) >= 1L) as.integer(args[[1L]]) else 5L
sizes <- if (length(args) >= 2L) {
  as.integer(strsplit(args[[2L]], ",", fixed = TRUE)[[1L]])
} else c(25L, 100L, 225L)
budgets <- if (length(args) >= 3L) {
  as.integer(strsplit(args[[3L]], ",", fixed = TRUE)[[1L]])
} else c(5L, 10L, 20L, 40L)
case_filter <- if (length(args) >= 4L) {
  strsplit(args[[4L]], ",", fixed = TRUE)[[1L]]
} else "gaussian"
archive_sha256 <- if (length(args) >= 5L) args[[5L]] else NA_character_
output_path <- if (length(args) >= 6L) args[[6L]] else
  paste0("docs/validation/data/neighbor_prediction_", case_filter, ".csv")
nsim <- if (length(args) >= 7L) as.integer(args[[7L]]) else 20000L
first_replicate <- if (length(args) >= 8L) as.integer(args[[8L]]) else 1L

if (is.na(replicates) || replicates < 1L || anyNA(sizes) ||
    !length(sizes) || any(!sizes %in% c(25L, 100L, 225L)) ||
    anyNA(budgets) || !length(budgets) || any(budgets < 1L) ||
    anyNA(nsim) || length(nsim) != 1L || nsim < 1L) {
  stop("Invalid replicate, grid-size, neighbor-budget or draw-count argument.",
       call. = FALSE)
}
valid_cases <- c("gaussian", "binomial", "poisson", "spacetime",
                 "grouped_slope")
if (length(case_filter) != 1L || any(!case_filter %in% valid_cases)) {
  stop("Run one NNP case per invocation; case must be a known case ID.",
       call. = FALSE)
}
if (is.na(first_replicate) || first_replicate < 1L ||
    first_replicate + replicates - 1L > 5L) {
  stop("Replicate indices must be within the registered range 1--5.",
       call. = FALSE)
}
partial_path <- paste0(output_path, ".partial")
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
if (first_replicate == 1L) {
  if (file.exists(output_path) || file.exists(partial_path)) {
    stop("Output or checkpoint already exists; choose another output path.",
         call. = FALSE)
  }
} else {
  if (!file.exists(partial_path) || file.exists(output_path)) {
    stop("A matching partial checkpoint is required to resume.",
         call. = FALSE)
  }
  checkpoint <- utils::read.csv(partial_path, stringsAsFactors = FALSE)
  if (nrow(checkpoint) == 0L ||
      any(checkpoint$case != case_filter) ||
      max(checkpoint$replicate) != first_replicate - 1L ||
      min(checkpoint$replicate) != 1L) {
    stop("The checkpoint does not match the requested case/replicate range.",
         call. = FALSE)
  }
}
for (file in list.files("R", full.names = TRUE)) source(file)

.sp_nnp_grid <- function(side) {
  expand.grid(east = seq(0, 5, length.out = side),
              north = seq(0, 5, length.out = side),
              KEEP.OUT.ATTRS = FALSE)
}

.sp_nnp_standardize <- function(x, center, scale) (x - center) / scale

.sp_nnp_capture <- function(expression) {
  captured <- character()
  value <- tryCatch(withCallingHandlers(expression, warning = function(warning) {
    captured <<- c(captured, conditionMessage(warning))
    invokeRestart("muffleWarning")
  }), error = identity)
  list(value = value, warnings = paste(unique(captured), collapse = " | "))
}

.sp_nnp_make_case <- function(case, seed) {
  smoothness <- 0.5
  spatial_variance <- 0.8
  scale <- 2
  nugget_variance <- 0.2
  temporal_scale <- 1.5
  fixed <- c(intercept = 0.3, slope = 0.4)
  slope_covariance <- matrix(c(0.3, 0.05, 0.05, 0.08), 2L)

  if (case %in% c("spacetime")) {
    xy <- expand.grid(east = 0:3, north = 0:3,
                      KEEP.OUT.ATTRS = FALSE)
    train <- xy[rep(seq_len(nrow(xy)), each = 3L), , drop = FALSE]
    train$time <- rep(c(0, 1, 2), times = nrow(xy))
  } else {
    train <- expand.grid(east = 0:5, north = 0:5,
                         KEEP.OUT.ATTRS = FALSE)
  }
  z_raw <- train$east + 0.5 * train$north
  z_center <- mean(z_raw)
  z_scale <- stats::sd(z_raw)
  train$z <- .sp_nnp_standardize(z_raw, z_center, z_scale)

  if (case == "grouped_slope") {
    train$group <- paste((train$east %/% 2) + 1L,
                         (train$north %/% 2) + 1L, sep = "_")
  }
  distance <- .sp_dist(as.matrix(train[c("east", "north")]))
  spatial_cor <- .sp_matern_cor_unchecked(distance, scale, smoothness)
  if (case == "spacetime") {
    time_lag <- abs(outer(train$time, train$time, "-"))
    correlation <- spatial_cor * exp(-time_lag / temporal_scale)
  } else {
    correlation <- spatial_cor
  }
  covariance <- spatial_variance * correlation
  diag(covariance) <- diag(covariance) + nugget_variance
  if (case == "grouped_slope") {
    design <- cbind(1, train$z)
    group_index <- match(train$group, unique(train$group))
    for (group in unique(group_index)) {
      index <- which(group_index == group)
      covariance[index, index] <- covariance[index, index] +
        design[index, , drop = FALSE] %*% slope_covariance %*%
        t(design[index, , drop = FALSE])
    }
  }

  set.seed(seed)
  latent <- as.vector(t(chol(covariance)) %*%
                       stats::rnorm(nrow(train)))
  eta <- fixed[[1L]] + fixed[[2L]] * train$z + latent
  family <- switch(case,
    gaussian = stats::gaussian(),
    binomial = stats::binomial(),
    poisson = stats::poisson(),
    spacetime = stats::gaussian(),
    grouped_slope = stats::gaussian()
  )
  train$y <- switch(case,
    gaussian = eta,
    binomial = stats::rbinom(nrow(train), 10L, stats::plogis(eta)),
    poisson = stats::rpois(nrow(train), exp(eta)),
    spacetime = eta,
    grouped_slope = eta
  )
  if (case == "binomial") train$trials <- 10L

  random_slope <- if (case == "grouped_slope") {
    sp_random_slope("group", "z", slope_covariance)
  } else NULL
  temporal <- if (case == "spacetime") {
    sp_temporal(~ time, temporal_scale, "years")
  } else NULL
  control <- sp_control(
    smoothness = smoothness, spatial_variance = spatial_variance,
    scale = scale, nugget_variance = nugget_variance,
    estimate = FALSE, distance_units = "unit"
  )
  fit_result <- .sp_nnp_capture(
    sp_fit(if (case == "binomial") cbind(y, trials - y) ~ z else y ~ z,
      train, family = family, coords = ~ east + north,
      time = temporal, random_slope = random_slope, control = control)
  )
  if (inherits(fit_result$value, "error")) stop(fit_result$value)
  list(fit = fit_result$value, training = train, case = case,
       fit_warnings = fit_result$warnings,
       z_center = z_center, z_scale = z_scale,
       slope_covariance = slope_covariance)
}

.sp_nnp_make_prediction <- function(case_data, side) {
  newdata <- .sp_nnp_grid(side)
  raw <- newdata$east + 0.5 * newdata$north
  newdata$z <- .sp_nnp_standardize(raw, case_data$z_center,
                                    case_data$z_scale)
  case <- case_data$case
  if (case == "spacetime") {
    newdata$time <- rep(c(0, 1, 2), length.out = nrow(newdata))
  }
  if (case == "grouped_slope") {
    train_xy <- as.matrix(case_data$training[c("east", "north")])
    new_xy <- as.matrix(newdata[c("east", "north")])
    nearest <- vapply(seq_len(nrow(new_xy)), function(i) {
      distance <- rowSums(sweep(train_xy, 2L, new_xy[i, ], "-")^2)
      which.min(distance)
    }, integer(1))
    newdata$group <- case_data$training$group[nearest]
  }
  midpoint <- 2.5
  quadrant <- ifelse(newdata$east < midpoint, "west", "east")
  quadrant <- paste(quadrant,
    ifelse(newdata$north < midpoint, "south", "north"), sep = "_")
  target_groups <- c("overall", "west_south", "west_north",
                     "east_south", "east_north")
  allocation_index <- rep(seq_len(nrow(newdata)), 5L)
  allocation_group <- rep(target_groups, each = nrow(newdata))
  keep <- allocation_group == "overall" |
    allocation_group == quadrant[rep(seq_len(nrow(newdata)), 5L)]
  allocation <- structure(list(
    index = allocation_index[keep], weights = rep(1, sum(keep)),
    groups = allocation_group[keep], point_count = nrow(newdata)
  ), class = "spmodels_polygon_weights")
  list(data = newdata, groups = target_groups, allocation = allocation)
}

.sp_nnp_target_metrics <- function(draws, threshold) {
  do.call(rbind, lapply(seq_len(ncol(draws)), function(column) {
    value <- draws[, column]
    deviation <- stats::sd(value)
    sd_influence <- if (is.finite(deviation) && deviation > 0) {
      ((value - mean(value))^2 - deviation^2) / (2 * deviation)
    } else rep(0, length(value))
    data.frame(target_id = colnames(draws)[[column]],
      mean = mean(value), sd = deviation,
      threshold_probability = mean(value > threshold),
      mean_mcse = deviation / sqrt(nrow(draws)),
      sd_mcse = stats::sd(sd_influence) / sqrt(nrow(draws)),
      threshold_mcse = sqrt(mean(value > threshold) *
        (1 - mean(value > threshold)) / nrow(draws)),
      stringsAsFactors = FALSE)
  }))
}

results <- list()
result_index <- 0L
for (case_index in match(case_filter, valid_cases)) {
  case <- valid_cases[[case_index]]
  for (replicate in seq.int(first_replicate,
                            length.out = replicates)) {
    replicate_start_index <- result_index + 1L
    response_seed <- 430000L + replicate
    case_data <- tryCatch(.sp_nnp_make_case(case, response_seed),
                          error = identity)
    for (side_count in sizes) {
      prediction <- if (inherits(case_data, "error")) NULL else
        .sp_nnp_make_prediction(case_data, as.integer(sqrt(side_count)))
      if (!is.null(prediction) && nrow(prediction$data) != side_count) {
        stop("Requested grid size is not a square side count.", call. = FALSE)
      }
      threshold <- switch(case, gaussian = 0, binomial = 0.5,
                          poisson = 1, spacetime = 0, grouped_slope = 0)
      type <- if (case %in% c("binomial", "poisson")) "response" else "link"
      exact <- NULL
      exact_target <- NULL
      exact_capture <- NULL
      exact_target_capture <- NULL
      exact_error <- ""
      exact_elapsed <- NA_real_
      fit_warnings <- if (inherits(case_data, "error"))
        conditionMessage(case_data) else case_data$fit_warnings
      if (!inherits(case_data, "error")) {
        exact_start <- proc.time()[["elapsed"]]
        exact_capture <- .sp_nnp_capture(sp_predict_draws(case_data$fit, prediction$data,
          nsim = nsim, target = "latent", type = type,
          seed = response_seed))
        exact <- exact_capture$value
        exact_elapsed <- proc.time()[["elapsed"]] - exact_start
        if (inherits(exact, "error")) {
          exact_error <- conditionMessage(exact)
          exact <- NULL
        } else {
          exact_target_capture <- .sp_nnp_capture(sp_predict_target(exact,
            weights = prediction$allocation, threshold = threshold))
          exact_target <- exact_target_capture$value
          if (inherits(exact_target, "error")) {
            exact_error <- conditionMessage(exact_target)
            exact_target <- NULL
          }
        }
      } else {
        exact_error <- conditionMessage(case_data)
      }
      prediction_warnings <- if (!is.null(exact_capture) &&
          nzchar(exact_capture$warnings)) exact_capture$warnings else ""
      if (!is.null(exact_target_capture) &&
          nzchar(exact_target_capture$warnings)) {
        prediction_warnings <- paste(prediction_warnings,
          exact_target_capture$warnings, sep = " | ")
      }
      exact_covariance <- if (is.null(exact)) NULL else
        exact$predictive_covariance
      exact_metrics <- if (is.null(exact_target)) NULL else
        .sp_nnp_target_metrics(exact_target$draws, threshold)

      orderings <- c("input", "maximin")
      config_budgets <- unique(pmin(budgets, side_count - 1L))
      for (ordering in orderings) for (budget in config_budgets) {
        covariance_result <- NULL
        covariance_elapsed <- NA_real_
        target_result <- NULL
        target_elapsed <- NA_real_
        errors <- character()
        if (is.null(exact) || is.null(exact_target)) {
          errors <- exact_error
        } else {
          covariance_start <- proc.time()[["elapsed"]]
          covariance_capture <- .sp_nnp_capture(sp_neighbor_covariance(
            case_data$fit, prediction$data, neighbors = budget,
            ordering = ordering, target = "latent", reference = TRUE
          ))
          covariance_result <- covariance_capture$value
          if (nzchar(covariance_capture$warnings)) {
            prediction_warnings <- paste(prediction_warnings,
              covariance_capture$warnings, sep = " | ")
          }
          covariance_elapsed <- proc.time()[["elapsed"]] - covariance_start
          if (inherits(covariance_result, "error")) {
            errors <- c(errors, conditionMessage(covariance_result))
            covariance_result <- NULL
          }
          target_start <- proc.time()[["elapsed"]]
          target_capture <- .sp_nnp_capture(sp_predict_target_neighbor(
            case_data$fit, prediction$data, weights = prediction$allocation,
            nsim = nsim, neighbors = budget, ordering = ordering,
            target = "latent", type = type, threshold = threshold,
            seed = response_seed
          ))
          target_result <- target_capture$value
          if (nzchar(target_capture$warnings)) {
            prediction_warnings <- paste(prediction_warnings,
              target_capture$warnings, sep = " | ")
          }
          target_elapsed <- proc.time()[["elapsed"]] - target_start
          if (inherits(target_result, "error")) {
            errors <- c(errors, conditionMessage(target_result))
            target_result <- NULL
          }
        }
        approximate_metrics <- if (is.null(target_result)) NULL else
          .sp_nnp_target_metrics(target_result$draws, threshold)
        exact_boundary_equal <- if (is.null(target_result) ||
            budget < side_count - 1L || ordering != "input") NA else
          isTRUE(all.equal(unname(target_result$draws),
            unname(exact_target$draws), tolerance = 2e-8,
            check.attributes = FALSE))
        covariance_errors <- if (is.null(covariance_result)) {
          c(relative_frobenius = NA_real_,
            maximum_absolute_covariance = NA_real_,
            maximum_absolute_marginal_variance = NA_real_)
        } else covariance_result$errors
        for (target_id in c("overall", "west_south", "west_north",
                            "east_south", "east_north")) {
          result_index <- result_index + 1L
          exact_row <- if (is.null(exact_metrics)) NULL else
            exact_metrics[exact_metrics$target_id == target_id, , drop = FALSE]
          approx_row <- if (is.null(approximate_metrics)) NULL else
            approximate_metrics[approximate_metrics$target_id == target_id,
                                , drop = FALSE]
          mean_difference <- if (is.null(exact_row) || is.null(approx_row))
            NA_real_ else approx_row$mean - exact_row$mean
          sd_difference <- if (is.null(exact_row) || is.null(approx_row))
            NA_real_ else approx_row$sd - exact_row$sd
          probability_difference <- if (is.null(exact_row) ||
              is.null(approx_row)) NA_real_ else
            approx_row$threshold_probability -
              exact_row$threshold_probability
          combined_mean_mcse <- if (is.null(exact_row) || is.null(approx_row))
            NA_real_ else sqrt(exact_row$mean_mcse^2 +
                                 approx_row$mean_mcse^2)
          combined_sd_mcse <- if (is.null(exact_row) || is.null(approx_row))
            NA_real_ else sqrt(exact_row$sd_mcse^2 +
                                 approx_row$sd_mcse^2)
          combined_probability_mcse <- if (is.null(exact_row) ||
              is.null(approx_row)) NA_real_ else sqrt(
                exact_row$threshold_mcse^2 + approx_row$threshold_mcse^2)
          results[[result_index]] <- data.frame(
            archive_sha256 = archive_sha256,
            r_version = as.character(getRversion()),
            platform = R.version$platform,
            operating_system = unname(Sys.info()[["sysname"]]),
            blas = unname(extSoftVersion()[["BLAS"]]),
            case = case, case_index = case_index,
            grid_size = side_count, replicate = replicate,
            response_seed = response_seed, nsim = nsim,
            fit_warnings = fit_warnings,
            fit_error = if (inherits(case_data, "error"))
              conditionMessage(case_data) else "",
            exact_prediction_error = exact_error,
            prediction_warnings = prediction_warnings,
            exact_elapsed_seconds = exact_elapsed,
            ordering = ordering, neighbors = budget,
            exact_boundary = budget >= side_count - 1L,
            exact_boundary_same_seed_equal = exact_boundary_equal,
            target_id = target_id, threshold = threshold, type = type,
            covariance_elapsed_seconds = covariance_elapsed,
            target_elapsed_seconds = target_elapsed,
            relative_frobenius_covariance_error =
              covariance_errors[["relative_frobenius"]],
            maximum_absolute_covariance_error =
              covariance_errors[["maximum_absolute_covariance"]],
            maximum_absolute_marginal_variance_error =
              covariance_errors[["maximum_absolute_marginal_variance"]],
            exact_target_mean = if (is.null(exact_row)) NA_real_ else
              exact_row$mean,
            approximate_target_mean = if (is.null(approx_row)) NA_real_ else
              approx_row$mean,
            target_mean_difference = mean_difference,
            combined_mean_mcse = combined_mean_mcse,
            exact_target_sd = if (is.null(exact_row)) NA_real_ else
              exact_row$sd,
            approximate_target_sd = if (is.null(approx_row)) NA_real_ else
              approx_row$sd,
            target_sd_difference = sd_difference,
            combined_sd_mcse = combined_sd_mcse,
            exact_threshold_probability = if (is.null(exact_row))
              NA_real_ else exact_row$threshold_probability,
            approximate_threshold_probability = if (is.null(approx_row))
              NA_real_ else approx_row$threshold_probability,
            threshold_probability_difference = probability_difference,
            combined_threshold_mcse = combined_probability_mcse,
            full_point_draw_bytes = if (is.null(exact)) NA_real_ else
              as.numeric(object.size(exact$draws)),
            streamed_target_draw_bytes = if (is.null(target_result))
              NA_real_ else as.numeric(object.size(target_result$draws)),
            peak_active_point_vectors = if (is.null(target_result))
              NA_integer_ else target_result$peak_active_locations,
            error = paste(unique(errors[nzchar(errors)]), collapse = " | "),
            stringsAsFactors = FALSE
          )
        }
        cat("completed case=", case, " grid=", side_count,
            " seed=", response_seed, " order=", ordering,
            " neighbors=", budget, "\n", sep = "")
      }
    }
    replicate_rows <- do.call(rbind,
      results[seq.int(replicate_start_index, result_index)])
    append_checkpoint <- file.exists(partial_path)
    utils::write.table(replicate_rows, partial_path, sep = ",", quote = TRUE,
      row.names = FALSE, col.names = !append_checkpoint,
      append = append_checkpoint, na = "")
    cat("checkpointed case=", case, " replicate=", replicate,
        " to ", partial_path, "\n", sep = "")
  }
}

if (!file.rename(partial_path, output_path)) {
  stop("Could not finalize NNP output from its checkpoint.", call. = FALSE)
}
cat("finalized ", file.info(output_path)$size, " bytes at ", output_path,
    "\n", sep = "")
