.sp_subset_prediction_design <- function(design, index) {
  answer <- design
  row_fields <- c("x", "coordinates", "data", "random_z", "random_slope",
                  "temporal_values", "intervention_multiplier", "offset")
  for (name in intersect(row_fields, names(answer))) {
    value <- answer[[name]]
    if (is.null(value)) next
    if (name == "random_slope") {
      answer[[name]] <- lapply(value, function(component) {
        component$group_values <- component$group_values[index]
        component$design <- component$design[index, , drop = FALSE]
        component
      })
    } else if (name == "data") {
      answer[[name]] <- value[index, , drop = FALSE]
    } else if (is.matrix(value)) {
      unseen_levels <- attr(value, "unseen_levels")
      answer[[name]] <- value[index, , drop = FALSE]
      if (!is.null(unseen_levels)) {
        attr(answer[[name]], "unseen_levels") <- lapply(
          unseen_levels, `[`, index
        )
      }
    } else {
      answer[[name]] <- value[index]
    }
  }
  answer
}

.sp_prediction_response_draws <- function(object, link_draws, design,
                                          target, type,
                                          intervention_multiplier_draws = NULL) {
  if (type == "link" || object$family$family == "gaussian") {
    return(link_draws)
  }
  if (object$family$family == "binomial") {
    probability <- stats::plogis(link_draws)
    if (!is.null(intervention_multiplier_draws)) {
      probability <- probability * intervention_multiplier_draws
    } else if (!is.null(design$intervention_multiplier)) {
      probability <- sweep(probability, 2L,
                           design$intervention_multiplier, "*")
    }
    if (target == "latent") return(probability)
    return(matrix(stats::rbinom(length(probability), 1L, probability),
                  nrow = nrow(probability)))
  }
  intensity <- exp(link_draws)
  if (target == "latent") return(intensity)
  matrix(stats::rpois(length(intensity), intensity),
         nrow = nrow(intensity))
}

.sp_maximin_prediction_order <- function(coordinates) {
  count <- nrow(coordinates)
  if (count <= 1L) return(seq_len(count))
  answer <- integer(count)
  selected <- rep(FALSE, count)
  minimum_distance <- rep(Inf, count)
  answer[[1L]] <- 1L
  selected[[1L]] <- TRUE
  for (step in 2:count) {
    latest <- answer[[step - 1L]]
    distance <- sqrt(rowSums(sweep(
      coordinates, 2L, coordinates[latest, ], "-"
    )^2))
    minimum_distance <- pmin(minimum_distance, distance)
    minimum_distance[selected] <- -Inf
    next_index <- which.max(minimum_distance)
    answer[[step]] <- next_index
    selected[[next_index]] <- TRUE
  }
  answer
}

.sp_psd_solve <- function(covariance, right) {
  covariance <- (covariance + t(covariance)) / 2
  root <- tryCatch(chol(covariance), error = function(e) NULL)
  if (!is.null(root)) return(.sp_solve_chol(root, right))
  decomposition <- eigen(covariance, symmetric = TRUE)
  scale <- max(1, max(abs(diag(covariance))))
  tolerance <- 1000 * .Machine$double.eps * scale * nrow(covariance)
  if (min(decomposition$values) < -tolerance) {
    stop("The predecessor covariance is not positive semidefinite.",
         call. = FALSE)
  }
  positive <- decomposition$values > tolerance
  solution <- if (any(positive)) {
    vectors <- decomposition$vectors[, positive, drop = FALSE]
    as.vector(vectors %*% (crossprod(vectors, right) /
      decomposition$values[positive]))
  } else rep(0, nrow(covariance))
  residual <- as.vector(covariance %*% solution - right)
  residual_tolerance <- sqrt(.Machine$double.eps) *
    max(1, max(abs(right))) * nrow(covariance)
  if (max(abs(residual)) > residual_tolerance) {
    stop("The current cross-covariance is incompatible with singular predecessors.",
         call. = FALSE)
  }
  solution
}

#' Bounded-memory approximate joint predictive draws
#'
#' Draw sequentially from nearest-predecessor conditionals of the fitted
#' Gaussian predictive distribution. This retains an auditable local
#' dependence approximation without constructing the complete prediction-site
#' covariance matrix. If every predecessor is retained, the factorization is
#' exact up to numerical precision.
#'
#' @param object A fitted spatial `spmodels_fit` object.
#' @param newdata Optional prediction data containing predictors and point
#'   coordinates. The fitted locations are used when this is `NULL`.
#' @param nsim Positive number of joint draws.
#' @param neighbors Positive maximum number of earlier prediction locations in
#'   each conditional distribution.
#' @param ordering Use input-row order or a deterministic maximin order that
#'   starts at row one and repeatedly selects the row farthest from the selected
#'   set. Ignored when `prediction_order` is supplied.
#' @param prediction_order Optional permutation of prediction-row positions.
#'   The input order is used by default and results are returned in input order.
#' @param target Draw the smooth latent target or future observation-site
#'   target.
#' @param type Draw on the response or linear-predictor scale.
#' @param seed Optional integer seed local to this call.
#' @param parameter_uncertainty Include the fitted local Gaussian parameter
#'   approximation? Each parameter draw rebuilds the fitted state and its
#'   complete nearest-predecessor factorization.
#' @param unseen Policy for grouping levels absent during fitting.
#'
#' @return An `spmodels_prediction_draws` object. Its full predictive covariance
#'   is `NULL`; `conditional_variance`, `conditional_neighbors` and
#'   `prediction_order` define the sequential approximation.
#' @export
sp_predict_draws_neighbor <- function(
    object, newdata = NULL, nsim = 1000L, neighbors = 30L,
    ordering = c("input", "maximin"), prediction_order = NULL,
    target = c("latent", "observation"),
    type = c("response", "link"), seed = NULL,
    parameter_uncertainty = FALSE,
    unseen = c("error", "marginal")) {
  if (!inherits(object, "spmodels_fit") || !isTRUE(object$spatial)) {
    stop("`object` must be a fitted spatial `spmodels_fit` model.",
         call. = FALSE)
  }
  integer_scalar <- function(value, name) {
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        !is.finite(value) || value < 1 || value != as.integer(value)) {
      stop(sprintf("`%s` must be one positive integer.", name), call. = FALSE)
    }
    as.integer(value)
  }
  nsim <- integer_scalar(nsim, "nsim")
  neighbors <- integer_scalar(neighbors, "neighbors")
  target <- match.arg(target)
  type <- match.arg(type)
  unseen <- match.arg(unseen)
  ordering <- match.arg(ordering)
  if (!is.logical(parameter_uncertainty) ||
      length(parameter_uncertainty) != 1L || is.na(parameter_uncertainty)) {
    stop("`parameter_uncertainty` must be TRUE or FALSE.", call. = FALSE)
  }
  design <- .sp_prediction_design(object, newdata, unseen)
  location_count <- nrow(design$x)
  if (!location_count) stop("Prediction data must contain at least one row.",
                            call. = FALSE)
  if (is.null(prediction_order)) {
    prediction_order <- if (ordering == "input") seq_len(location_count) else
      .sp_maximin_prediction_order(design$coordinates)
  } else {
    if (!is.numeric(prediction_order) || length(prediction_order) !=
        location_count || anyNA(prediction_order) ||
        any(!is.finite(prediction_order)) ||
        any(prediction_order != as.integer(prediction_order)) ||
        !setequal(as.integer(prediction_order), seq_len(location_count))) {
      stop("`prediction_order` must be a permutation of prediction-row positions.",
           call. = FALSE)
    }
    prediction_order <- as.integer(prediction_order)
    ordering <- "user"
  }

  conditional_neighbors <- vector("list", location_count)
  for (step in seq_len(location_count)) {
    current <- prediction_order[[step]]
    previous <- if (step == 1L) integer() else
      prediction_order[seq_len(step - 1L)]
    selected <- previous
    if (length(previous) > neighbors) {
      distance <- sqrt(rowSums(sweep(
        design$coordinates[previous, , drop = FALSE], 2L,
        design$coordinates[current, ], "-"
      )^2))
      selected <- previous[order(distance, seq_along(previous))[
        seq_len(neighbors)
      ]]
    }
    conditional_neighbors[[current]] <- selected
  }

  sequential_draws <- function(current_object, standard) {
    current_nsim <- nrow(standard)
    current_draws <- matrix(NA_real_, current_nsim, location_count)
    current_mean <- current_variance <- rep(NA_real_, location_count)
    for (step in seq_len(location_count)) {
      current <- prediction_order[[step]]
      selected <- conditional_neighbors[[current]]
      local_index <- c(selected, current)
      local_design <- .sp_subset_prediction_design(design, local_index)
      local_distribution <- if (current_object$family$family == "gaussian") {
        .sp_joint_gaussian(current_object, local_design, target)
      } else {
        .sp_joint_nongaussian(current_object, local_design, target)
      }
      local_count <- length(local_index)
      current_mean[[current]] <- local_distribution$mean[[local_count]]
      if (!length(selected)) {
        conditional_mean <- rep(local_distribution$mean[[1L]], current_nsim)
        variance <- local_distribution$covariance[1L, 1L]
      } else {
        neighbor_covariance <- local_distribution$covariance[
          -local_count, -local_count, drop = FALSE
        ]
        cross <- local_distribution$covariance[-local_count, local_count]
        coefficient <- tryCatch(
          .sp_psd_solve(neighbor_covariance, cross),
          error = function(e) stop(sprintf(
            "Predecessor covariance at prediction row %d is unusable: %s",
            current, conditionMessage(e)), call. = FALSE)
        )
        centered <- sweep(current_draws[, selected, drop = FALSE], 2L,
                          local_distribution$mean[-local_count], "-")
        conditional_mean <- local_distribution$mean[[local_count]] +
          as.vector(centered %*% coefficient)
        variance <- local_distribution$covariance[local_count, local_count] -
          sum(cross * coefficient)
      }
      scale <- max(1, abs(local_distribution$covariance[local_count,
                                                         local_count]))
      tolerance <- 1000 * .Machine$double.eps * scale * local_count
      if (!is.finite(variance) || variance < -tolerance) {
        stop(sprintf(
          "Prediction row %d has an invalid conditional variance.", current),
          call. = FALSE)
      }
      variance <- max(0, variance)
      current_variance[[current]] <- variance
      current_draws[, current] <- conditional_mean +
        sqrt(variance) * standard[, step]
    }
    list(draws = current_draws, mean = current_mean,
         variance = current_variance)
  }

  seed_scope <- .sp_seed_scope(seed)
  if (!is.null(seed_scope)) {
    on.exit(.sp_restore_seed(seed_scope), add = TRUE)
    set.seed(seed_scope$seed)
  }
  parameter_draws <- intervention_multiplier_draws <- NULL
  if (parameter_uncertainty) {
    specification <- .sp_parameter_specification(object)
    parameter_root <- .sp_covariance_root(specification$covariance)
    parameter_z <- matrix(stats::rnorm(nsim * length(specification$values)),
                          nrow = nsim)
    parameter_draws <- if (parameter_root$orientation == "upper") {
      sweep(parameter_z %*% parameter_root$root, 2L,
            specification$values, "+")
    } else {
      sweep(parameter_z %*% t(parameter_root$root), 2L,
            specification$values, "+")
    }
    link_draws <- matrix(NA_real_, nsim, location_count)
    predictive_mean <- conditional_variance <- rep(0, location_count)
    if (specification$intervention_count > 0L) {
      intervention_multiplier_draws <- matrix(NA_real_, nsim, location_count)
    }
    for (draw_index in seq_len(nsim)) {
      draw_object <- .sp_object_at_parameters(
        object, parameter_draws[draw_index, ], specification
      )
      if (is.null(draw_object)) {
        stop(sprintf("Parameter draw %d produced an unusable fitted state.",
                     draw_index), call. = FALSE)
      }
      current <- sequential_draws(
        draw_object,
        matrix(stats::rnorm(location_count), nrow = 1L)
      )
      link_draws[draw_index, ] <- current$draws
      predictive_mean <- predictive_mean + current$mean / nsim
      conditional_variance <- conditional_variance +
        current$variance / nsim
      if (!is.null(intervention_multiplier_draws)) {
        intervention_multiplier_draws[draw_index, ] <-
          .sp_intervention_multiplier(draw_object$intervention, design$data)
      }
    }
    predictive_mean <- colMeans(link_draws)
  } else {
    current <- sequential_draws(
      object, matrix(stats::rnorm(nsim * location_count), nrow = nsim)
    )
    link_draws <- current$draws
    predictive_mean <- current$mean
    conditional_variance <- current$variance
  }

  draws <- .sp_prediction_response_draws(
    object, link_draws, design, target, type,
    intervention_multiplier_draws = intervention_multiplier_draws
  )
  location_names <- row.names(if (is.null(newdata)) object$model else newdata)
  if (is.null(location_names)) {
    location_names <- paste0("location_", seq_len(location_count))
  }
  colnames(draws) <- location_names
  names(predictive_mean) <- names(conditional_variance) <- location_names
  names(conditional_neighbors) <- location_names
  structure(list(
    draws = draws, predictive_mean = predictive_mean,
    predictive_covariance = NULL,
    conditional_variance = conditional_variance,
    conditional_neighbors = conditional_neighbors,
    prediction_order = prediction_order,
    ordering = ordering,
    neighbors = min(neighbors, location_count - 1L),
    approximation = "nearest-predecessor Gaussian",
    parameter_draws = parameter_draws,
    intervention_multiplier_draws = intervention_multiplier_draws,
    parameter_uncertainty = parameter_uncertainty,
    family = object$family$family, target = target, type = type,
    seed = seed, nsim = nsim, call = match.call()
  ), class = "spmodels_prediction_draws")
}

#' Stream bounded-memory approximate predictive targets
#'
#' Generate nearest-predecessor predictive draws and aggregate declared
#' weighted targets immediately. Prediction-row draw vectors are retained only
#' until their final use in the neighbor graph.
#'
#' @param object A fitted spatial `spmodels_fit` object.
#' @param newdata Optional point prediction data.
#' @param weights Finite non-negative weights matching prediction rows, or an
#'   object returned by [sp_polygon_weights()] or [sp_polygon_area_weights()].
#' @param groups Optional target identifier per prediction row. It requires
#'   numeric `weights`; the default is one target.
#' @param nsim Positive number of draws.
#' @param neighbors Positive predecessor budget.
#' @param ordering Input or deterministic maximin order.
#' @param prediction_order Optional explicit permutation of prediction rows.
#' @param target Latent or future observation-site target.
#' @param type Response or link scale.
#' @param threshold Optional finite target exceedance threshold.
#' @param probs Finite target quantile probabilities in zero--one.
#' @param seed Optional integer seed local to the call.
#' @param unseen Policy for grouping levels absent during fitting.
#'
#' @return An `spmodels_prediction_target` containing target-level draws and
#'   summaries, sparse normalized weight records, neighbor audit metadata and
#'   the peak number of point-draw vectors retained simultaneously.
#' @export
sp_predict_target_neighbor <- function(
    object, newdata = NULL, weights, groups = NULL, nsim = 1000L,
    neighbors = 30L, ordering = c("input", "maximin"),
    prediction_order = NULL, target = c("latent", "observation"),
    type = c("response", "link"), threshold = NULL,
    probs = c(0.025, 0.5, 0.975), seed = NULL,
    unseen = c("error", "marginal")) {
  if (!inherits(object, "spmodels_fit") || !isTRUE(object$spatial)) {
    stop("`object` must be a fitted spatial `spmodels_fit` model.",
         call. = FALSE)
  }
  integer_scalar <- function(value, name) {
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        !is.finite(value) || value < 1 || value != as.integer(value)) {
      stop(sprintf("`%s` must be one positive integer.", name), call. = FALSE)
    }
    as.integer(value)
  }
  nsim <- integer_scalar(nsim, "nsim")
  neighbors <- integer_scalar(neighbors, "neighbors")
  ordering <- match.arg(ordering)
  target <- match.arg(target)
  type <- match.arg(type)
  unseen <- match.arg(unseen)
  if (!is.numeric(probs) || !length(probs) || anyNA(probs) ||
      any(!is.finite(probs)) || any(probs < 0 | probs > 1)) {
    stop("`probs` must contain finite probabilities between zero and one.",
         call. = FALSE)
  }
  probs <- sort(unique(probs))
  if (!is.null(threshold) &&
      (!is.numeric(threshold) || length(threshold) != 1L ||
       is.na(threshold) || !is.finite(threshold))) {
    stop("`threshold` must be NULL or one finite number.", call. = FALSE)
  }
  design <- .sp_prediction_design(object, newdata, unseen)
  location_count <- nrow(design$x)
  if (!location_count) stop("Prediction data must contain at least one row.",
                            call. = FALSE)

  polygon_weights <- NULL
  if (inherits(weights, "spmodels_polygon_weights")) {
    if (!is.null(groups)) {
      stop("Do not supply `groups` with polygon weights.", call. = FALSE)
    }
    if (weights$point_count != location_count) {
      stop("Polygon weights do not match the prediction locations.",
           call. = FALSE)
    }
    polygon_weights <- weights
    allocation_index <- weights$index
    allocation_weight <- weights$weights
    allocation_group <- as.character(weights$groups)
  } else {
    if (!is.numeric(weights) || length(weights) != location_count ||
        anyNA(weights) || any(!is.finite(weights)) || any(weights < 0)) {
      stop("`weights` must be finite non-negative values matching prediction locations.",
           call. = FALSE)
    }
    if (is.null(groups)) groups <- rep("all", location_count)
    if (length(groups) != location_count || anyNA(groups)) {
      stop("`groups` must be non-missing and match prediction locations.",
           call. = FALSE)
    }
    allocation_index <- seq_len(location_count)
    allocation_weight <- as.numeric(weights)
    allocation_group <- as.character(groups)
  }
  group_values <- unique(allocation_group)
  group_index <- match(allocation_group, group_values)
  totals <- as.numeric(rowsum(allocation_weight, group_index,
                              reorder = FALSE))
  if (length(totals) != length(group_values) || any(totals <= 0)) {
    stop("Every group must have positive total weight.", call. = FALSE)
  }
  normalized_weight <- allocation_weight / totals[group_index]
  allocation_by_location <- split(seq_along(allocation_index),
                                  allocation_index)

  if (is.null(prediction_order)) {
    prediction_order <- if (ordering == "input") seq_len(location_count) else
      .sp_maximin_prediction_order(design$coordinates)
  } else {
    if (!is.numeric(prediction_order) || length(prediction_order) !=
        location_count || anyNA(prediction_order) ||
        any(!is.finite(prediction_order)) ||
        any(prediction_order != as.integer(prediction_order)) ||
        !setequal(as.integer(prediction_order), seq_len(location_count))) {
      stop("`prediction_order` must be a permutation of prediction-row positions.",
           call. = FALSE)
    }
    prediction_order <- as.integer(prediction_order)
    ordering <- "user"
  }
  conditional_neighbors <- vector("list", location_count)
  position <- integer(location_count)
  position[prediction_order] <- seq_len(location_count)
  last_use <- position
  for (step in seq_len(location_count)) {
    current <- prediction_order[[step]]
    previous <- if (step == 1L) integer() else
      prediction_order[seq_len(step - 1L)]
    selected <- previous
    if (length(previous) > neighbors) {
      distance <- sqrt(rowSums(sweep(
        design$coordinates[previous, , drop = FALSE], 2L,
        design$coordinates[current, ], "-"
      )^2))
      selected <- previous[order(distance, seq_along(previous))[
        seq_len(neighbors)
      ]]
    }
    conditional_neighbors[[current]] <- selected
    if (length(selected)) last_use[selected] <- pmax(last_use[selected], step)
  }
  release_at <- split(seq_len(location_count), last_use)

  seed_scope <- .sp_seed_scope(seed)
  if (!is.null(seed_scope)) {
    on.exit(.sp_restore_seed(seed_scope), add = TRUE)
    set.seed(seed_scope$seed)
  }
  active <- new.env(parent = emptyenv(), hash = TRUE)
  peak_active <- 0L
  conditional_variance <- rep(NA_real_, location_count)
  target_draws <- matrix(0, nsim, length(group_values),
                         dimnames = list(NULL, group_values))
  for (step in seq_len(location_count)) {
    current <- prediction_order[[step]]
    selected <- conditional_neighbors[[current]]
    local_index <- c(selected, current)
    local_design <- .sp_subset_prediction_design(design, local_index)
    local_distribution <- if (object$family$family == "gaussian") {
      .sp_joint_gaussian(object, local_design, target)
    } else {
      .sp_joint_nongaussian(object, local_design, target)
    }
    local_count <- length(local_index)
    if (!length(selected)) {
      conditional_mean <- rep(local_distribution$mean[[1L]], nsim)
      variance <- local_distribution$covariance[1L, 1L]
    } else {
      neighbor_covariance <- local_distribution$covariance[
        -local_count, -local_count, drop = FALSE
      ]
      cross <- local_distribution$covariance[-local_count, local_count]
      coefficient <- .sp_psd_solve(neighbor_covariance, cross)
      predecessor_draws <- do.call(cbind, lapply(
        selected, function(index) get(as.character(index), active,
                                      inherits = FALSE)
      ))
      centered <- sweep(predecessor_draws, 2L,
                        local_distribution$mean[-local_count], "-")
      conditional_mean <- local_distribution$mean[[local_count]] +
        as.vector(centered %*% coefficient)
      variance <- local_distribution$covariance[local_count, local_count] -
        sum(cross * coefficient)
    }
    scale <- max(1, abs(local_distribution$covariance[local_count,
                                                       local_count]))
    tolerance <- 1000 * .Machine$double.eps * scale * local_count
    if (!is.finite(variance) || variance < -tolerance) {
      stop(sprintf("Prediction row %d has an invalid conditional variance.",
                   current), call. = FALSE)
    }
    variance <- max(0, variance)
    conditional_variance[[current]] <- variance
    link_draw <- conditional_mean + sqrt(variance) * stats::rnorm(nsim)
    assign(as.character(current), link_draw, active)
    peak_active <- max(peak_active, length(ls(active, all.names = TRUE)))

    point_draw <- link_draw
    if (type == "response" && object$family$family == "binomial") {
      point_draw <- stats::plogis(link_draw)
      if (!is.null(design$intervention_multiplier)) {
        point_draw <- point_draw * design$intervention_multiplier[[current]]
      }
      if (target == "observation") {
        point_draw <- stats::rbinom(nsim, 1L, point_draw)
      }
    } else if (type == "response" && object$family$family == "poisson") {
      point_draw <- exp(link_draw)
      if (target == "observation") point_draw <- stats::rpois(nsim, point_draw)
    }
    allocation <- allocation_by_location[[as.character(current)]]
    if (length(allocation)) {
      for (entry in allocation) {
        target_column <- group_index[[entry]]
        target_draws[, target_column] <- target_draws[, target_column] +
          normalized_weight[[entry]] * point_draw
      }
    }
    releasing <- release_at[[as.character(step)]]
    if (length(releasing)) {
      rm(list = as.character(releasing), envir = active)
    }
  }

  quantiles <- t(apply(target_draws, 2L, stats::quantile,
                       probs = probs, names = FALSE, type = 8))
  colnames(quantiles) <- paste0("q_", format(probs, trim = TRUE,
                                              scientific = FALSE))
  summary <- data.frame(
    target = group_values, mean = colMeans(target_draws),
    sd = apply(target_draws, 2L, stats::sd), row.names = NULL,
    check.names = FALSE
  )
  summary <- cbind(summary, quantiles)
  if (!is.null(threshold)) {
    summary$exceedance_probability <- colMeans(target_draws > threshold)
  }
  sparse_weights <- data.frame(
    location = allocation_index, target = allocation_group,
    weight = normalized_weight, stringsAsFactors = FALSE
  )
  structure(list(
    draws = target_draws, summary = summary, weights = sparse_weights,
    groups = group_values, threshold = threshold, probs = probs,
    polygon_weights = polygon_weights,
    family = object$family$family, target = target, type = type,
    approximation = "nearest-predecessor Gaussian",
    conditional_neighbors = conditional_neighbors,
    conditional_variance = conditional_variance,
    prediction_order = prediction_order, ordering = ordering,
    neighbors = min(neighbors, location_count - 1L),
    peak_active_locations = peak_active,
    seed = seed, nsim = nsim, call = match.call()
  ), class = "spmodels_prediction_target")
}

#' Inspect the covariance implied by nearest-predecessor prediction
#'
#' Reconstruct the complete covariance implied by the sequential conditional
#' approximation for a manageable prediction set. This is a diagnostic that
#' deliberately allocates a full matrix; use the draw or target-streaming
#' functions for bounded-memory production work.
#'
#' @param object A fitted spatial `spmodels_fit` object.
#' @param newdata Optional point prediction data.
#' @param neighbors Positive predecessor budget.
#' @param ordering Input or deterministic maximin order.
#' @param prediction_order Optional explicit permutation of prediction rows.
#' @param target Latent or future observation-site target.
#' @param unseen Policy for grouping levels absent during fitting.
#' @param reference Also construct the ordinary exact fitted predictive
#'   covariance and report approximation errors?
#'
#' @return An object of class `spmodels_neighbor_covariance` containing the
#'   implied covariance, predictive mean, conditional coefficients and
#'   variances, graph/order metadata and optional exact reference/errors.
#' @export
sp_neighbor_covariance <- function(
    object, newdata = NULL, neighbors = 30L,
    ordering = c("input", "maximin"), prediction_order = NULL,
    target = c("latent", "observation"),
    unseen = c("error", "marginal"), reference = FALSE) {
  if (!inherits(object, "spmodels_fit") || !isTRUE(object$spatial)) {
    stop("`object` must be a fitted spatial `spmodels_fit` model.",
         call. = FALSE)
  }
  if (!is.numeric(neighbors) || length(neighbors) != 1L ||
      is.na(neighbors) || !is.finite(neighbors) || neighbors < 1 ||
      neighbors != as.integer(neighbors)) {
    stop("`neighbors` must be one positive integer.", call. = FALSE)
  }
  neighbors <- as.integer(neighbors)
  ordering <- match.arg(ordering)
  target <- match.arg(target)
  unseen <- match.arg(unseen)
  if (!is.logical(reference) || length(reference) != 1L || is.na(reference)) {
    stop("`reference` must be TRUE or FALSE.", call. = FALSE)
  }
  design <- .sp_prediction_design(object, newdata, unseen)
  location_count <- nrow(design$x)
  if (!location_count) stop("Prediction data must contain at least one row.",
                            call. = FALSE)
  if (is.null(prediction_order)) {
    prediction_order <- if (ordering == "input") seq_len(location_count) else
      .sp_maximin_prediction_order(design$coordinates)
  } else {
    if (!is.numeric(prediction_order) || length(prediction_order) !=
        location_count || anyNA(prediction_order) ||
        any(!is.finite(prediction_order)) ||
        any(prediction_order != as.integer(prediction_order)) ||
        !setequal(as.integer(prediction_order), seq_len(location_count))) {
      stop("`prediction_order` must be a permutation of prediction-row positions.",
           call. = FALSE)
    }
    prediction_order <- as.integer(prediction_order)
    ordering <- "user"
  }

  covariance <- matrix(0, location_count, location_count)
  predictive_mean <- conditional_variance <- rep(NA_real_, location_count)
  conditional_neighbors <- conditional_coefficients <-
    vector("list", location_count)
  for (step in seq_len(location_count)) {
    current <- prediction_order[[step]]
    previous <- if (step == 1L) integer() else
      prediction_order[seq_len(step - 1L)]
    selected <- previous
    if (length(previous) > neighbors) {
      distance <- sqrt(rowSums(sweep(
        design$coordinates[previous, , drop = FALSE], 2L,
        design$coordinates[current, ], "-"
      )^2))
      selected <- previous[order(distance, seq_along(previous))[
        seq_len(neighbors)
      ]]
    }
    local_index <- c(selected, current)
    local_design <- .sp_subset_prediction_design(design, local_index)
    local_distribution <- if (object$family$family == "gaussian") {
      .sp_joint_gaussian(object, local_design, target)
    } else {
      .sp_joint_nongaussian(object, local_design, target)
    }
    local_count <- length(local_index)
    predictive_mean[[current]] <- local_distribution$mean[[local_count]]
    conditional_neighbors[[current]] <- selected
    if (!length(selected)) {
      coefficient <- numeric()
      variance <- local_distribution$covariance[1L, 1L]
      marginal_variance <- variance
    } else {
      neighbor_covariance <- local_distribution$covariance[
        -local_count, -local_count, drop = FALSE
      ]
      cross <- local_distribution$covariance[-local_count, local_count]
      coefficient <- .sp_psd_solve(neighbor_covariance, cross)
      variance <- local_distribution$covariance[local_count, local_count] -
        sum(cross * coefficient)
      approximate_cross <- as.vector(crossprod(
        coefficient, covariance[selected, previous, drop = FALSE]
      ))
      covariance[current, previous] <- approximate_cross
      covariance[previous, current] <- approximate_cross
      marginal_variance <- variance + sum(coefficient *
        (covariance[selected, selected, drop = FALSE] %*% coefficient))
    }
    scale <- max(1, abs(local_distribution$covariance[local_count,
                                                       local_count]))
    tolerance <- 1000 * .Machine$double.eps * scale * local_count
    if (!is.finite(variance) || variance < -tolerance ||
        !is.finite(marginal_variance) || marginal_variance < -tolerance) {
      stop(sprintf("Prediction row %d has an invalid implied variance.",
                   current), call. = FALSE)
    }
    variance <- max(0, variance)
    covariance[current, current] <- max(0, marginal_variance)
    conditional_variance[[current]] <- variance
    conditional_coefficients[[current]] <- coefficient
  }
  location_names <- row.names(if (is.null(newdata)) object$model else newdata)
  if (is.null(location_names)) {
    location_names <- paste0("location_", seq_len(location_count))
  }
  dimnames(covariance) <- list(location_names, location_names)
  names(predictive_mean) <- names(conditional_variance) <- location_names
  names(conditional_neighbors) <- names(conditional_coefficients) <-
    location_names

  exact <- errors <- NULL
  if (reference) {
    exact_distribution <- if (object$family$family == "gaussian") {
      .sp_joint_gaussian(object, design, target)
    } else {
      .sp_joint_nongaussian(object, design, target)
    }
    exact <- exact_distribution$covariance
    difference <- covariance - exact
    errors <- c(
      relative_frobenius = sqrt(sum(difference^2) / sum(exact^2)),
      maximum_absolute_covariance = max(abs(difference)),
      maximum_absolute_marginal_variance =
        max(abs(diag(covariance) - diag(exact)))
    )
  }
  structure(list(
    mean = predictive_mean, covariance = covariance,
    conditional_variance = conditional_variance,
    conditional_neighbors = conditional_neighbors,
    conditional_coefficients = conditional_coefficients,
    prediction_order = prediction_order, ordering = ordering,
    neighbors = min(neighbors, location_count - 1L), target = target,
    family = object$family$family, reference_covariance = exact,
    errors = errors, call = match.call()
  ), class = "spmodels_neighbor_covariance")
}

#' @export
print.spmodels_neighbor_covariance <- function(x, ...) {
  cat("Nearest-predecessor predictive covariance diagnostic\n")
  cat("Locations:", nrow(x$covariance), " Neighbors:", x$neighbors,
      " Ordering:", x$ordering, "\n")
  if (!is.null(x$errors)) print(x$errors)
  invisible(x)
}

