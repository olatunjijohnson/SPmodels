.sp_has_fork_workers <- function() {
  platform <- paste(R.version$platform, Sys.info()[["sysname"]])
  .Platform$OS.type != "windows" &&
    !grepl("emscripten|wasm", platform, ignore.case = TRUE)
}

.sp_resolve_folds <- function(folds, n, seed) {
  if (is.numeric(folds) && length(folds) == 1L && !is.na(folds) &&
      is.finite(folds) && folds == as.integer(folds)) {
    k <- as.integer(folds)
    if (k < 2L || k > n) stop("A fold count must be between 2 and the number of observations.", call. = FALSE)
    seed_scope <- .sp_seed_scope(seed)
    if (!is.null(seed_scope)) {
      on.exit(.sp_restore_seed(seed_scope), add = TRUE)
      set.seed(seed_scope$seed)
    }
    return(sample(rep(seq_len(k), length.out = n)))
  }
  if (length(folds) != n || anyNA(folds)) {
    stop("`folds` must be a fold count or one non-missing identifier per observation.", call. = FALSE)
  }
  labels <- unique(folds)
  if (length(labels) < 2L) stop("At least two folds are required.", call. = FALSE)
  folds
}

#' Integrated pointwise predictive log scores
#'
#' Evaluate negative log predictive densities at held-out outcomes. Gaussian
#' scores use their predictive normal distribution. Binomial and Poisson
#' scores integrate the response probability over the approximate Gaussian
#' predictive distribution on the link scale.
#'
#' @param object A fitted `spmodels_fit` model.
#' @param newdata Prediction data containing all required predictors,
#'   coordinates and optional time/intervention fields.
#' @param observed One finite observed response per prediction row. Binomial
#'   values are success counts.
#' @param trials Optional non-negative integer binomial trial totals. It must
#'   be omitted for other families and defaults to one for binomial models.
#' @param unseen Policy for grouping levels absent from a training fit.
#'
#' @return A numeric vector of negative log predictive densities.
#' @export
sp_predictive_log_score <- function(object, newdata, observed, trials = NULL,
                                    unseen = c("error", "marginal")) {
  if (!inherits(object, "spmodels_fit")) {
    stop("`object` must be a fitted `spmodels_fit` model.", call. = FALSE)
  }
  if (!is.data.frame(newdata)) {
    stop("`newdata` must be a data frame.", call. = FALSE)
  }
  if (!is.numeric(observed) || length(observed) != nrow(newdata) ||
      anyNA(observed) || any(!is.finite(observed))) {
    stop("`observed` must contain one finite numeric value per prediction row.",
         call. = FALSE)
  }
  unseen <- match.arg(unseen)
  family <- object$family$family
  if (family == "gaussian") {
    if (!is.null(trials)) stop("`trials` is only used for binomial models.",
                               call. = FALSE)
    prediction <- stats::predict(object, newdata, se.fit = TRUE,
      target = "observation", type = "response", unseen = unseen)
    return(-stats::dnorm(observed, prediction$fit, prediction$se.fit,
                         log = TRUE))
  }
  if (family == "binomial") {
    if (is.null(trials)) trials <- rep(1, length(observed))
    if (!is.numeric(trials) || length(trials) != length(observed) ||
        anyNA(trials) || any(!is.finite(trials)) || any(trials < 0) ||
        any(trials != floor(trials)) || any(observed < 0) ||
        any(observed != floor(observed)) || any(observed > trials)) {
      stop("Binomial `observed` and `trials` must be compatible non-negative integer counts.",
           call. = FALSE)
    }
  } else {
    if (!is.null(trials)) stop("`trials` is only used for binomial models.",
                               call. = FALSE)
    if (any(observed < 0) || any(observed != floor(observed))) {
      stop("Poisson `observed` values must be non-negative integer counts.",
           call. = FALSE)
    }
  }
  link <- stats::predict(object, newdata, se.fit = TRUE,
    target = "observation", type = "link", unseen = unseen)
  multiplier <- if (family == "binomial") {
    value <- .sp_intervention_multiplier(object$intervention, newdata)
    if (is.null(value)) rep(1, length(observed)) else value
  } else NULL
  log_probability <- vapply(seq_along(observed), function(i) {
    log_response_probability <- function(eta) {
      if (family == "binomial") {
        return(stats::dbinom(observed[[i]], trials[[i]],
          multiplier[[i]] * stats::plogis(eta), log = TRUE))
      }
      answer <- rep(-Inf, length(eta))
      usable <- eta < log(.Machine$double.xmax)
      answer[usable] <- observed[[i]] * eta[usable] - exp(eta[usable]) -
        lgamma(observed[[i]] + 1)
      answer
    }
    if (link$se.fit[[i]] <= 100 * .Machine$double.eps) {
      return(log_response_probability(link$fit[[i]]))
    }
    log_integrand <- function(z) {
      eta <- link$fit[[i]] + link$se.fit[[i]] * z
      log_response_probability(eta) + stats::dnorm(z, log = TRUE)
    }
    centre_bounds <- (-link$fit[[i]] + c(-20, 20)) / link$se.fit[[i]]
    search_bounds <- range(c(-20, 20, centre_bounds))
    optimum <- tryCatch(stats::optimize(
      function(z) {
        value <- -log_integrand(z)
        value[!is.finite(value)] <- .Machine$double.xmax
        value
      }, search_bounds
    ), error = function(e) NULL)
    if (is.null(optimum) || !is.finite(optimum$objective)) return(NA_real_)
    maximum <- -optimum$objective
    scaled <- tryCatch(stats::integrate(
      function(z) exp(log_integrand(z) - maximum),
      optimum$minimum - 12, optimum$minimum + 12,
      rel.tol = 1e-8, subdivisions = 200L
    )$value, error = function(e) NA_real_)
    if (!is.finite(scaled) || scaled <= 0) return(NA_real_)
    maximum + log(scaled)
  }, numeric(1))
  if (anyNA(log_probability) || any(!is.finite(log_probability))) {
    stop("Predictive probability integration failed or underflowed.",
         call. = FALSE)
  }
  -log_probability
}

#' Construct geometry-driven spatial folds
#'
#' Cluster fitted coordinate locations into spatially compact validation folds.
#' Repeated coordinates are assigned together so one physical location cannot
#' leak across training and held-out sets.
#'
#' @param object A spatial `spmodels_fit` object.
#' @param folds Number of spatial folds.
#' @param seed Optional local seed for k-means initialization.
#' @param nstart Positive number of random k-means starts.
#'
#' @return An integer vector of fold identifiers in fitted-row order, with
#'   class `spmodels_spatial_folds` and canonical cluster centres attached.
#' @export
sp_spatial_folds <- function(object, folds = 5L, seed = NULL,
                             nstart = 20L) {
  if (!inherits(object, "spmodels_fit") || !isTRUE(object$spatial) ||
      is.null(object$coordinates)) {
    stop("`object` must be a fitted spatial `spmodels_fit` model.",
         call. = FALSE)
  }
  integer_scalar <- function(value, name, minimum) {
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        !is.finite(value) || value < minimum || value != as.integer(value)) {
      stop(sprintf("`%s` must be one integer not less than %d.",
                   name, minimum), call. = FALSE)
    }
    as.integer(value)
  }
  folds <- integer_scalar(folds, "folds", 2L)
  nstart <- integer_scalar(nstart, "nstart", 1L)
  coordinates <- object$coordinates
  key <- paste(format(coordinates[, 1L], digits = 17, scientific = TRUE),
               format(coordinates[, 2L], digits = 17, scientific = TRUE),
               sep = "\r")
  first <- !duplicated(key)
  unique_coordinates <- coordinates[first, , drop = FALSE]
  unique_key <- key[first]
  if (folds > nrow(unique_coordinates)) {
    stop("`folds` cannot exceed the number of unique coordinate locations.",
         call. = FALSE)
  }
  seed_scope <- .sp_seed_scope(seed)
  if (!is.null(seed_scope)) {
    on.exit(.sp_restore_seed(seed_scope), add = TRUE)
    set.seed(seed_scope$seed)
  }
  if (folds == nrow(unique_coordinates)) {
    cluster <- seq_len(folds)
    centres <- unique_coordinates
  } else {
    clustering <- stats::kmeans(unique_coordinates, centers = folds,
                                nstart = nstart)
    cluster <- clustering$cluster
    centres <- clustering$centers
  }
  canonical_order <- order(-centres[, 2L], centres[, 1L])
  canonical_label <- integer(folds)
  canonical_label[canonical_order] <- seq_len(folds)
  unique_fold <- canonical_label[cluster]
  resolved <- unique_fold[match(key, unique_key)]
  canonical_centres <- centres[canonical_order, , drop = FALSE]
  row.names(canonical_centres) <- as.character(seq_len(folds))
  structure(as.integer(resolved), class = c("spmodels_spatial_folds", "integer"),
            centres = canonical_centres, seed = seed, nstart = nstart,
            call = match.call())
}

#' @export
print.spmodels_spatial_folds <- function(x, ...) {
  cat("SPmodels geometry-driven spatial folds\n")
  print(table(fold = as.integer(x)))
  invisible(x)
}

#' Construct regular spatial-block folds
#'
#' Partition fitted coordinates into axis-aligned blocks, then allocate whole
#' blocks to validation folds. No block or repeated coordinate can be split
#' between training and test data.
#'
#' @param object A spatial `spmodels_fit` object.
#' @param block_size One positive block width, or widths for the two coordinate
#'   axes, in fitted distance units.
#' @param folds Number of validation folds, not exceeding occupied blocks.
#' @param origin Optional two-coordinate grid origin. The coordinate-wise
#'   minima are used by default.
#' @param assignment Allocate canonically ordered blocks round-robin, or
#'   randomly permute blocks before balanced allocation.
#' @param seed Optional local seed used for random block allocation.
#'
#' @return An integer vector in fitted-row order with block metadata and class
#'   `spmodels_block_folds`.
#' @export
sp_block_folds <- function(object, block_size, folds = 5L, origin = NULL,
                           assignment = c("systematic", "random"),
                           seed = NULL) {
  if (!inherits(object, "spmodels_fit") || !isTRUE(object$spatial) ||
      is.null(object$coordinates)) {
    stop("`object` must be a fitted spatial `spmodels_fit` model.",
         call. = FALSE)
  }
  if (!is.numeric(block_size) || !length(block_size) ||
      length(block_size) > 2L || anyNA(block_size) ||
      any(!is.finite(block_size)) || any(block_size <= 0)) {
    stop("`block_size` must contain one or two finite positive widths.",
         call. = FALSE)
  }
  block_size <- rep(as.numeric(block_size), length.out = 2L)
  coordinates <- object$coordinates
  if (is.null(origin)) origin <- apply(coordinates, 2L, min)
  if (!is.numeric(origin) || length(origin) != 2L || anyNA(origin) ||
      any(!is.finite(origin))) {
    stop("`origin` must be NULL or two finite coordinates.", call. = FALSE)
  }
  if (!is.numeric(folds) || length(folds) != 1L || is.na(folds) ||
      !is.finite(folds) || folds < 2L || folds != as.integer(folds)) {
    stop("`folds` must be one integer not less than 2.", call. = FALSE)
  }
  assignment <- match.arg(assignment)
  block_index <- floor(sweep(coordinates, 2L, origin, "-") /
                         rep(block_size, each = nrow(coordinates)))
  key <- paste(block_index[, 1L], block_index[, 2L], sep = "\r")
  unique_index <- !duplicated(key)
  blocks <- block_index[unique_index, , drop = FALSE]
  block_key <- key[unique_index]
  folds <- as.integer(folds)
  if (folds > nrow(blocks)) {
    stop("`folds` cannot exceed the number of occupied spatial blocks.",
         call. = FALSE)
  }
  canonical <- order(-blocks[, 2L], blocks[, 1L])
  allocation_order <- canonical
  if (assignment == "random") {
    seed_scope <- .sp_seed_scope(seed)
    if (!is.null(seed_scope)) {
      on.exit(.sp_restore_seed(seed_scope), add = TRUE)
      set.seed(seed_scope$seed)
    }
    allocation_order <- sample(canonical)
  }
  block_fold <- integer(nrow(blocks))
  block_fold[allocation_order] <- rep(seq_len(folds),
                                      length.out = nrow(blocks))
  resolved <- block_fold[match(key, block_key)]
  metadata <- data.frame(
    block_x = blocks[, 1L], block_y = blocks[, 2L], fold = block_fold,
    rows = as.integer(tabulate(match(key, block_key), nbins = nrow(blocks)))
  )
  metadata <- metadata[canonical, , drop = FALSE]
  row.names(metadata) <- NULL
  structure(as.integer(resolved),
    class = c("spmodels_block_folds", "integer"), blocks = metadata,
    block_size = block_size, origin = as.numeric(origin),
    assignment = assignment, seed = seed,
    distance_units = object$distance_units, call = match.call())
}

#' @export
print.spmodels_block_folds <- function(x, ...) {
  cat("SPmodels regular spatial-block folds\n")
  cat("Block size:", paste(format(attr(x, "block_size")), collapse = " x "),
      attr(x, "distance_units"), " Assignment:", attr(x, "assignment"), "\n")
  print(table(fold = as.integer(x)))
  invisible(x)
}

#' Construct group-constrained validation folds
#'
#' Allocate every level of a declared fitted-data grouping variable wholly to
#' one validation fold. A greedy allocation balances row counts while
#' preventing group leakage.
#'
#' @param object A fitted `spmodels_fit` object.
#' @param group Name of one grouping column in the fitted data.
#' @param folds Number of folds, not exceeding the number of observed groups.
#' @param assignment Deterministic size-first allocation, or a locally seeded
#'   randomized group order followed by the same least-loaded allocation.
#' @param seed Optional local seed for randomized assignment.
#'
#' @return An integer fitted-row fold vector of class
#'   `spmodels_group_folds`, with group allocation metadata.
#' @export
sp_group_folds <- function(object, group, folds = 5L,
                           assignment = c("balanced", "random"),
                           seed = NULL) {
  if (!inherits(object, "spmodels_fit")) {
    stop("`object` must be a fitted `spmodels_fit` model.", call. = FALSE)
  }
  if (!is.character(group) || length(group) != 1L || is.na(group) ||
      !nzchar(group) || !group %in% names(object$data_used)) {
    stop("`group` must name one column in the fitted data.", call. = FALSE)
  }
  value <- object$data_used[[group]]
  if (is.list(value) || is.matrix(value) || anyNA(value)) {
    stop("The grouping column must be atomic and contain no missing values.",
         call. = FALSE)
  }
  label <- as.character(value)
  if (any(!nzchar(label))) {
    stop("Group labels must be non-empty.", call. = FALSE)
  }
  levels_seen <- unique(label)
  if (!is.numeric(folds) || length(folds) != 1L || is.na(folds) ||
      !is.finite(folds) || folds < 2L || folds != as.integer(folds)) {
    stop("`folds` must be one integer not less than 2.", call. = FALSE)
  }
  folds <- as.integer(folds)
  if (folds > length(levels_seen)) {
    stop("`folds` cannot exceed the number of observed groups.", call. = FALSE)
  }
  assignment <- match.arg(assignment)
  group_index <- match(label, levels_seen)
  group_size <- tabulate(group_index, nbins = length(levels_seen))
  order_groups <- order(-group_size, seq_along(levels_seen))
  if (assignment == "random") {
    seed_scope <- .sp_seed_scope(seed)
    if (!is.null(seed_scope)) {
      on.exit(.sp_restore_seed(seed_scope), add = TRUE)
      set.seed(seed_scope$seed)
    }
    order_groups <- sample(seq_along(levels_seen))
  }
  group_fold <- integer(length(levels_seen))
  fold_rows <- integer(folds)
  for (index in order_groups) {
    destination <- which.min(fold_rows)
    group_fold[[index]] <- destination
    fold_rows[[destination]] <- fold_rows[[destination]] + group_size[[index]]
  }
  metadata <- data.frame(
    group = levels_seen, rows = group_size, fold = group_fold,
    stringsAsFactors = FALSE
  )
  structure(as.integer(group_fold[group_index]),
    class = c("spmodels_group_folds", "integer"), groups = metadata,
    group = group, assignment = assignment, seed = seed,
    call = match.call())
}

#' @export
print.spmodels_group_folds <- function(x, ...) {
  cat("SPmodels group-constrained validation folds\n")
  cat("Grouping variable:", attr(x, "group"),
      " Assignment:", attr(x, "assignment"), "\n")
  print(table(fold = as.integer(x)))
  invisible(x)
}

#' Construct environmental-space validation folds
#'
#' Cluster unique fitted locations by declared numeric environmental variables.
#' Repeated coordinates stay together and must have consistent environmental
#' values.
#'
#' @param object A spatial `spmodels_fit` object.
#' @param variables Character names or a one-sided additive formula selecting
#'   numeric environmental variables from the fitted data.
#' @param folds Number of environmental folds.
#' @param scale Standardize variables before clustering?
#' @param seed Optional local k-means seed.
#' @param nstart Positive number of k-means starts.
#'
#' @return An integer fitted-row fold vector of class
#'   `spmodels_environmental_folds` with clustering metadata.
#' @export
sp_environmental_folds <- function(object, variables, folds = 5L,
                                   scale = TRUE, seed = NULL,
                                   nstart = 20L) {
  if (!inherits(object, "spmodels_fit") || !isTRUE(object$spatial) ||
      is.null(object$coordinates)) {
    stop("`object` must be a fitted spatial `spmodels_fit` model.",
         call. = FALSE)
  }
  if (inherits(variables, "formula")) {
    if (length(variables) != 2L) {
      stop("`variables` formula must be one-sided.", call. = FALSE)
    }
    labels <- attr(stats::terms(variables), "term.labels")
    variables <- all.vars(variables)
    if (!length(labels) || !identical(labels, variables)) {
      stop("`variables` supports names joined by `+`, without transformations.",
           call. = FALSE)
    }
  }
  if (!is.character(variables) || !length(variables) || anyNA(variables) ||
      any(!nzchar(variables)) || anyDuplicated(variables)) {
    stop("`variables` must select one or more unique named columns.",
         call. = FALSE)
  }
  missing <- setdiff(variables, names(object$data_used))
  if (length(missing)) {
    stop(sprintf("Environmental variables are missing: %s.",
                 paste(missing, collapse = ", ")), call. = FALSE)
  }
  environmental <- object$data_used[, variables, drop = FALSE]
  if (inherits(environmental, "sf")) environmental <- sf::st_drop_geometry(environmental)
  if (any(!vapply(environmental, is.numeric, logical(1))) ||
      anyNA(environmental) || any(!is.finite(as.matrix(environmental)))) {
    stop("Environmental variables must be finite, non-missing and numeric.",
         call. = FALSE)
  }
  if (!is.logical(scale) || length(scale) != 1L || is.na(scale)) {
    stop("`scale` must be TRUE or FALSE.", call. = FALSE)
  }
  integer_scalar <- function(value, name, minimum) {
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        !is.finite(value) || value < minimum || value != as.integer(value)) {
      stop(sprintf("`%s` must be one integer not less than %d.", name, minimum),
           call. = FALSE)
    }
    as.integer(value)
  }
  folds <- integer_scalar(folds, "folds", 2L)
  nstart <- integer_scalar(nstart, "nstart", 1L)
  key <- paste(format(object$coordinates[, 1L], digits = 17, scientific = TRUE),
               format(object$coordinates[, 2L], digits = 17, scientific = TRUE),
               sep = "\r")
  groups <- split(seq_along(key), key)
  inconsistent <- vapply(groups, function(index) {
    any(vapply(environmental, function(value)
      length(unique(value[index])) > 1L, logical(1)))
  }, logical(1))
  if (any(inconsistent)) {
    stop("Repeated coordinates must have identical environmental values.",
         call. = FALSE)
  }
  first <- !duplicated(key)
  unique_environmental <- as.matrix(environmental[first, , drop = FALSE])
  unique_key <- key[first]
  if (folds > nrow(unique_environmental)) {
    stop("`folds` cannot exceed the number of unique coordinate locations.",
         call. = FALSE)
  }
  centre <- rep(0, ncol(unique_environmental))
  spread <- rep(1, ncol(unique_environmental))
  clustering_values <- unique_environmental
  if (scale) {
    centre <- colMeans(unique_environmental)
    spread <- apply(unique_environmental, 2L, stats::sd)
    if (any(!is.finite(spread)) || any(spread <= 0)) {
      stop("Scaled environmental variables must vary across unique locations.",
           call. = FALSE)
    }
    clustering_values <- sweep(
      sweep(unique_environmental, 2L, centre, "-"), 2L, spread, "/"
    )
  }
  seed_scope <- .sp_seed_scope(seed)
  if (!is.null(seed_scope)) {
    on.exit(.sp_restore_seed(seed_scope), add = TRUE)
    set.seed(seed_scope$seed)
  }
  clustering <- stats::kmeans(clustering_values, centers = folds,
                              nstart = nstart)
  canonical <- do.call(order, as.data.frame(clustering$centers))
  canonical_label <- integer(folds)
  canonical_label[canonical] <- seq_len(folds)
  unique_fold <- canonical_label[clustering$cluster]
  resolved <- unique_fold[match(key, unique_key)]
  centres <- clustering$centers[canonical, , drop = FALSE]
  colnames(centres) <- variables
  row.names(centres) <- as.character(seq_len(folds))
  structure(as.integer(resolved),
    class = c("spmodels_environmental_folds", "integer"),
    centres = centres, variables = variables, scaled = scale,
    scaling_centre = stats::setNames(centre, variables),
    scaling_scale = stats::setNames(spread, variables), seed = seed,
    nstart = nstart, call = match.call())
}

#' @export
print.spmodels_environmental_folds <- function(x, ...) {
  cat("SPmodels environmental-space folds\n")
  cat("Variables:", paste(attr(x, "variables"), collapse = ", "),
      " Scaled:", attr(x, "scaled"), "\n")
  print(table(fold = as.integer(x)))
  invisible(x)
}

#' Construct buffered spatial validation folds
#'
#' Start from geometry-driven test folds and declare an exclusion distance
#' around every held-out location. [sp_cross_validate()] removes training rows
#' inside that buffer separately for each fold.
#'
#' @param object A fitted spatial `spmodels_fit` object.
#' @param folds Number of geometry-driven folds.
#' @param buffer Finite positive exclusion distance in the fitted coordinate
#'   units.
#' @param seed Optional local seed for fold construction.
#' @param nstart Positive number of k-means starts.
#'
#' @return An integer vector of fold identifiers with class
#'   `spmodels_buffered_folds` and its buffer metadata.
#' @export
sp_buffered_folds <- function(object, folds = 5L, buffer, seed = NULL,
                              nstart = 20L) {
  if (!is.numeric(buffer) || length(buffer) != 1L || is.na(buffer) ||
      !is.finite(buffer) || buffer <= 0) {
    stop("`buffer` must be one finite positive distance.", call. = FALSE)
  }
  answer <- sp_spatial_folds(object, folds = folds, seed = seed,
                             nstart = nstart)
  attr(answer, "buffer") <- as.numeric(buffer)
  attr(answer, "distance_units") <- object$distance_units
  class(answer) <- c("spmodels_buffered_folds", class(answer))
  answer
}

#' @export
print.spmodels_buffered_folds <- function(x, ...) {
  cat("SPmodels buffered spatial folds\n")
  print(table(fold = as.integer(x)))
  cat("Training exclusion distance:", format(attr(x, "buffer")),
      attr(x, "distance_units"), "\n")
  invisible(x)
}

#' Refit fixed-fold spatial cross-validation
#'
#' Refit an SPmodels model after withholding each declared fold, then predict
#' the held-out observations without using them during fitting.
#'
#' @param object A fitted `spmodels_fit` object retaining its aligned source
#'   data and controls.
#' @param folds Either an integer number of balanced randomly assigned folds,
#'   one fixed fold identifier per fitted observation, or an object returned by
#'   [sp_buffered_folds()].
#' @param seed Optional local seed used only when `folds` is a count.
#' @param unseen Policy for held-out grouping levels absent from a training
#'   fold. The safe default errors; `"marginal"` integrates new zero-mean
#'   grouped effects.
#'
#' @return An object of class `spmodels_cross_validation` with row predictions,
#'   overall metrics, fold metrics and the resolved fold vector.
#' @export
sp_cross_validate <- function(object, folds = 5L, seed = NULL,
                              unseen = c("error", "marginal")) {
  if (!inherits(object, "spmodels_fit")) {
    stop("`object` must be a fitted `spmodels_fit` model.", call. = FALSE)
  }
  if (isTRUE(object$external_offset)) {
    stop("Cross-validation currently requires offsets to be written in the formula.",
         call. = FALSE)
  }
  unseen <- match.arg(unseen)
  n <- object$nobs
  fold_id <- .sp_resolve_folds(folds, n, seed)
  buffered <- inherits(folds, "spmodels_buffered_folds")
  buffer <- if (buffered) attr(folds, "buffer") else NULL
  prediction <- standard_error <- rep(NA_real_, n)
  integrated_log_score <- rep(NA_real_, n)
  training_size <- stats::setNames(integer(length(unique(fold_id))),
                                   sort(unique(fold_id)))
  for (fold in sort(unique(fold_id))) {
    test <- fold_id == fold
    train <- !test
    if (buffered) {
      training_index <- which(train)
      test_index <- which(test)
      distance_to_test <- .sp_dist(object$coordinates[training_index, , drop = FALSE],
                                   object$coordinates[test_index, , drop = FALSE])
      excluded <- apply(distance_to_test <= buffer, 1L, any)
      train[training_index[excluded]] <- FALSE
    }
    training_size[[as.character(fold)]] <- sum(train)
    if (!any(train)) {
      stop(sprintf("Fold %s has no training rows after buffer exclusion.",
                   as.character(fold)), call. = FALSE)
    }
    train_data <- object$data_used[train, , drop = FALSE]
    test_data <- object$data_used[test, , drop = FALSE]
    train_weights <- if (object$family$family == "binomial" &&
      identical(object$response_encoding, "proportion"))
      object$trials[train] else NULL
    fold_fit <- sp_fit(
      object$formula, train_data, family = object$family,
      coords = if (inherits(train_data, "sf")) NULL else object$coords,
      spatial = object$spatial, time = object$temporal, random = object$random,
      random_slope = object$random_slope,
      intervention = object$intervention,
      weights = train_weights,
      control = object$control
    )
    if (identical(object$likelihood_method, "Monte Carlo ML")) {
      fold_fit <- sp_mcml(
        fold_fit, control = object$mcml$control, seed = object$mcml$seed
      )
    }
    fold_prediction <- stats::predict(
      fold_fit, test_data, se.fit = TRUE, target = "observation",
      type = "response", unseen = unseen
    )
    prediction[test] <- fold_prediction$fit
    standard_error[test] <- fold_prediction$se.fit
    if (object$family$family != "gaussian") {
      integrated_log_score[test] <- sp_predictive_log_score(
        fold_fit, test_data, object$y[test],
        trials = if (object$family$family == "binomial")
          object$trials[test] else NULL,
        unseen = unseen
      )
    }
  }
  observed <- if (object$family$family == "binomial") {
    ifelse(object$trials > 0, object$y / object$trials, NA_real_)
  } else object$y
  error <- observed - prediction
  log_score <- if (object$family$family == "gaussian") {
    -stats::dnorm(object$y, prediction, standard_error, log = TRUE)
  } else integrated_log_score
  rows <- data.frame(
    row = row.names(object$model), fold = fold_id, observed = observed,
    predicted = prediction, se.fit = standard_error, error = error,
    squared_error = error^2, absolute_error = abs(error),
    log_score = log_score, stringsAsFactors = FALSE
  )
  summarise <- function(index) c(
    n = sum(index), rmse = sqrt(mean(rows$squared_error[index], na.rm = TRUE)),
    mae = mean(rows$absolute_error[index], na.rm = TRUE),
    mean_log_score = mean(rows$log_score[index], na.rm = TRUE)
  )
  summary <- as.data.frame(as.list(summarise(rep(TRUE, n))))
  fold_summary <- do.call(rbind, lapply(sort(unique(fold_id)), function(fold) {
    data.frame(fold = fold, as.list(summarise(fold_id == fold)),
               row.names = NULL)
  }))
  structure(list(
    predictions = rows, summary = summary, fold_summary = fold_summary,
    folds = fold_id, family = object$family$family,
    training_size = training_size, buffer = buffer, unseen = unseen,
    call = match.call()
  ), class = "spmodels_cross_validation")
}

.sp_empirical_crps <- function(draws, truth) {
  b <- length(draws)
  ordered <- sort(draws)
  coefficient <- 2 * seq_len(b) - b - 1
  mean(abs(draws - truth)) - sum(coefficient * ordered) / b^2
}

#' Assess predictive draws against known or observed values
#'
#' Calculate location-level errors, empirical CRPS and interval coverage from
#' joint predictive draws, plus overall recovery summaries.
#'
#' @param object An object returned by [sp_predict_draws()].
#' @param truth One finite value for every prediction location.
#' @param interval Central interval probability strictly between zero and one.
#' @param threshold Optional finite threshold for predicted exceedance
#'   probabilities and observed threshold indicators.
#'
#' @return An object of class `spmodels_simulation_assessment` containing
#'   location-level and overall summaries.
#' @export
sp_simulation_assess <- function(object, truth, interval = 0.95,
                                 threshold = NULL) {
  if (!inherits(object, "spmodels_prediction_draws")) {
    stop("`object` must be returned by `sp_predict_draws()`.", call. = FALSE)
  }
  if (!is.numeric(truth) || length(truth) != ncol(object$draws) ||
      anyNA(truth) || any(!is.finite(truth))) {
    stop("`truth` must contain one finite value per prediction location.",
         call. = FALSE)
  }
  if (!is.numeric(interval) || length(interval) != 1L || is.na(interval) ||
      !is.finite(interval) || interval <= 0 || interval >= 1) {
    stop("`interval` must be one number strictly between zero and one.",
         call. = FALSE)
  }
  if (!is.null(threshold) && (!is.numeric(threshold) ||
      length(threshold) != 1L || is.na(threshold) || !is.finite(threshold))) {
    stop("`threshold` must be NULL or one finite number.", call. = FALSE)
  }
  alpha <- (1 - interval) / 2
  means <- colMeans(object$draws)
  limits <- t(apply(object$draws, 2L, stats::quantile,
                    probs = c(alpha, 1 - alpha), names = FALSE, type = 8))
  error <- means - truth
  location <- data.frame(
    location = if (is.null(colnames(object$draws)))
      paste0("location_", seq_along(truth)) else colnames(object$draws),
    truth = truth, mean = means, error = error,
    squared_error = error^2, absolute_error = abs(error),
    lower = limits[, 1L], upper = limits[, 2L],
    covered = truth >= limits[, 1L] & truth <= limits[, 2L],
    interval_width = limits[, 2L] - limits[, 1L],
    crps = vapply(seq_along(truth), function(j)
      .sp_empirical_crps(object$draws[, j], truth[j]), numeric(1)),
    stringsAsFactors = FALSE
  )
  if (!is.null(threshold)) {
    location$exceedance_probability <- colMeans(object$draws > threshold)
    location$truth_exceeds <- truth > threshold
  }
  summary <- data.frame(
    locations = length(truth), bias = mean(error),
    rmse = sqrt(mean(error^2)), mae = mean(abs(error)),
    coverage = mean(location$covered),
    mean_interval_width = mean(location$interval_width),
    mean_crps = mean(location$crps)
  )
  structure(list(
    locations = location, summary = summary, interval = interval,
    threshold = threshold, family = object$family, target = object$target,
    type = object$type, call = match.call()
  ), class = "spmodels_simulation_assessment")
}

#' @export
print.spmodels_cross_validation <- function(x, ...) {
  cat("SPmodels fixed-fold cross-validation\n")
  print(x$summary, row.names = FALSE)
  invisible(x)
}

#' @export
print.spmodels_simulation_assessment <- function(x, ...) {
  cat("SPmodels predictive-draw assessment\n")
  print(x$summary, row.names = FALSE)
  invisible(x)
}

.sp_simulation_data <- function(object, simulated_response) {
  data <- object$data_used
  response <- object$formula[[2L]]
  if (is.symbol(response)) {
    response_name <- as.character(response)
    if (object$family$family == "binomial" &&
        identical(object$response_encoding, "proportion")) {
      data[[response_name]] <- ifelse(object$trials > 0,
                                      simulated_response / object$trials, 0)
    } else {
      data[[response_name]] <- simulated_response
    }
    return(data)
  }
  if (is.call(response) && identical(response[[1L]], quote(cbind)) &&
      length(response) == 3L && is.symbol(response[[2L]]) &&
      object$family$family == "binomial") {
    data[[as.character(response[[2L]])]] <- simulated_response
    if (is.symbol(response[[3L]])) {
      data[[as.character(response[[3L]])]] <-
        object$trials - simulated_response
    }
    return(data)
  }
  stop("Simulation refitting supports a named response or a two-column binomial response.",
       call. = FALSE)
}

#' Repeated fitted-model simulation and parameter recovery
#'
#' Simulate response datasets from a fitted model, refit the same declared
#' model to every replicate, and summarize parameter recovery.
#'
#' @param object A fitted `spmodels_fit` object.
#' @param nsim Positive number of simulation/refit replicates.
#' @param seed Optional local seed for response simulation and any MCML refits.
#' @param level Wald interval probability strictly between zero and one.
#' @param on_error Either record failed refits and continue or stop on failure.
#'   Sequential execution stops immediately; already-dispatched parallel work
#'   may finish before the error is returned.
#' @param workers Positive number of replicate workers. Values above one use
#'   forked processes on non-Windows systems. Stochastic MCML refits then
#'   require an explicit `seed`.
#'
#' @return An object of class `spmodels_simulation_recovery` containing the
#'   simulated responses, replicate estimates, standard errors, convergence,
#'   failures and parameter-level recovery summary.
#' @export
sp_simulation_recovery <- function(object, nsim = 100L, seed = NULL,
                                   level = 0.95,
                                   on_error = c("record", "stop"),
                                   workers = 1L) {
  if (!inherits(object, "spmodels_fit")) {
    stop("`object` must be a fitted `spmodels_fit` model.", call. = FALSE)
  }
  if (isTRUE(object$external_offset)) {
    stop("Simulation refitting currently requires offsets to be written in the formula.",
         call. = FALSE)
  }
  if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) ||
      !is.finite(nsim) || nsim < 1 || nsim != as.integer(nsim)) {
    stop("`nsim` must be one positive integer.", call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      !is.finite(level) || level <= 0 || level >= 1) {
    stop("`level` must be one number strictly between zero and one.",
         call. = FALSE)
  }
  on_error <- match.arg(on_error)
  if (!is.numeric(workers) || length(workers) != 1L || is.na(workers) ||
      !is.finite(workers) || workers < 1L || workers != as.integer(workers)) {
    stop("`workers` must be one positive integer.", call. = FALSE)
  }
  workers <- as.integer(workers)
  if (workers > 1L && !.sp_has_fork_workers()) {
    stop(paste0("Fork-based parallel recovery replicates are unavailable ",
                "on Windows and WebAssembly; use `workers = 1`."),
         call. = FALSE)
  }
  if (workers > 1L && identical(object$likelihood_method, "Monte Carlo ML") &&
      is.null(seed)) {
    stop("Parallel MCML recovery requires an explicit `seed`.", call. = FALSE)
  }
  simulations <- stats::simulate(object, nsim = as.integer(nsim), seed = seed)
  slope_truth <- .sp_random_slope_estimated_parameters(
    object$random_slope_specification
  )
  covariance_names <- .sp_fitted_covariance_names(object)
  truth <- c(object$coefficients,
    if (length(covariance_names)) object$covariance[covariance_names] else NULL,
    if (object$covariance_estimated && length(object$random_variance))
      stats::setNames(object$random_variance,
                      paste0("random::", names(object$random_variance))) else NULL,
    if (length(slope_truth)) stats::setNames(
      slope_truth, paste0("random_slope::", names(slope_truth))
    ) else NULL,
    if (isTRUE(object$temporal_scale_estimated))
      c(temporal_scale = object$temporal_scale) else NULL,
    if (isTRUE(object$intervention_estimated))
      c(alpha = object$intervention$alpha,
        decay = object$intervention$decay) else NULL)
  parameter_names <- names(truth)
  estimates <- standard_errors <- matrix(
    NA_real_, nrow = as.integer(nsim), ncol = length(truth),
    dimnames = list(paste0("sim_", seq_len(as.integer(nsim))),
                    parameter_names)
  )
  convergence <- rep(NA_integer_, as.integer(nsim))
  failures <- rep(NA_character_, as.integer(nsim))
  run_replicate <- function(simulation_index) {
    attempt <- tryCatch({
      simulated_data <- .sp_simulation_data(
        object, simulations[[simulation_index]]
      )
      simulation_weights <- if (object$family$family == "binomial" &&
        identical(object$response_encoding, "proportion"))
        object$trials else NULL
      refit <- sp_fit(
        object$formula, simulated_data, family = object$family,
        coords = if (inherits(simulated_data, "sf")) NULL else object$coords,
        spatial = object$spatial, time = object$temporal,
        random = object$random, random_slope = object$random_slope,
        intervention = object$intervention,
        weights = simulation_weights,
        control = object$control
      )
      if (identical(object$likelihood_method, "Monte Carlo ML")) {
        refit_seed <- if (is.null(seed)) NULL else
          as.integer((as.double(seed) + simulation_index) %%
                       .Machine$integer.max)
        refit <- sp_mcml(refit, control = object$mcml$control,
                         seed = refit_seed)
      }
      refit
    }, error = identity)
    if (inherits(attempt, "error")) {
      if (on_error == "stop" && workers == 1L) stop(attempt)
      return(list(error = conditionMessage(attempt)))
    }
    fitted_parameters <- c(attempt$coefficients,
      if (length(covariance_names))
        attempt$covariance[covariance_names] else NULL,
      if (object$covariance_estimated && length(object$random_variance))
        stats::setNames(attempt$random_variance,
                        paste0("random::", names(attempt$random_variance))) else NULL,
      if (length(slope_truth)) stats::setNames(
        .sp_random_slope_estimated_parameters(
          attempt$random_slope_specification
        ), paste0("random_slope::", names(slope_truth))
      ) else NULL,
      if (isTRUE(object$temporal_scale_estimated))
        c(temporal_scale = attempt$temporal_scale) else NULL,
      if (isTRUE(object$intervention_estimated))
        c(alpha = attempt$intervention$alpha,
          decay = attempt$intervention$decay) else NULL)
    fitted_se <- c(sqrt(diag(attempt$coefficient_vcov)),
      if (length(covariance_names))
        sqrt(diag(attempt$covariance_vcov))[covariance_names] else NULL,
      if (object$covariance_estimated && length(object$random_variance))
        stats::setNames(sqrt(diag(attempt$random_variance_vcov)),
          paste0("random::", names(attempt$random_variance))) else NULL,
      if (length(slope_truth)) stats::setNames(
        sqrt(diag(attempt$random_slope_vcov))[names(slope_truth)],
        paste0("random_slope::", names(slope_truth))
      ) else NULL,
      if (isTRUE(object$temporal_scale_estimated))
        c(temporal_scale = attempt$temporal_scale_se) else NULL,
      if (isTRUE(object$intervention_estimated))
        sqrt(diag(attempt$intervention_vcov)) else NULL)
    list(error = NA_character_,
         estimates = fitted_parameters[parameter_names],
         standard_errors = fitted_se[parameter_names],
         convergence = if (is.null(attempt$optimisation)) 0L else
           attempt$optimisation$convergence)
  }
  replicate_results <- if (workers == 1L) {
    lapply(seq_len(as.integer(nsim)), run_replicate)
  } else {
    parallel::mclapply(seq_len(as.integer(nsim)), run_replicate,
      mc.cores = min(workers, as.integer(nsim)), mc.preschedule = FALSE,
      mc.set.seed = FALSE)
  }
  for (simulation_index in seq_len(as.integer(nsim))) {
    result <- replicate_results[[simulation_index]]
    if (inherits(result, "try-error")) {
      result <- list(error = as.character(result))
    }
    if (!is.na(result$error)) {
      failures[[simulation_index]] <- result$error
      if (on_error == "stop") stop(result$error, call. = FALSE)
    } else {
      estimates[simulation_index, ] <- result$estimates
      standard_errors[simulation_index, ] <- result$standard_errors
      convergence[[simulation_index]] <- result$convergence
    }
  }
  z <- stats::qnorm((1 + level) / 2)
  parameter_summary <- do.call(rbind, lapply(seq_along(truth), function(j) {
    valid <- is.finite(estimates[, j])
    valid_se <- valid & is.finite(standard_errors[, j])
    errors <- estimates[valid, j] - truth[j]
    data.frame(
      parameter = parameter_names[j], truth = unname(truth[j]),
      successful = sum(valid), mean = if (any(valid))
        mean(estimates[valid, j]) else NA_real_,
      bias = if (any(valid)) mean(errors) else NA_real_,
      rmse = if (any(valid)) sqrt(mean(errors^2)) else NA_real_,
      empirical_sd = if (sum(valid) > 1L)
        stats::sd(estimates[valid, j]) else NA_real_,
      mean_reported_se = if (any(valid_se))
        mean(standard_errors[valid_se, j]) else NA_real_,
      coverage = if (any(valid_se)) mean(
        truth[j] >= estimates[valid_se, j] - z * standard_errors[valid_se, j] &
        truth[j] <= estimates[valid_se, j] + z * standard_errors[valid_se, j]
      ) else NA_real_, row.names = NULL
    )
  }))
  structure(list(
    simulations = simulations, estimates = estimates,
    standard_errors = standard_errors, truth = truth,
    summary = parameter_summary, convergence = convergence,
    failures = failures, level = level, seed = seed, workers = workers,
    call = match.call()
  ), class = "spmodels_simulation_recovery")
}

#' @export
print.spmodels_simulation_recovery <- function(x, ...) {
  cat("SPmodels repeated simulation and parameter recovery\n")
  cat("Replicate workers:", x$workers, "\n")
  cat("Successful refits:", sum(is.na(x$failures)), "of",
      length(x$failures), "\n")
  print(x$summary, row.names = FALSE)
  invisible(x)
}

#' Construct a Cartesian grid of simulation-recovery scenarios
#'
#' Create explicitly named generative scenario objects by replacing recoverable
#' fitted-model parameters over a Cartesian product. The returned objects are
#' intended for [sp_simulation_scenarios()], not as refitted analyses.
#'
#' @param object A fitted `spmodels_fit` object that defines the design,
#'   response family and refitting specification.
#' @param values Named list of non-empty numeric vectors. Names must identify
#'   recoverable mean coefficients, covariance parameters, random-intercept
#'   variances written as `random::group`, estimated grouped covariance entries
#'   using their `random_slope::` model-table names, `temporal_scale`, or
#'   intervention parameters `alpha` and `decay` when estimated.
#' @param prefix Non-empty prefix used for deterministic scenario names.
#'
#' @return A named list of generative `spmodels_fit` scenario objects with a
#'   `grid` attribute recording the parameter combinations.
#' @export
sp_simulation_grid <- function(object, values, prefix = "scenario") {
  if (!inherits(object, "spmodels_fit")) {
    stop("`object` must be a fitted `spmodels_fit` model.", call. = FALSE)
  }
  if (!is.list(values) || !length(values) || is.null(names(values)) ||
      anyNA(names(values)) || any(!nzchar(names(values))) ||
      anyDuplicated(names(values))) {
    stop("`values` must be a non-empty named list with unique parameter names.",
         call. = FALSE)
  }
  if (!is.character(prefix) || length(prefix) != 1L || is.na(prefix) ||
      !nzchar(trimws(prefix))) {
    stop("`prefix` must be one non-empty character value.", call. = FALSE)
  }
  valid_values <- vapply(values, function(value) {
    is.numeric(value) && length(value) && !anyNA(value) &&
      all(is.finite(value))
  }, logical(1))
  if (!all(valid_values)) {
    stop("Every element of `values` must be a non-empty finite numeric vector.",
         call. = FALSE)
  }
  slope_parameters <- .sp_random_slope_estimated_parameters(
    object$random_slope_specification
  )
  recoverable <- c(names(object$coefficients),
    .sp_fitted_covariance_names(object),
    if (isTRUE(object$covariance_estimated) && length(object$random_variance))
      paste0("random::", names(object$random_variance)) else NULL,
    if (length(slope_parameters))
      paste0("random_slope::", names(slope_parameters)) else NULL,
    if (isTRUE(object$temporal_scale_estimated)) "temporal_scale" else NULL,
    if (isTRUE(object$intervention_estimated)) c("alpha", "decay") else NULL)
  duplicated_parameters <- unique(recoverable[
    duplicated(recoverable) | duplicated(recoverable, fromLast = TRUE)
  ])
  ambiguous <- intersect(names(values), duplicated_parameters)
  if (length(ambiguous)) {
    stop(sprintf("Grid parameter names are ambiguous in this model: %s.",
                 paste(ambiguous, collapse = ", ")), call. = FALSE)
  }
  unknown <- setdiff(names(values), recoverable)
  if (length(unknown)) {
    stop(sprintf("Grid parameters are not recoverable in this model: %s.",
                 paste(unknown, collapse = ", ")), call. = FALSE)
  }
  positive <- intersect(names(values), c(
    names(object$covariance), paste0("random::", names(object$random_variance)),
    paste0("random_slope::", names(slope_parameters)[
      grepl("::var(", names(slope_parameters), fixed = TRUE)
    ]), "temporal_scale", "decay"
  ))
  if (length(positive) && any(vapply(values[positive], function(x) any(x <= 0),
                                     logical(1)))) {
    stop("Variance, scale and decay grid values must be strictly positive.",
         call. = FALSE)
  }
  if ("alpha" %in% names(values) &&
      any(values$alpha <= 0 | values$alpha >= 1)) {
    stop("Estimated `alpha` grid values must lie strictly between zero and one.",
         call. = FALSE)
  }
  grid <- do.call(expand.grid, c(values, list(
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )))
  names(grid) <- names(values)
  rebuild <- function(row) {
    scenario <- object
    slope_updates <- numeric()
    for (name in names(row)) {
      value <- unname(row[[name]])
      if (name %in% names(scenario$coefficients)) {
        scenario$coefficients[[name]] <- value
      } else if (name %in% names(scenario$covariance)) {
        scenario$covariance[[name]] <- value
        scenario$control[[name]] <- value
      } else if (startsWith(name, "random::")) {
        group <- sub("^random::", "", name)
        scenario$random_variance[[group]] <- value
        scenario$control$random_variance <- scenario$random_variance
      } else if (startsWith(name, "random_slope::")) {
        slope_name <- sub("^random_slope::", "", name)
        slope_updates[[slope_name]] <- value
      } else if (name == "temporal_scale") {
        scenario$temporal_scale <- value
        scenario$temporal$scale <- value
      } else if (name %in% c("alpha", "decay")) {
        scenario$intervention[[name]] <- value
      }
    }
    if (length(slope_updates)) {
      scenario$random_slope_specification <- .sp_random_slope_set_parameters(
        scenario$random_slope_specification, slope_updates
      )
      scenario$random_slope <- .sp_random_slope_export(
        scenario$random_slope_specification
      )
      scenario$random_slope_parameterization <-
        .sp_random_slope_parameterization(
          scenario$random_slope_specification
        )
    }
    temporal_lag <- .sp_temporal_lag(scenario$temporal_values)
    random_z <- if (is.null(scenario$random_specification)) NULL else
      scenario$random_specification$z
    random_group <- if (is.null(scenario$random_specification)) NULL else
      scenario$random_specification$column_group
    distance <- .sp_dist(scenario$coordinates)
    if (scenario$family$family == "gaussian") {
      state <- .sp_gaussian_state(
        scenario$covariance, distance, scenario$x,
        scenario$y - scenario$offset, scenario$smoothness,
        random_z = random_z, random_variance = scenario$random_variance,
        random_column_group = random_group,
        random_slope_specification = scenario$random_slope_specification,
        temporal_lag = temporal_lag,
        temporal_scale = scenario$temporal_scale
      )
      if (is.null(state)) stop("A grid scenario covariance is not positive definite.",
                               call. = FALSE)
      state$beta <- scenario$coefficients
      scenario$fitted_state <- state
    } else {
      multiplier <- .sp_intervention_multiplier(
        scenario$intervention, scenario$data_used
      )
      state <- .sp_laplace_state(
        scenario$coefficients, scenario$covariance, distance, scenario$x,
        scenario$offset, scenario$y, scenario$trials, scenario$family,
        scenario$smoothness, multiplier, temporal_lag,
        scenario$temporal_scale, random_z, scenario$random_variance,
        random_group, scenario$random_slope_specification
      )
      if (is.null(state)) {
        stop("A grid scenario could not form a valid latent Gaussian state.",
             call. = FALSE)
      }
      scenario$intervention_multiplier <- multiplier
      scenario$fitted_state <- state
      scenario$latent_mode <- state$mode
    }
    attr(scenario, "simulation_truth") <- unlist(row, use.names = TRUE)
    scenario
  }
  scenarios <- lapply(seq_len(nrow(grid)), function(i) {
    rebuild(as.list(grid[i, , drop = FALSE]))
  })
  width <- nchar(nrow(grid))
  names(scenarios) <- sprintf(paste0(trimws(prefix), "_%0", width, "d"),
                              seq_len(nrow(grid)))
  structure(scenarios, grid = grid,
            class = c("spmodels_simulation_grid", "list"))
}

#' Run parameter-recovery experiments across fitted scenarios
#'
#' Apply [sp_simulation_recovery()] to each named fitted model with explicit,
#' reproducible scenario seeds, retaining every individual result and a
#' combined parameter summary.
#'
#' @param objects Named non-empty list of fitted `spmodels_fit` objects.
#' @param nsim One positive replicate count for every scenario, or one count
#'   per scenario. A named vector is reordered to the scenario names.
#' @param seed Optional base integer seed. Scenario `i` uses `seed + i - 1`
#'   modulo R's maximum integer.
#' @param level Wald interval probability strictly between zero and one.
#' @param on_error Passed to [sp_simulation_recovery()].
#' @param workers Positive number of scenario workers. Values above one use
#'   forked processes on non-Windows systems; Windows currently requires one.
#'
#' @return An object of class `spmodels_simulation_scenarios` containing the
#'   named recovery results, combined summary and scenario metadata.
#' @export
sp_simulation_scenarios <- function(objects, nsim = 100L, seed = NULL,
                                    level = 0.95,
                                    on_error = c("record", "stop"),
                                    workers = 1L) {
  if (!is.list(objects) || !length(objects) ||
      any(!vapply(objects, inherits, logical(1), what = "spmodels_fit"))) {
    stop("`objects` must be a non-empty list of fitted `spmodels_fit` models.",
         call. = FALSE)
  }
  scenario <- names(objects)
  if (is.null(scenario) || anyNA(scenario) || any(!nzchar(scenario)) ||
      anyDuplicated(scenario)) {
    stop("`objects` must have unique, non-empty scenario names.", call. = FALSE)
  }
  valid_count <- function(value) is.numeric(value) && !is.na(value) &&
    is.finite(value) && value >= 1 && value == as.integer(value)
  if (length(nsim) == 1L) {
    nsim <- rep(nsim, length(objects))
  } else if (!is.null(names(nsim))) {
    if (!setequal(names(nsim), scenario) || anyDuplicated(names(nsim))) {
      stop("Named `nsim` values must match the scenario names exactly.",
           call. = FALSE)
    }
    nsim <- nsim[scenario]
  }
  if (length(nsim) != length(objects) ||
      any(!vapply(as.list(nsim), valid_count, logical(1)))) {
    stop("`nsim` must supply positive integer counts for all scenarios.",
         call. = FALSE)
  }
  nsim <- as.integer(nsim)
  names(nsim) <- scenario
  if (!is.numeric(workers) || length(workers) != 1L || is.na(workers) ||
      !is.finite(workers) || workers < 1L || workers != as.integer(workers)) {
    stop("`workers` must be one positive integer.", call. = FALSE)
  }
  workers <- as.integer(workers)
  if (workers > 1L && !.sp_has_fork_workers()) {
    stop(paste0("Fork-based parallel scenarios are unavailable on Windows ",
                "and WebAssembly; use `workers = 1`."),
         call. = FALSE)
  }
  seed_scope <- .sp_seed_scope(seed)
  scenario_seed <- if (is.null(seed_scope)) {
    rep(NA_integer_, length(objects))
  } else {
    as.integer((as.double(seed_scope$seed) + seq_along(objects) - 1) %%
      .Machine$integer.max)
  }
  names(scenario_seed) <- scenario
  on_error <- match.arg(on_error)
  run_scenario <- function(i) {
    sp_simulation_recovery(
      objects[[i]], nsim = nsim[[i]],
      seed = if (is.na(scenario_seed[[i]])) NULL else scenario_seed[[i]],
      level = level, on_error = on_error
    )
  }
  results <- if (workers == 1L) {
    lapply(seq_along(objects), run_scenario)
  } else {
    parallel::mclapply(seq_along(objects), run_scenario,
      mc.cores = min(workers, length(objects)), mc.preschedule = FALSE,
      mc.set.seed = FALSE)
  }
  names(results) <- scenario
  summary <- do.call(rbind, lapply(seq_along(results), function(i) {
    data.frame(scenario = scenario[[i]], results[[i]]$summary,
               row.names = NULL, check.names = FALSE)
  }))
  metadata <- data.frame(
    scenario = scenario, family = vapply(objects, function(x) x$family$family,
                                          character(1)),
    replicates = unname(nsim), seed = unname(scenario_seed),
    successful_refits = vapply(results, function(x) sum(is.na(x$failures)),
                                integer(1)),
    stringsAsFactors = FALSE
  )
  structure(list(
    results = results, summary = summary, scenarios = metadata,
    level = level, base_seed = seed, workers = workers, call = match.call()
  ), class = "spmodels_simulation_scenarios")
}

#' @export
print.spmodels_simulation_scenarios <- function(x, ...) {
  cat("SPmodels multi-scenario simulation recovery\n")
  cat("Scenario workers:", x$workers, "\n")
  print(x$scenarios, row.names = FALSE)
  cat("\nParameter recovery:\n")
  print(x$summary, row.names = FALSE)
  invisible(x)
}

#' Precision-directed simulation recovery
#'
#' Increase the number of fitted-model simulation/refit replicates until every
#' recovered parameter's Monte Carlo standard error for mean bias meets a
#' declared target, or an exact replicate cap is evaluated.
#'
#' @param object A fitted `spmodels_fit` object.
#' @param initial_nsim Initial positive replicate count.
#' @param max_nsim Maximum replicate count, not smaller than `initial_nsim`.
#' @param growth Multiplicative increase strictly greater than one.
#' @param target_mcse Positive target for every parameter's bias MCSE.
#' @param min_success Minimum successful refits required per parameter.
#' @param seed Optional base seed; attempt seeds advance deterministically.
#' @param level,on_error Passed to [sp_simulation_recovery()].
#' @param workers Positive number of replicate workers passed to every complete
#'   recovery attempt.
#'
#' @return An object of class `spmodels_simulation_precision` containing the
#'   final recovery result, attempt history, achieved MCSE values and status.
#' @export
sp_simulation_precision <- function(
    object, initial_nsim = 20L, max_nsim = 1000L, growth = 2,
    target_mcse, min_success = 10L, seed = NULL, level = 0.95,
    on_error = c("record", "stop"), workers = 1L) {
  positive_integer <- function(value, name) {
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        !is.finite(value) || value < 1L || value != as.integer(value)) {
      stop(sprintf("`%s` must be one positive integer.", name), call. = FALSE)
    }
    as.integer(value)
  }
  initial_nsim <- positive_integer(initial_nsim, "initial_nsim")
  max_nsim <- positive_integer(max_nsim, "max_nsim")
  min_success <- positive_integer(min_success, "min_success")
  if (max_nsim < initial_nsim) {
    stop("`max_nsim` must not be smaller than `initial_nsim`.", call. = FALSE)
  }
  if (!is.numeric(growth) || length(growth) != 1L || is.na(growth) ||
      !is.finite(growth) || growth <= 1) {
    stop("`growth` must be one finite number greater than one.", call. = FALSE)
  }
  if (!is.numeric(target_mcse) || length(target_mcse) != 1L ||
      is.na(target_mcse) || !is.finite(target_mcse) || target_mcse <= 0) {
    stop("`target_mcse` must be one finite positive number.", call. = FALSE)
  }
  seed_scope <- .sp_seed_scope(seed)
  on_error <- match.arg(on_error)
  workers <- positive_integer(workers, "workers")
  if (workers > 1L && !.sp_has_fork_workers()) {
    stop(paste0("Fork-based parallel recovery replicates are unavailable ",
                "on Windows and WebAssembly; use `workers = 1`."),
         call. = FALSE)
  }
  if (workers > 1L && identical(object$likelihood_method, "Monte Carlo ML") &&
      is.null(seed)) {
    stop("Parallel MCML recovery requires an explicit `seed`.", call. = FALSE)
  }
  current <- initial_nsim
  attempt <- 0L
  history <- list()
  repeat {
    attempt <- attempt + 1L
    attempt_seed <- if (is.null(seed_scope)) NULL else
      as.integer((as.double(seed_scope$seed) + attempt - 1) %%
                   .Machine$integer.max)
    result <- sp_simulation_recovery(
      object, nsim = current, seed = attempt_seed, level = level,
      on_error = on_error, workers = workers
    )
    successful <- colSums(is.finite(result$estimates))
    empirical_sd <- apply(result$estimates, 2L, function(value) {
      value <- value[is.finite(value)]
      if (length(value) > 1L) stats::sd(value) else NA_real_
    })
    mcse <- empirical_sd / sqrt(successful)
    parameter_pass <- successful >= min_success & is.finite(mcse) &
      mcse <= target_mcse
    passed <- length(parameter_pass) > 0L && all(parameter_pass)
    history[[attempt]] <- data.frame(
      attempt = attempt, replicates = current,
      seed = if (is.null(attempt_seed)) NA_integer_ else attempt_seed,
      minimum_successful = min(successful),
      maximum_mcse = if (any(is.finite(mcse))) max(mcse, na.rm = TRUE) else
        NA_real_, passed = passed
    )
    if (passed || current == max_nsim) break
    next_count <- min(max_nsim, as.integer(ceiling(current * growth)))
    if (next_count <= current) next_count <- min(max_nsim, current + 1L)
    current <- next_count
  }
  structure(list(
    result = result, history = do.call(rbind, history), mcse = mcse,
    successful = successful, parameter_pass = parameter_pass, passed = passed,
    target_mcse = target_mcse, min_success = min_success,
    initial_nsim = initial_nsim, max_nsim = max_nsim, growth = growth,
    base_seed = seed, workers = workers, call = match.call()
  ), class = "spmodels_simulation_precision")
}

#' @export
print.spmodels_simulation_precision <- function(x, ...) {
  cat("SPmodels precision-directed simulation recovery\n")
  cat("Replicate workers:", x$workers, "\n")
  cat("Target met:", x$passed, " Target bias MCSE:",
      format(x$target_mcse), "\n")
  print(x$history, row.names = FALSE)
  invisible(x)
}
