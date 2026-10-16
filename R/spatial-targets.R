#' Construct a regular prediction grid inside a projected boundary
#'
#' Create point centres on a regular grid and retain centres covered by a
#' polygon boundary. When only point data are supplied, their convex hull is
#' used as the boundary.
#'
#' @param data A projected `sf` object containing points or polygonal geometry.
#' @param cellsize One positive grid spacing or two positive x/y spacings, in
#'   the coordinate units of `data`.
#' @param boundary Optional projected polygon `sf` or `sfc` object. It must use
#'   the same CRS as `data`.
#'
#' @return A point `sf` object with stable row-major `grid_id` values.
#' @export
sp_grid <- function(data, cellsize, boundary = NULL) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package `sf` is required by `sp_grid()`.", call. = FALSE)
  }
  if (!inherits(data, "sf") || !nrow(data)) {
    stop("`data` must be a non-empty `sf` object.", call. = FALSE)
  }
  crs <- sf::st_crs(data)
  if (is.na(crs) || isTRUE(sf::st_is_longlat(data))) {
    stop("`data` must have a known projected coordinate reference system.",
         call. = FALSE)
  }
  if (!is.numeric(cellsize) || !length(cellsize) || length(cellsize) > 2L ||
      anyNA(cellsize) || any(!is.finite(cellsize)) || any(cellsize <= 0)) {
    stop("`cellsize` must contain one or two finite positive numbers.",
         call. = FALSE)
  }
  if (length(cellsize) == 1L) cellsize <- rep(cellsize, 2L)
  geometry_type <- unique(as.character(sf::st_geometry_type(data)))
  if (is.null(boundary)) {
    if (all(geometry_type %in% c("POINT", "MULTIPOINT"))) {
      boundary_geometry <- sf::st_convex_hull(sf::st_union(sf::st_geometry(data)))
    } else if (all(geometry_type %in% c("POLYGON", "MULTIPOLYGON"))) {
      boundary_geometry <- sf::st_union(sf::st_geometry(data))
    } else {
      stop("`data` must contain only point or only polygonal geometry.",
           call. = FALSE)
    }
  } else {
    if (!inherits(boundary, c("sf", "sfc"))) {
      stop("`boundary` must be an `sf` or `sfc` polygon object.", call. = FALSE)
    }
    if (!isTRUE(sf::st_crs(boundary) == crs)) {
      stop("`boundary` and `data` must use the same CRS.", call. = FALSE)
    }
    boundary_geometry <- sf::st_union(sf::st_geometry(boundary))
  }
  boundary_type <- as.character(sf::st_geometry_type(boundary_geometry))
  if (!all(boundary_type %in% c("POLYGON", "MULTIPOLYGON")) ||
      as.numeric(sf::st_area(boundary_geometry)) <= 0) {
    stop("The grid boundary must have positive polygonal area.", call. = FALSE)
  }
  centres <- sf::st_make_grid(boundary_geometry, cellsize = cellsize,
                              what = "centers", square = TRUE)
  covered <- lengths(sf::st_intersects(centres, boundary_geometry)) > 0L
  centres <- centres[covered]
  if (!length(centres)) {
    stop("No grid centres fall inside the boundary at this `cellsize`.",
         call. = FALSE)
  }
  coordinates <- sf::st_coordinates(centres)
  ordering <- order(-coordinates[, 2L], coordinates[, 1L])
  centres <- centres[ordering]
  sf::st_sf(grid_id = sprintf("grid_%06d", seq_along(centres)),
            geometry = centres)
}

#' Assign point prediction weights to polygon regions
#'
#' Match projected point locations to polygon identifiers and attach equal or
#' supplied non-negative point weights for regional draw aggregation.
#'
#' @param points A projected point `sf` object in prediction-row order.
#' @param polygons A projected polygon `sf` object using the same CRS.
#' @param region One column name in `polygons` containing non-missing region
#'   identifiers.
#' @param weight Optional numeric point weights or one numeric column name in
#'   `points`. Equal raw weights are used by default.
#' @param unmatched Either fail for points outside every polygon or drop them
#'   from subsequent aggregation.
#'
#' @return An object of class `spmodels_polygon_weights` accepted directly as
#'   the `weights` argument of [sp_predict_target()].
#' @export
sp_polygon_weights <- function(points, polygons, region, weight = NULL,
                               unmatched = c("error", "drop")) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package `sf` is required by `sp_polygon_weights()`.", call. = FALSE)
  }
  unmatched <- match.arg(unmatched)
  if (!inherits(points, "sf") || !nrow(points) ||
      any(!as.character(sf::st_geometry_type(points)) %in% "POINT")) {
    stop("`points` must be a non-empty point `sf` object.", call. = FALSE)
  }
  if (!inherits(polygons, "sf") || !nrow(polygons) ||
      any(!as.character(sf::st_geometry_type(polygons)) %in%
          c("POLYGON", "MULTIPOLYGON"))) {
    stop("`polygons` must be a non-empty polygon `sf` object.", call. = FALSE)
  }
  point_crs <- sf::st_crs(points)
  polygon_crs <- sf::st_crs(polygons)
  if (is.na(point_crs) || isTRUE(sf::st_is_longlat(points)) ||
      !isTRUE(point_crs == polygon_crs)) {
    stop("Points and polygons must use the same known projected CRS.",
         call. = FALSE)
  }
  if (!is.character(region) || length(region) != 1L || is.na(region) ||
      !region %in% names(polygons)) {
    stop("`region` must name one column in `polygons`.", call. = FALSE)
  }
  region_values <- polygons[[region]]
  if (anyNA(region_values)) {
    stop("Polygon region identifiers must be non-missing.", call. = FALSE)
  }
  point_weights <- if (is.null(weight)) {
    rep(1, nrow(points))
  } else if (is.character(weight) && length(weight) == 1L &&
             !is.na(weight) && weight %in% names(points)) {
    points[[weight]]
  } else weight
  if (!is.numeric(point_weights) || length(point_weights) != nrow(points) ||
      anyNA(point_weights) || any(!is.finite(point_weights)) ||
      any(point_weights < 0)) {
    stop("Point weights must be finite non-negative values matching `points`.",
         call. = FALSE)
  }
  membership <- sf::st_within(points, polygons)
  assigned <- lapply(membership, function(index) unique(as.character(
    region_values[index]
  )))
  overlap <- lengths(assigned) > 1L
  if (any(overlap)) {
    stop("A point belongs to polygons with multiple region identifiers.",
         call. = FALSE)
  }
  outside <- lengths(assigned) == 0L
  if (any(outside) && unmatched == "error") {
    stop("Some points fall outside every polygon; use `unmatched = \"drop\"` to exclude them.",
         call. = FALSE)
  }
  index <- which(!outside)
  groups <- vapply(assigned[index], `[[`, character(1), 1L)
  structure(list(
    index = index, weights = as.numeric(point_weights[index]), groups = groups,
    region = region, unmatched = which(outside), crs = point_crs,
    point_count = nrow(points), call = match.call()
  ), class = "spmodels_polygon_weights")
}

#' Construct exact fractional cell-area polygon weights
#'
#' Intersect projected polygon prediction cells with target regions and use
#' exact intersection areas, optionally multiplied by a cell-level density, as
#' raw aggregation weights. A cell may contribute to several regions.
#'
#' @param cells A non-empty projected polygon `sf` object in prediction-row
#'   order.
#' @param polygons A projected polygon `sf` object using the same CRS.
#' @param region One column name in `polygons` containing non-missing region
#'   identifiers.
#' @param density Optional finite non-negative numeric value per cell, or one
#'   numeric column name in `cells`. Weights are area times density.
#' @param unmatched Policy for cells with no positive-area region intersection.
#' @param empty_regions Policy for declared region identifiers with no
#'   positive-area cell intersection.
#'
#' @return An `spmodels_polygon_weights` object accepted by
#'   [sp_predict_target()], with intersection areas and cell fractions retained.
#' @export
sp_polygon_area_weights <- function(
    cells, polygons, region, density = NULL,
    unmatched = c("error", "drop"),
    empty_regions = c("error", "drop")) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package `sf` is required by `sp_polygon_area_weights()`.",
         call. = FALSE)
  }
  unmatched <- match.arg(unmatched)
  empty_regions <- match.arg(empty_regions)
  polygon_types <- c("POLYGON", "MULTIPOLYGON")
  if (!inherits(cells, "sf") || !nrow(cells) ||
      any(!as.character(sf::st_geometry_type(cells)) %in% polygon_types) ||
      any(sf::st_is_empty(cells))) {
    stop("`cells` must be a non-empty polygon `sf` object without empty geometries.",
         call. = FALSE)
  }
  if (!inherits(polygons, "sf") || !nrow(polygons) ||
      any(!as.character(sf::st_geometry_type(polygons)) %in% polygon_types) ||
      any(sf::st_is_empty(polygons))) {
    stop("`polygons` must be a non-empty polygon `sf` object without empty geometries.",
         call. = FALSE)
  }
  cell_crs <- sf::st_crs(cells)
  if (is.na(cell_crs) || isTRUE(sf::st_is_longlat(cells)) ||
      !isTRUE(cell_crs == sf::st_crs(polygons))) {
    stop("Cells and polygons must use the same known projected CRS.",
         call. = FALSE)
  }
  if (!is.character(region) || length(region) != 1L || is.na(region) ||
      !region %in% names(polygons)) {
    stop("`region` must name one column in `polygons`.", call. = FALSE)
  }
  region_values <- as.character(polygons[[region]])
  if (anyNA(region_values) || any(!nzchar(region_values))) {
    stop("Polygon region identifiers must be non-missing and non-empty.",
         call. = FALSE)
  }
  overlaps <- sf::st_overlaps(polygons)
  if (any(lengths(overlaps) > 0L)) {
    stop("Target polygons must not overlap each other.",
         call. = FALSE)
  }
  cell_density <- if (is.null(density)) {
    rep(1, nrow(cells))
  } else if (is.character(density) && length(density) == 1L &&
             !is.na(density) && density %in% names(cells)) {
    cells[[density]]
  } else density
  if (!is.numeric(cell_density) || length(cell_density) != nrow(cells) ||
      anyNA(cell_density) || any(!is.finite(cell_density)) ||
      any(cell_density < 0)) {
    stop("Cell densities must be finite non-negative values matching `cells`.",
         call. = FALSE)
  }
  cell_area <- as.numeric(sf::st_area(cells))
  if (any(!is.finite(cell_area)) || any(cell_area <= 0)) {
    stop("Every prediction cell must have positive finite area.", call. = FALSE)
  }
  cell_input <- sf::st_sf(.cell_index = seq_len(nrow(cells)),
                          geometry = sf::st_geometry(cells))
  region_input <- sf::st_sf(.region_id = region_values,
                            geometry = sf::st_geometry(polygons))
  intersection <- suppressWarnings(sf::st_intersection(cell_input,
                                                         region_input))
  if (nrow(intersection)) {
    intersection_area <- as.numeric(sf::st_area(intersection))
    positive <- is.finite(intersection_area) & intersection_area > 0
    intersection <- intersection[positive, , drop = FALSE]
    intersection_area <- intersection_area[positive]
  } else {
    intersection_area <- numeric()
  }
  if (!nrow(intersection)) {
    stop("No prediction cell has a positive-area region intersection.",
         call. = FALSE)
  }
  index <- intersection$.cell_index
  groups <- as.character(intersection$.region_id)
  matched_cells <- sort(unique(index))
  unmatched_index <- setdiff(seq_len(nrow(cells)), matched_cells)
  if (length(unmatched_index) && unmatched == "error") {
    stop("Some cells have no positive-area region intersection; use `unmatched = \"drop\"` to exclude them.",
         call. = FALSE)
  }
  declared_regions <- unique(region_values)
  empty <- setdiff(declared_regions, unique(groups))
  if (length(empty) && empty_regions == "error") {
    stop(sprintf("Regions without positive-area cell intersections: %s.",
                 paste(empty, collapse = ", ")), call. = FALSE)
  }
  fractions <- if (length(index)) intersection_area / cell_area[index] else
    numeric()
  raw_weights <- if (length(index))
    intersection_area * cell_density[index] else numeric()
  structure(list(
    index = as.integer(index), weights = raw_weights, groups = groups,
    region = region, unmatched = unmatched_index, empty_regions = empty,
    empty_region_policy = empty_regions, crs = cell_crs,
    point_count = nrow(cells), intersection_area = intersection_area,
    cell_fraction = fractions, allocation = "fractional_area",
    call = match.call()
  ), class = "spmodels_polygon_weights")
}

#' Summarise point, grid, or grouped targets from joint draws
#'
#' Calculate draw-based means, standard deviations, quantiles and optional
#' exceedance probabilities. Optional weights form a weighted mean within each
#' supplied group while preserving joint uncertainty.
#'
#' @param object An object returned by [sp_predict_draws()].
#' @param weights Optional finite non-negative weight for every prediction
#'   location. Weights are normalised within each group.
#' @param groups Optional group identifier for every location. It requires
#'   `weights`; the default is one aggregate target.
#' @param threshold Optional finite scalar used to calculate the probability
#'   that each target exceeds it.
#' @param probs Numeric quantile probabilities between zero and one.
#'
#' @return An object of class `spmodels_prediction_target` containing the
#'   target-level draws, summary table and complete target definition.
#' @export
sp_predict_target <- function(object, weights = NULL, groups = NULL,
                              threshold = NULL,
                              probs = c(0.025, 0.5, 0.975)) {
  if (!inherits(object, "spmodels_prediction_draws")) {
    stop("`object` must be returned by `sp_predict_draws()`.", call. = FALSE)
  }
  draws <- object$draws
  polygon_weights <- NULL
  if (inherits(weights, "spmodels_polygon_weights")) {
    if (!is.null(groups)) {
      stop("Do not supply `groups` when `weights` comes from `sp_polygon_weights()`.",
           call. = FALSE)
    }
    if (weights$point_count != ncol(draws)) {
      stop("Polygon weights do not match the prediction locations.",
           call. = FALSE)
    }
    polygon_weights <- weights
    draws <- draws[, weights$index, drop = FALSE]
    groups <- weights$groups
    weights <- weights$weights
  }
  locations <- ncol(draws)
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
  if (is.null(weights)) {
    if (!is.null(groups)) {
      stop("`groups` requires `weights`.", call. = FALSE)
    }
    target_draws <- draws
    target_id <- colnames(draws)
    if (is.null(target_id)) target_id <- paste0("location_", seq_len(locations))
    normalised_weights <- NULL
    group_values <- NULL
  } else {
    if (!is.numeric(weights) || length(weights) != locations || anyNA(weights) ||
        any(!is.finite(weights)) || any(weights < 0)) {
      stop("`weights` must be finite non-negative values matching prediction locations.",
           call. = FALSE)
    }
    if (is.null(groups)) groups <- rep("all", locations)
    if (length(groups) != locations || anyNA(groups)) {
      stop("`groups` must be non-missing and match prediction locations.",
           call. = FALSE)
    }
    group_values <- unique(as.character(groups))
    weight_matrix <- matrix(0, nrow = locations, ncol = length(group_values),
                            dimnames = list(colnames(draws), group_values))
    for (j in seq_along(group_values)) {
      selected <- as.character(groups) == group_values[j]
      total <- sum(weights[selected])
      if (total <= 0) {
        stop("Every group must have positive total weight.", call. = FALSE)
      }
      weight_matrix[selected, j] <- weights[selected] / total
    }
    target_draws <- draws %*% weight_matrix
    target_id <- group_values
    normalised_weights <- weight_matrix
  }
  quantiles <- t(apply(target_draws, 2L, stats::quantile,
                       probs = probs, names = FALSE, type = 8))
  quantile_names <- paste0("q_", format(probs, trim = TRUE,
                                         scientific = FALSE))
  colnames(quantiles) <- quantile_names
  summary <- data.frame(
    target = target_id, mean = colMeans(target_draws),
    sd = apply(target_draws, 2L, stats::sd), row.names = NULL,
    check.names = FALSE
  )
  summary <- cbind(summary, quantiles)
  if (!is.null(threshold)) {
    summary$exceedance_probability <- colMeans(target_draws > threshold)
  }
  structure(list(
    draws = target_draws, summary = summary, weights = normalised_weights,
    groups = group_values, threshold = threshold, probs = probs,
    polygon_weights = polygon_weights,
    family = object$family, target = object$target, type = object$type,
    call = match.call()
  ), class = "spmodels_prediction_target")
}

#' @export
print.spmodels_prediction_target <- function(x, ...) {
  cat("Prediction target summary\n")
  if (!is.null(x$approximation)) {
    cat("Approximation:", x$approximation, "with at most", x$neighbors,
        "predecessors\n")
    cat("Peak active point-draw vectors:", x$peak_active_locations, "\n")
  }
  print(x$summary, row.names = FALSE)
  invisible(x)
}

#' Attach prediction-target summaries to spatial geometry
#'
#' Join the unrounded numerical summary from [sp_predict_target()] to point or
#' polygon geometry without changing the calculations.
#'
#' @param object A prediction-target object.
#' @param geometry A non-empty `sf` object.
#' @param id Optional column name in `geometry` matched to
#'   `object$summary$target`. When omitted, rows must already have the same
#'   length and order.
#'
#' @return An `sf` object of class `spmodels_target_sf` containing the target
#'   summary and geometry.
#' @export
sp_target_geometry <- function(object, geometry, id = NULL) {
  if (!inherits(object, "spmodels_prediction_target")) {
    stop("`object` must be returned by `sp_predict_target()`.", call. = FALSE)
  }
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package `sf` is required by `sp_target_geometry()`.", call. = FALSE)
  }
  if (!inherits(geometry, "sf") || !nrow(geometry) ||
      any(sf::st_is_empty(geometry))) {
    stop("`geometry` must be a non-empty `sf` object without empty geometries.",
         call. = FALSE)
  }
  summary <- object$summary
  if (is.null(id)) {
    if (nrow(geometry) != nrow(summary)) {
      stop("Without `id`, geometry rows must match target rows in number and order.",
           call. = FALSE)
    }
    ordering <- seq_len(nrow(summary))
  } else {
    if (!is.character(id) || length(id) != 1L || is.na(id) ||
        !id %in% names(geometry)) {
      stop("`id` must name one column in `geometry`.", call. = FALSE)
    }
    geometry_id <- as.character(geometry[[id]])
    target_id <- as.character(summary$target)
    if (anyNA(geometry_id) || anyDuplicated(geometry_id) ||
        anyDuplicated(target_id)) {
      stop("Geometry and target identifiers must be non-missing and unique.",
           call. = FALSE)
    }
    ordering <- match(geometry_id, target_id)
    if (anyNA(ordering) || length(geometry_id) != length(target_id) ||
        any(!target_id %in% geometry_id)) {
      stop("Geometry and prediction target identifiers must match exactly.",
           call. = FALSE)
    }
  }
  attributes_to_add <- summary[ordering, , drop = FALSE]
  non_geometry <- sf::st_drop_geometry(geometry)
  duplicate_names <- intersect(names(non_geometry), names(attributes_to_add))
  if (length(duplicate_names)) {
    stop(sprintf("Geometry already contains target-summary column%s: %s.",
                 if (length(duplicate_names) == 1L) "" else "s",
                 paste(duplicate_names, collapse = ", ")), call. = FALSE)
  }
  out <- cbind(geometry, attributes_to_add)
  class(out) <- c("spmodels_target_sf", class(out))
  attr(out, "prediction_family") <- object$family
  attr(out, "prediction_target") <- object$target
  attr(out, "prediction_scale") <- object$type
  out
}

#' @export
plot.spmodels_target_sf <- function(x, value = "mean", ...) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !value %in% names(x) || !is.numeric(x[[value]])) {
    stop("`value` must name one numeric target-summary column.", call. = FALSE)
  }
  NextMethod("plot", x[value], ...)
  invisible(x)
}
