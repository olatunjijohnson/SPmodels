#' Summarise pairwise spatial distances
#'
#' Calculate an auditable summary of all unordered pairs, with an explicit
#' policy for repeated coordinate locations.
#'
#' @param data A fitted `spmodels_fit`, a projected point `sf` object, or a data
#'   frame containing coordinate columns.
#' @param coords A one-sided two-column coordinate formula for ordinary data
#'   frames. It is ignored for a fitted model and must be omitted for `sf`.
#' @param distance_units Required coordinate units for an ordinary data frame.
#'   Fitted objects and projected `sf` inputs supply their recorded units.
#' @param duplicates Whether zero-distance pairs are included, excluded, or
#'   treated as an error.
#' @param probs Finite quantile probabilities between zero and one.
#'
#' @return An object of class `spmodels_distance_summary` containing retained
#'   pair distances, a one-row summary and the duplicate policy.
#' @export
sp_distance_summary <- function(data, coords = NULL, distance_units = NULL,
                                duplicates = c("include", "exclude", "error"),
                                probs = c(0, 0.25, 0.5, 0.75, 1)) {
  duplicates <- match.arg(duplicates)
  if (!is.numeric(probs) || !length(probs) || anyNA(probs) ||
      any(!is.finite(probs)) || any(probs < 0 | probs > 1)) {
    stop("`probs` must contain finite probabilities between zero and one.",
         call. = FALSE)
  }
  probs <- sort(unique(probs))
  if (inherits(data, "spmodels_fit")) {
    if (!isTRUE(data$spatial) || is.null(data$coordinates)) {
      stop("`data` must be a fitted spatial model.", call. = FALSE)
    }
    if (!is.null(coords) || !is.null(distance_units)) {
      stop("Do not supply `coords` or `distance_units` for a fitted model.",
           call. = FALSE)
    }
    coordinates <- data$coordinates
    units <- data$distance_units
  } else if (inherits(data, "sf")) {
    if (!is.null(coords) || !is.null(distance_units)) {
      stop("Do not supply `coords` or `distance_units` for `sf` input.",
           call. = FALSE)
    }
    coordinate_info <- .sp_coordinates(NULL, data)
    coordinates <- coordinate_info$values
    units <- coordinate_info$inferred_units
  } else {
    if (!is.data.frame(data)) {
      stop("`data` must be a fitted model, point `sf`, or data frame.",
           call. = FALSE)
    }
    coordinate_info <- .sp_coordinates(coords, data)
    coordinates <- coordinate_info$values
    if (!is.character(distance_units) || length(distance_units) != 1L ||
        is.na(distance_units) || !nzchar(trimws(distance_units))) {
      stop("Ordinary data-frame input requires explicit `distance_units`.",
           call. = FALSE)
    }
    units <- trimws(distance_units)
  }
  if (nrow(coordinates) < 2L) {
    stop("At least two coordinate rows are required.", call. = FALSE)
  }
  distance_matrix <- .sp_dist(coordinates)
  all_distances <- distance_matrix[lower.tri(distance_matrix)]
  same_location <- outer(seq_len(nrow(coordinates)),
                         seq_len(nrow(coordinates)), Vectorize(function(i, j) {
    coordinates[i, 1L] == coordinates[j, 1L] &&
      coordinates[i, 2L] == coordinates[j, 2L]
  }))
  repeated <- same_location[lower.tri(same_location)]
  duplicate_pairs <- sum(repeated)
  if (duplicate_pairs && duplicates == "error") {
    stop(sprintf("Found %d zero-distance coordinate pair%s.", duplicate_pairs,
                 if (duplicate_pairs == 1L) "" else "s"), call. = FALSE)
  }
  retained <- if (duplicates == "exclude") all_distances[!repeated] else
    all_distances
  if (!length(retained)) {
    stop("No pair distances remain under the duplicate-location policy.",
         call. = FALSE)
  }
  quantiles <- stats::quantile(retained, probs = probs, names = FALSE,
                               type = 8)
  quantile_names <- paste0("q_", format(probs, trim = TRUE,
                                          scientific = FALSE))
  summary <- data.frame(
    locations = nrow(coordinates), total_pairs = length(all_distances),
    duplicate_pairs = duplicate_pairs, included_pairs = length(retained),
    mean = mean(retained), sd = if (length(retained) > 1L)
      stats::sd(retained) else NA_real_,
    stringsAsFactors = FALSE, check.names = FALSE
  )
  summary <- cbind(summary, stats::setNames(as.list(quantiles),
                                            quantile_names))
  structure(list(
    distances = retained, summary = summary, units = units,
    duplicate_policy = duplicates, probs = probs, call = match.call()
  ), class = "spmodels_distance_summary")
}

#' @export
print.spmodels_distance_summary <- function(x, ...) {
  cat("SPmodels pairwise distance summary\n")
  cat("Units:", if (is.null(x$units)) "unknown" else x$units,
      " Duplicate policy:", x$duplicate_policy, "\n")
  print(x$summary, row.names = FALSE)
  invisible(x)
}

.sp_utm_zone <- function(longitude, latitude) {
  zone <- pmin(60L, pmax(1L, floor((longitude + 180) / 6) + 1L))
  norway <- latitude >= 56 & latitude < 64 & longitude >= 3 & longitude < 12
  zone[norway] <- 32L
  svalbard <- latitude >= 72 & latitude < 84
  zone[svalbard & longitude >= 0 & longitude < 9] <- 31L
  zone[svalbard & longitude >= 9 & longitude < 21] <- 33L
  zone[svalbard & longitude >= 21 & longitude < 33] <- 35L
  zone[svalbard & longitude >= 33 & longitude < 42] <- 37L
  as.integer(zone)
}

#' Choose a projected UTM coordinate reference system
#'
#' Select a WGS84 UTM zone and hemisphere from the complete coordinates of an
#' `sf` object. Zone and equator crossings are handled by explicit policies;
#' no projection is silently selected for a multi-zone study by default.
#'
#' @param data A non-empty `sf` or `sfc` object with a known CRS.
#' @param zone_policy What to do when coordinates occupy more than one UTM
#'   zone: error, or use the bounding-box centre.
#' @param hemisphere_policy What to do when coordinates span the equator:
#'   error, or use the bounding-box centre.
#'
#' @return An object of class `spmodels_crs_choice` containing the selected
#'   `sf::crs`, EPSG code, zone, hemisphere and audit metadata.
#' @export
sp_utm_crs <- function(data, zone_policy = c("error", "centre"),
                       hemisphere_policy = c("error", "centre")) {
  zone_policy <- match.arg(zone_policy)
  hemisphere_policy <- match.arg(hemisphere_policy)
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package `sf` is required by `sp_utm_crs()`.", call. = FALSE)
  }
  if (!inherits(data, c("sf", "sfc")) || length(sf::st_geometry(data)) == 0L) {
    stop("`data` must be a non-empty `sf` or `sfc` object.", call. = FALSE)
  }
  if (is.na(sf::st_crs(data))) {
    stop("`data` must have a known coordinate reference system.",
         call. = FALSE)
  }
  geographic <- sf::st_transform(sf::st_geometry(data), 4326)
  if (any(sf::st_is_empty(geographic))) {
    stop("`data` must not contain empty geometries.", call. = FALSE)
  }
  coordinates <- sf::st_coordinates(geographic)[, 1:2, drop = FALSE]
  longitude <- coordinates[, 1L]
  latitude <- coordinates[, 2L]
  if (any(!is.finite(longitude)) || any(!is.finite(latitude))) {
    stop("Transformed longitude and latitude must be finite.", call. = FALSE)
  }
  if (any(latitude < -80 | latitude > 84)) {
    stop("UTM is defined here only between 80 degrees south and 84 degrees north.",
         call. = FALSE)
  }
  if (diff(range(longitude)) > 180) {
    stop("Geometries crossing the antimeridian require an explicit local projection.",
         call. = FALSE)
  }
  zones <- sort(unique(.sp_utm_zone(longitude, latitude)))
  bbox <- sf::st_bbox(geographic)
  centre_longitude <- mean(c(bbox[["xmin"]], bbox[["xmax"]]))
  centre_latitude <- mean(c(bbox[["ymin"]], bbox[["ymax"]]))
  if (length(zones) > 1L && zone_policy == "error") {
    stop(sprintf("Coordinates occupy multiple UTM zones: %s.",
                 paste(zones, collapse = ", ")), call. = FALSE)
  }
  spans_equator <- min(latitude) < 0 && max(latitude) >= 0
  if (spans_equator && hemisphere_policy == "error") {
    stop("Coordinates span the equator; select a hemisphere policy explicitly.",
         call. = FALSE)
  }
  zone <- if (length(zones) == 1L) zones else
    .sp_utm_zone(centre_longitude, centre_latitude)
  hemisphere <- if (centre_latitude >= 0) "north" else "south"
  epsg <- if (hemisphere == "north") 32600L + zone else 32700L + zone
  structure(list(
    crs = sf::st_crs(epsg), epsg = epsg, zone = zone,
    hemisphere = hemisphere, centre_longitude = centre_longitude,
    centre_latitude = centre_latitude, zones_present = zones,
    spans_equator = spans_equator, zone_policy = zone_policy,
    hemisphere_policy = hemisphere_policy, call = match.call()
  ), class = "spmodels_crs_choice")
}

#' @export
print.spmodels_crs_choice <- function(x, ...) {
  cat("SPmodels UTM CRS choice\n")
  cat("EPSG:", x$epsg, " Zone:", x$zone,
      " Hemisphere:", x$hemisphere, "\n")
  if (length(x$zones_present) > 1L) {
    cat("Input zones:", paste(x$zones_present, collapse = ", "),
        " (centre policy used)\n")
  }
  invisible(x)
}
