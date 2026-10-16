if (file.exists(file.path("R", "geometry-utils.R"))) {
  for (file in list.files("R", full.names = TRUE)) source(file)
} else {
  library(SPmodels)
}

close <- function(a, b, tolerance = 1e-12) {
  isTRUE(all.equal(a, b, tolerance = tolerance, check.attributes = FALSE))
}
expect_error <- function(expression, pattern) {
  error <- tryCatch({ force(expression); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}

locations <- data.frame(
  east = c(0, 3, 0, 0), north = c(0, 4, 0, 4)
)
all_pairs <- c(5, 0, 5, 4, 3, 4)
included <- sp_distance_summary(
  locations, coords = ~ east + north, distance_units = "km",
  duplicates = "include", probs = c(0, 0.5, 1)
)
excluded <- sp_distance_summary(
  locations, coords = ~ east + north, distance_units = "km",
  duplicates = "exclude", probs = c(0, 0.5, 1)
)
stopifnot(close(sort(included$distances), sort(all_pairs)),
          included$summary$total_pairs == choose(4, 2),
          included$summary$duplicate_pairs == 1L,
          included$summary$included_pairs == 6L,
          excluded$summary$included_pairs == 5L,
          excluded$summary$q_0 == 3,
          excluded$summary$q_1 == 5,
          identical(included$units, "km"))
expect_error(sp_distance_summary(
  locations, ~ east + north, "km", duplicates = "error"
), "zero-distance")
expect_error(sp_distance_summary(locations, ~ east + north),
             "distance_units")

fit <- sp_fit(
  c(1, 2, 1.5, 2.5) ~ east, locations,
  coords = ~ east + north,
  control = sp_control(
    spatial_variance = 1, scale = 2, nugget_variance = 0.2,
    estimate = FALSE, distance_units = "km"
  )
)
fit_distances <- sp_distance_summary(fit)
stopifnot(close(fit_distances$distances, included$distances),
          identical(fit_distances$units, "km"))

if (requireNamespace("sf", quietly = TRUE)) {
  manchester <- sf::st_as_sf(
    data.frame(lon = c(-2.4, -2.1), lat = c(53.3, 53.6)),
    coords = c("lon", "lat"), crs = 4326
  )
  north_choice <- sp_utm_crs(manchester)
  stopifnot(north_choice$zone == 30L,
            north_choice$hemisphere == "north",
            north_choice$epsg == 32630L,
            isTRUE(sf::st_crs(north_choice$crs) == sf::st_crs(32630)))

  southern <- sf::st_as_sf(
    data.frame(lon = c(18.2, 18.6), lat = c(-34.2, -33.8)),
    coords = c("lon", "lat"), crs = 4326
  )
  south_choice <- sp_utm_crs(southern)
  stopifnot(south_choice$zone == 34L,
            south_choice$hemisphere == "south",
            south_choice$epsg == 32734L)

  # Exact six-degree boundaries use the eastern zone; longitude 180 is capped
  # at zone 60. Norway and Svalbard use their standard exceptional zones.
  boundary <- sf::st_as_sf(
    data.frame(lon = c(0, 5), lat = c(10, 10)),
    coords = c("lon", "lat"), crs = 4326
  )
  stopifnot(sp_utm_crs(boundary)$zone == 31L)
  norway <- sf::st_as_sf(data.frame(lon = 6, lat = 60),
                         coords = c("lon", "lat"), crs = 4326)
  svalbard <- sf::st_as_sf(data.frame(lon = 15, lat = 78),
                           coords = c("lon", "lat"), crs = 4326)
  dateline <- sf::st_as_sf(data.frame(lon = 180, lat = 10),
                           coords = c("lon", "lat"), crs = 4326)
  stopifnot(sp_utm_crs(norway)$zone == 32L,
            sp_utm_crs(svalbard)$zone == 33L,
            sp_utm_crs(dateline)$zone == 60L)

  multi_zone <- sf::st_as_sf(
    data.frame(lon = c(-1, 7), lat = c(50, 50)),
    coords = c("lon", "lat"), crs = 4326
  )
  expect_error(sp_utm_crs(multi_zone), "multiple UTM zones")
  centre_choice <- sp_utm_crs(multi_zone, zone_policy = "centre")
  stopifnot(length(centre_choice$zones_present) == 2L,
            centre_choice$zone == 31L)

  equator <- sf::st_as_sf(
    data.frame(lon = c(20, 20), lat = c(-1, 1)),
    coords = c("lon", "lat"), crs = 4326
  )
  expect_error(sp_utm_crs(equator), "span the equator")
  equator_choice <- sp_utm_crs(
    equator, hemisphere_policy = "centre"
  )
  stopifnot(equator_choice$hemisphere == "north",
            isTRUE(equator_choice$spans_equator))

  polar <- sf::st_as_sf(data.frame(lon = 0, lat = 85),
                        coords = c("lon", "lat"), crs = 4326)
  expect_error(sp_utm_crs(polar), "only between")
  antimeridian <- sf::st_as_sf(
    data.frame(lon = c(179, -179), lat = c(10, 10)),
    coords = c("lon", "lat"), crs = 4326
  )
  expect_error(sp_utm_crs(antimeridian), "antimeridian")

  projected <- sf::st_transform(manchester, north_choice$crs)
  projected_distances <- sp_distance_summary(projected)
  stopifnot(identical(projected_distances$units, "metre"),
            projected_distances$summary$total_pairs == 1L)
}

cat("Distance and CRS utility tests passed.\n")
