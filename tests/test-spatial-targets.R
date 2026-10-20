if (file.exists(file.path("R", "spatial-targets.R"))) {
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

# Target summaries are reconstructed directly from a fixed joint-draw matrix.
draw_matrix <- matrix(c(
  0.1, 0.2, 0.4, 0.5,
  0.3, 0.4, 0.6, 0.7,
  0.8, 0.7, 0.9, 1.0
), nrow = 4, ncol = 3,
dimnames = list(NULL, c("a", "b", "c")))
draw_object <- structure(list(
  draws = draw_matrix, family = "binomial", target = "latent",
  type = "response"
), class = "spmodels_prediction_draws")
point_target <- sp_predict_target(draw_object, threshold = 0.55,
                                  probs = c(0.25, 0.5, 0.75))
stopifnot(close(point_target$summary$mean, colMeans(draw_matrix)),
          close(point_target$summary$sd, apply(draw_matrix, 2, sd)),
          close(point_target$summary$exceedance_probability,
                colMeans(draw_matrix > 0.55)),
          identical(as.data.frame(point_target), point_target$summary))

weights <- c(1, 3, 2)
groups <- c("west", "west", "east")
group_target <- sp_predict_target(draw_object, weights = weights,
                                  groups = groups, threshold = 0.5)
west_expected <- 0.25 * draw_matrix[, 1] + 0.75 * draw_matrix[, 2]
east_expected <- draw_matrix[, 3]
stopifnot(close(group_target$draws,
                cbind(west = west_expected, east = east_expected)),
          close(colSums(group_target$weights), c(1, 1)))
expect_error(sp_predict_target(draw_object, groups = groups),
             "requires `weights`")
expect_error(sp_predict_target(draw_object, weights = c(0, 0, 1),
                               groups = groups), "positive total")

if (requireNamespace("sf", quietly = TRUE)) {
  points <- sf::st_as_sf(
    data.frame(x = c(0, 2, 0, 2), y = c(0, 0, 2, 2)),
    coords = c("x", "y"), crs = 3857
  )
  grid <- sp_grid(points, cellsize = 0.5)
  stopifnot(inherits(grid, "sf"), nrow(grid) == 16L,
            identical(grid$grid_id,
                      sprintf("grid_%06d", seq_len(nrow(grid)))))
  coordinates <- sf::st_coordinates(grid)
  ordering <- order(-coordinates[, 2], coordinates[, 1])
  stopifnot(identical(ordering, seq_len(nrow(grid))))
  polygons <- sf::st_sf(
    region = c("west", "east"),
    geometry = sf::st_sfc(
      sf::st_polygon(list(matrix(c(0, 0, 1, 0, 1, 2, 0, 2, 0, 0),
                                 ncol = 2, byrow = TRUE))),
      sf::st_polygon(list(matrix(c(1, 0, 2, 0, 2, 2, 1, 2, 1, 0),
                                 ncol = 2, byrow = TRUE))),
      crs = 3857
    )
  )
  grid$population <- seq_len(nrow(grid))
  polygon_weights <- sp_polygon_weights(
    grid, polygons, region = "region", weight = "population"
  )
  grid_draws <- structure(list(
    draws = matrix(seq_len(4 * nrow(grid)), nrow = 4,
                   dimnames = list(NULL, grid$grid_id)),
    family = "gaussian", target = "latent", type = "response"
  ), class = "spmodels_prediction_draws")
  polygon_target <- sp_predict_target(grid_draws, weights = polygon_weights)
  direct_west <- grid_draws$draws[, polygon_weights$index[
    polygon_weights$groups == "west"], drop = FALSE] %*%
    (polygon_weights$weights[polygon_weights$groups == "west"] /
       sum(polygon_weights$weights[polygon_weights$groups == "west"]))
  stopifnot(nrow(polygon_target$summary) == 2L,
            close(polygon_target$draws[, "west"], as.vector(direct_west)))
  expect_error(sp_grid(sf::st_transform(points, 4326), 0.01), "projected")
  expect_error(sp_polygon_weights(sf::st_transform(grid, 4326), polygons,
                                  "region"), "same known projected CRS")

  # Exact polygon-cell intersections split cells across regions and use
  # intersection area (optionally times density) inside each joint draw.
  cells <- sf::st_sf(
    cell_id = c("cell_1", "cell_2"), density = c(2, 1),
    geometry = sf::st_sfc(
      sf::st_polygon(list(matrix(c(0,0, 1,0, 1,1, 0,1, 0,0),
                                 ncol = 2, byrow = TRUE))),
      sf::st_polygon(list(matrix(c(1,0, 2,0, 2,1, 1,1, 1,0),
                                 ncol = 2, byrow = TRUE))),
      crs = 3857
    )
  )
  split_regions <- sf::st_sf(
    region = c("west", "east"),
    geometry = sf::st_sfc(
      sf::st_polygon(list(matrix(c(0,0, .5,0, .5,1, 0,1, 0,0),
                                 ncol = 2, byrow = TRUE))),
      sf::st_polygon(list(matrix(c(.5,0, 2,0, 2,1, .5,1, .5,0),
                                 ncol = 2, byrow = TRUE))),
      crs = 3857
    )
  )
  area_weights <- sp_polygon_area_weights(
    cells, split_regions, region = "region"
  )
  stopifnot(identical(area_weights$index, c(1L, 1L, 2L)),
            close(area_weights$intersection_area, c(.5, .5, 1)),
            close(area_weights$cell_fraction, c(.5, .5, 1)),
            identical(area_weights$groups, c("west", "east", "east")))
  cell_draws <- structure(list(
    draws = matrix(c(1,2,3,4, 10,20,30,40), nrow = 4,
                   dimnames = list(NULL, cells$cell_id)),
    family = "gaussian", target = "latent", type = "response"
  ), class = "spmodels_prediction_draws")
  area_target <- sp_predict_target(cell_draws, weights = area_weights)
  stopifnot(close(area_target$draws[, "west"], cell_draws$draws[, 1]),
            close(area_target$draws[, "east"],
                  cell_draws$draws[, 1] / 3 +
                    2 * cell_draws$draws[, 2] / 3))
  mapped_target <- sp_target_geometry(
    area_target, split_regions[c(2, 1), ], id = "region"
  )
  stopifnot(inherits(mapped_target, "spmodels_target_sf"),
            identical(mapped_target$target, c("east", "west")),
            close(mapped_target$mean,
                  area_target$summary$mean[match(c("east", "west"),
                                                  area_target$summary$target)]),
            identical(attr(mapped_target, "prediction_family"), "gaussian"))
  expect_error(sp_target_geometry(area_target, split_regions[1, ],
                                  id = "region"), "match exactly")
  density_weights <- sp_polygon_area_weights(
    cells, split_regions, region = "region", density = "density"
  )
  density_target <- sp_predict_target(cell_draws, weights = density_weights)
  stopifnot(close(density_target$draws[, "east"],
                  .5 * cell_draws$draws[, 1] +
                    .5 * cell_draws$draws[, 2]))

  regions_with_empty <- rbind(
    split_regions,
    sf::st_sf(region = "empty", geometry = sf::st_sfc(
      sf::st_polygon(list(matrix(c(3,0, 4,0, 4,1, 3,1, 3,0),
                                 ncol = 2, byrow = TRUE))), crs = 3857
    ))
  )
  expect_error(sp_polygon_area_weights(
    cells, regions_with_empty, "region"
  ), "Regions without")
  dropped_empty <- sp_polygon_area_weights(
    cells, regions_with_empty, "region", empty_regions = "drop"
  )
  stopifnot(identical(dropped_empty$empty_regions, "empty"))

  outside_cell <- cells[2, ]
  sf::st_geometry(outside_cell) <- sf::st_set_crs(
    sf::st_geometry(outside_cell) + c(10, 0), sf::st_crs(cells)
  )
  cells_with_outside <- rbind(cells[1, ], outside_cell)
  expect_error(sp_polygon_area_weights(
    cells_with_outside, split_regions, "region"
  ), "no positive-area")
  dropped_cell <- sp_polygon_area_weights(
    cells_with_outside, split_regions, "region", unmatched = "drop",
    empty_regions = "drop"
  )
  stopifnot(identical(dropped_cell$unmatched, 2L))

  collinear <- sf::st_as_sf(data.frame(x = 0:2, y = 0:2),
                            coords = c("x", "y"), crs = 3857)
  expect_error(sp_grid(collinear, 0.5), "positive polygonal area")
}

cat("Spatial grid and prediction-target tests passed.\n")
