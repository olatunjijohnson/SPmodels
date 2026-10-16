if (file.exists(file.path("R", "spacetime.R"))) {
  for (file in list.files("R", full.names = TRUE)) source(file)
} else library(SPmodels)

expect_error <- function(expr) stopifnot(inherits(try(force(expr), silent = TRUE),
                                                  "try-error"))
coordinates <- cbind(east = c(0, 1, 0, 2), north = c(0, 0, 1, 1))
times <- c(0, 0, 2, 3)
distance <- as.matrix(dist(coordinates))
time_lag <- abs(outer(times, times, "-"))
correlation <- spacetime_cor(distance, time_lag, spatial_scale = 1.4,
                             smoothness = .5, temporal_scale = 2)
expected <- exp(-distance / 1.4) * exp(-time_lag / 2)
stopifnot(all.equal(correlation, expected, tolerance = 1e-14),
          all(diag(correlation) == 1),
          isSymmetric(correlation),
          min(eigen(correlation, symmetric = TRUE, only.values = TRUE)$values) >
            -1e-12,
          all.equal(spacetime_cor(distance, time_lag * 3, 1.4, .5, 6),
                    correlation, tolerance = 1e-14),
          all.equal(spacetime_cor(distance, matrix(0, 4, 4), 1.4, .5, 2),
                    exp(-distance / 1.4), tolerance = 1e-14))
expect_error(spacetime_cor(distance, time_lag[-1, ], 1, .5, 1))
expect_error(spacetime_cor(distance, time_lag, 1, .5, 0))
expect_error(spacetime_cor(distance, replace(time_lag, 1, -1), 1, .5, 1))
cat("Space-time correlation tests passed.\n")
