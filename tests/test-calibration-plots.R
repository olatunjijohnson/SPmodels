if (file.exists(file.path("R", "calibration-plots.R"))) {
  for (file in list.files("R", full.names = TRUE)) source(file)
} else {
  library(SPmodels)
}

expect_error <- function(expr) {
  value <- try(force(expr), silent = TRUE)
  stopifnot(inherits(value, "try-error"))
}

draw_matrix <- matrix(c(
  -2, -1,  0,
  -1,  0,  1,
   0,  1,  2,
   1,  2,  3,
   2,  3,  4
), nrow = 5L, byrow = TRUE)
draws <- structure(list(
  draws = draw_matrix, family = "gaussian", target = "observation",
  type = "response"
), class = "spmodels_prediction_draws")
truth <- c(0, 1.5, 5)
levels <- c(0.5, 0.8)
calibration <- sp_calibration(draws, truth, levels = levels,
                              confidence = 0.9)

# Independently reconstruct the type-8 intervals and Wilson bounds.
expected_coverage <- vapply(levels, function(level) {
  alpha <- (1 - level) / 2
  limits <- t(apply(draw_matrix, 2L, quantile,
                    probs = c(alpha, 1 - alpha), names = FALSE, type = 8))
  mean(truth >= limits[, 1L] & truth <= limits[, 2L])
}, numeric(1))
stopifnot(all.equal(calibration$calibration$observed, expected_coverage,
                    tolerance = 1e-14))
z <- qnorm(0.95)
n <- length(truth)
denominator <- 1 + z^2 / n
expected_centre <- (expected_coverage + z^2 / (2 * n)) / denominator
expected_half <- z * sqrt(expected_coverage * (1 - expected_coverage) / n +
                            z^2 / (4 * n^2)) / denominator
stopifnot(all.equal(calibration$calibration$lower,
                    pmax(0, expected_centre - expected_half),
                    tolerance = 1e-14))
stopifnot(all.equal(calibration$calibration$upper,
                    pmin(1, expected_centre + expected_half),
                    tolerance = 1e-14))
stopifnot(identical(calibration$calibration$locations, rep(3L, 2L)))

expect_error(sp_calibration(list(), truth))
expect_error(sp_calibration(draws, truth[-1L]))
expect_error(sp_calibration(draws, truth, levels = c(0.8, 0.8)))
expect_error(sp_calibration(draws, truth, levels = 1))

# Plot methods must be views over the already recorded, unrounded values.
cross_validation <- structure(list(predictions = data.frame(
  fold = c("A", "A", "B"), observed = c(1, 2, 3),
  predicted = c(1.1, 1.8, 3.2), error = c(-0.1, 0.2, -0.2)
)), class = "spmodels_cross_validation")
assessment <- structure(list(locations = data.frame(
  truth = truth, mean = colMeans(draw_matrix),
  error = colMeans(draw_matrix) - truth,
  lower = apply(draw_matrix, 2L, min), upper = apply(draw_matrix, 2L, max)
)), class = "spmodels_simulation_assessment")
recovery <- structure(list(
  estimates = matrix(c(0.9, 1.1, 1.0), ncol = 1L,
                     dimnames = list(NULL, "intercept")),
  standard_errors = matrix(rep(0.1, 3L), ncol = 1L,
                           dimnames = list(NULL, "intercept")),
  truth = c(intercept = 1), level = 0.95
), class = "spmodels_simulation_recovery")

plot_file <- tempfile(fileext = ".pdf")
grDevices::pdf(plot_file)
stopifnot(identical(plot(calibration), calibration))
for (type in c("observed", "error", "fold"))
  stopifnot(identical(plot(cross_validation, type = type), cross_validation))
for (type in c("prediction", "error"))
  stopifnot(identical(plot(assessment, type = type), assessment))
stopifnot(identical(plot(recovery, parameter = "intercept"), recovery))
grDevices::dev.off()
stopifnot(file.exists(plot_file), file.info(plot_file)$size > 0)
unlink(plot_file)
expect_error(plot(recovery, parameter = "missing"))

cat("Calibration and validation/simulation display tests passed.\n")
