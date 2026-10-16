if (file.exists(file.path("R", "report-presentation.R"))) {
  for (file in list.files("R", full.names = TRUE)) source(file)
} else {
  library(SPmodels)
}

expect_error <- function(expr) {
  value <- try(force(expr), silent = TRUE)
  stopifnot(inherits(value, "try-error"))
}

raw <- data.frame(
  term = c("alpha", "beta", "gamma"),
  estimate = c(1.23456, -0.0004, NA_real_),
  p_value = c(0.12345, 0.0004, NA_real_),
  stringsAsFactors = FALSE
)
formatted <- sp_format_table(raw, digits = 2, p_digits = 3,
                             missing = "missing")
stopifnot(inherits(formatted, "spmodels_formatted_table"))
stopifnot(identical(formatted$estimate, c("1.23", "-0.00", "missing")))
stopifnot(identical(formatted$p_value, c("0.123", "<0.001", "missing")))
stopifnot(identical(attr(formatted, "raw_table"), raw))
captured_table <- capture.output(print(formatted))
stopifnot(any(grepl("alpha", captured_table)), any(grepl("<0.001", captured_table)))
expect_error(sp_format_table(raw, digits = -1))
expect_error(sp_format_table(raw, missing = NA_character_))
expect_error(sp_format_table(list()))

dat <- data.frame(
  y = c(1.2, 0.7, 2.1, 1.8, 3.0, 2.6),
  x = c(-1, 0, 0.5, 1, 1.5, 2),
  east = c(0, 1, 0, 2, 1, 3), north = c(0, 0, 1, 1, 2, 2)
)
control <- sp_control(
  smoothness = 0.5, spatial_variance = 1.4, scale = 1.7,
  nugget_variance = 0.3, estimate = FALSE, distance_units = "km"
)
fit <- sp_fit(y ~ x, dat, coords = ~ east + north, control = control)
newdata <- data.frame(x = c(0.25, 1.75), east = c(0.5, 2.5),
                      north = c(0.25, 1.5), row.names = c("west", "east"))
draws <- sp_predict_draws(fit, newdata, nsim = 20, seed = 920)
validation <- sp_cross_validate(fit, folds = c(1, 1, 2, 2, 3, 3))
report <- sp_report(fit, prediction = draws, validation = validation,
                    title = "Checked spatial analysis", digits = 4)
stopifnot(inherits(report, "spmodels_report"))
stopifnot(identical(report$tables$model, sp_model_table(fit)))
stopifnot(identical(report$tables$prediction, sp_prediction_table(draws)))
stopifnot(identical(report$tables$validation,
                    sp_validation_table(validation)))
stopifnot(identical(attr(report$formatted$model, "raw_table"),
                    report$tables$model))
stopifnot(report$metadata$family == "gaussian",
          report$metadata$observations == nrow(dat),
          report$metadata$distance_units == "km")
captured_report <- capture.output(print(report))
stopifnot(any(grepl("Checked spatial analysis", captured_report)),
          any(grepl("Model estimates", captured_report)),
          any(grepl("Predictions", captured_report)),
          any(grepl("Validation", captured_report)))

fold_report <- sp_report(fit, validation = validation,
                         validation_detail = "fold")
stopifnot(identical(fold_report$tables$validation,
                    sp_validation_table(validation, detail = "fold")))
expect_error(sp_report(list()))
expect_error(sp_report(fit, title = ""))
expect_error(sp_report(fit, prediction = list()))

cat("Formatted and assembled report tests passed.\n")
