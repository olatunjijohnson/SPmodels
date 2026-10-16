if (file.exists(file.path("R", "prediction-update.R"))) {
  for (file in list.files("R", full.names = TRUE)) source(file)
} else {
  library(SPmodels)
}

expect_error <- function(expr) {
  value <- try(force(expr), silent = TRUE)
  stopifnot(inherits(value, "try-error"))
}

dat <- data.frame(
  y = c(1.1, 1.7, 2.2, 2.9, 3.1, 3.8),
  x = c(-1, 0, 0.5, 1, 1.5, 2), exposure = c(1, 2, 1, 3, 2, 4),
  group = factor(c("low", "low", "high", "high", "low", "high")),
  east = c(0, 1, 0, 2, 1, 3), north = c(0, 0, 1, 1, 2, 2)
)
control <- sp_control(
  smoothness = 0.5, spatial_variance = 0.8, scale = 1.4,
  nugget_variance = 0.2, estimate = FALSE, distance_units = "km"
)
fit <- sp_fit(y ~ x + group + offset(log(exposure)), dat,
              coords = ~ east + north, control = control)
newdata <- data.frame(
  x = c(0.2, 0.7, 1.2, 1.7), exposure = c(2, 2, 3, 3),
  group = factor(c("low", "high", "low", "high"),
                 levels = levels(dat$group)),
  east = c(0.3, 0.8, 1.3, 1.8), north = c(0.2, 0.7, 1.2, 1.7),
  retained = LETTERS[1:4], row.names = c("site_a", "site_b", "site_c", "site_d")
)
coefficient_before <- coef(fit)
covariance_before <- fit$covariance
scenario <- sp_prediction_update(
  fit, newdata,
  changes = list(
    x = c(site_d = 2.4, site_b = 1.4),
    exposure = c(site_d = 6, site_b = 4)
  ),
  rows = c("site_b", "site_d"), se.fit = TRUE,
  target = "observation"
)

manual <- newdata
manual[c("site_b", "site_d"), "x"] <- c(1.4, 2.4)
manual[c("site_b", "site_d"), "exposure"] <- c(4, 6)
baseline_direct <- predict(fit, newdata, se.fit = TRUE,
                           target = "observation")
updated_direct <- predict(fit, manual, se.fit = TRUE,
                          target = "observation")
stopifnot(identical(scenario$data, manual))
stopifnot(identical(scenario$comparison$row, row.names(newdata)))
stopifnot(identical(scenario$comparison$changed,
                    c(FALSE, TRUE, FALSE, TRUE)))
stopifnot(all.equal(scenario$comparison$baseline, baseline_direct$fit,
                    tolerance = 1e-13))
stopifnot(all.equal(scenario$comparison$updated, updated_direct$fit,
                    tolerance = 1e-13))
stopifnot(all.equal(scenario$comparison$difference,
                    updated_direct$fit - baseline_direct$fit,
                    tolerance = 1e-13))
stopifnot(all.equal(scenario$comparison$baseline_se, baseline_direct$se.fit,
                    tolerance = 1e-13))
stopifnot(all.equal(scenario$comparison$updated_se, updated_direct$se.fit,
                    tolerance = 1e-13))
stopifnot(all.equal(scenario$baseline_offset, log(newdata$exposure),
                    tolerance = 1e-14))
stopifnot(all.equal(scenario$updated_offset, log(manual$exposure),
                    tolerance = 1e-14))
stopifnot(identical(scenario$baseline_design[c(1, 3), , drop = FALSE],
                    scenario$updated_design[c(1, 3), , drop = FALSE]))
stopifnot(identical(newdata$retained, scenario$data$retained))
stopifnot(all(scenario$preservation))
stopifnot(identical(coef(fit), coefficient_before),
          identical(fit$covariance, covariance_before))

# Scalar and logical-row updates preserve exact row order and factor levels.
factor_scenario <- sp_prediction_update(
  fit, newdata, changes = list(group = "high"),
  rows = c(TRUE, FALSE, TRUE, FALSE)
)
stopifnot(identical(as.character(factor_scenario$data$group),
                    c("high", "high", "high", "high")))
stopifnot(identical(levels(factor_scenario$data$group), levels(newdata$group)))

expect_error(sp_prediction_update(fit, newdata, list(missing = 1)))
expect_error(sp_prediction_update(fit, newdata, list(east = 10)))
expect_error(sp_prediction_update(fit, newdata, list(x = 1), rows = integer()))
expect_error(sp_prediction_update(fit, newdata, list(x = c(a = 1, b = 2)),
                                  rows = c("site_b", "site_d")))
expect_error(sp_prediction_update(fit, newdata, list(group = "new")))

external <- sp_fit(y ~ x, dat, coords = ~ east + north,
                   offset = log(dat$exposure), control = control)
expect_error(sp_prediction_update(external, newdata, list(x = 1)))

if (requireNamespace("sf", quietly = TRUE)) {
  sf_data <- sf::st_as_sf(dat, coords = c("east", "north"), crs = 32630)
  sf_fit <- sp_fit(y ~ x + offset(log(exposure)), sf_data, control = control)
  sf_new <- sf::st_as_sf(newdata, coords = c("east", "north"), crs = 32630)
  sf_scenario <- sp_prediction_update(sf_fit, sf_new,
                                      changes = list(exposure = 5), rows = 2L)
  stopifnot(inherits(sf_scenario$data, "sf"))
  stopifnot(identical(sf::st_geometry(sf_new),
                      sf::st_geometry(sf_scenario$data)))
  stopifnot(isTRUE(sf_scenario$preservation[["geometry"]]))
}

cat("Prediction scenario update tests passed.\n")
