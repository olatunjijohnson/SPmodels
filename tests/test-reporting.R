if (file.exists(file.path("R", "fit.R"))) {
  for (file in list.files("R", full.names = TRUE)) source(file)
} else {
  library(SPmodels)
}

expect_error <- function(expression, pattern) {
  error <- tryCatch({ force(expression); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}

dat <- data.frame(
  response = c(1.2, 0.7, 2.1, 1.8, 3.0, 2.6),
  exposure = c(-1, 0, 0.5, 1, 1.5, 2),
  east = c(0, 1, 0, 2, 1, 3), north = c(0, 0, 1, 1, 2, 2)
)
control <- sp_control(
  smoothness = 0.5, spatial_variance = 1.4, scale = 1.7,
  nugget_variance = 0.3, estimate = FALSE, distance_units = "km"
)
fit <- sp_fit(response ~ exposure, dat, coords = ~ east + north,
              control = control)

model_table <- sp_model_table(fit)
stopifnot(identical(names(model_table), c(
  "component", "term", "estimate", "std_error", "statistic", "p_value",
  "status", "units"
)))
stopifnot(identical(model_table$term,
  c("(Intercept)", "exposure", "spatial_variance", "scale",
    "nugget_variance", "smoothness")))
stopifnot(all.equal(model_table$estimate[1:2], unname(coef(fit)),
                    tolerance = 1e-12))
stopifnot(model_table$units[model_table$term == "scale"] == "km")
stopifnot(all(model_table$status[model_table$component == "covariance"] ==
                "fixed"))

new <- data.frame(exposure = c(0.25, 1.75), east = c(0.5, 2.5),
                  north = c(0.25, 1.5), row.names = c("west", "east"))
draws <- sp_predict_draws(fit, new, nsim = 40, seed = 812)
prediction_table <- sp_prediction_table(draws, threshold = 2,
                                         probs = c(0.1, 0.9))
stopifnot(identical(prediction_table$target, c("west", "east")))
stopifnot(identical(names(prediction_table), c(
  "target", "mean", "sd", "q_0.1", "q_0.9",
  "exceedance_probability", "family", "prediction_target", "scale"
)))
stopifnot(all.equal(prediction_table$mean, unname(colMeans(draws$draws)),
                    tolerance = 1e-12))
stopifnot(all.equal(prediction_table$exceedance_probability,
                    unname(colMeans(draws$draws > 2)), tolerance = 1e-12))

target <- sp_predict_target(draws, weights = c(1, 3), threshold = 2)
target_table <- sp_prediction_table(target)
stopifnot(nrow(target_table) == 1L, target_table$target == "all")
expect_error(sp_prediction_table(target, threshold = 2), "already defined")

folds <- c(1, 1, 2, 2, 3, 3)
cross_validation <- sp_cross_validate(fit, folds = folds)
overall_table <- sp_validation_table(cross_validation)
fold_table <- sp_validation_table(cross_validation, detail = "fold")
stopifnot(identical(names(overall_table),
  c("assessment", "scope", "id", "metric", "value")))
stopifnot(all(overall_table$assessment == "cross_validation"),
          all(overall_table$scope == "overall"),
          identical(unique(fold_table$id), c("1", "2", "3")))
stopifnot(all(c("rmse", "mae", "mean_log_score") %in%
                overall_table$metric))

assessment <- sp_simulation_assess(draws, truth = c(1.5, 2.5))
assessment_table <- sp_validation_table(assessment)
stopifnot(all(assessment_table$assessment == "predictive_assessment"),
          "coverage" %in% assessment_table$metric)
expect_error(sp_validation_table(assessment, detail = "fold"),
             "only available")

# The same application-facing schemas work for prevalence and count models.
non_gaussian <- data.frame(
  infected = c(0, 1, 1, 2, 2, 3, 1, 3),
  examined = c(2, 2, 3, 3, 3, 4, 2, 4),
  cases = c(0, 1, 1, 2, 2, 4, 3, 5),
  person_time = c(1, 2, 1, 3, 2, 4, 3, 5),
  exposure = c(-1, -0.6, -0.2, 0.1, 0.4, 0.8, 1.1, 1.5),
  east = c(0, 1, 0, 1, 2, 2, 3, 3),
  north = c(0, 0, 1, 1, 0, 1, 0, 1)
)
non_gaussian_control <- sp_control(
  spatial_variance = 0.7, scale = 1.4, nugget_variance = 0.2,
  estimate = FALSE, distance_units = "km"
)
prevalence_fit <- sp_fit(
  cbind(infected, examined - infected) ~ exposure, non_gaussian,
  family = binomial(), coords = ~ east + north,
  control = non_gaussian_control
)
count_fit <- sp_fit(
  cases ~ exposure + offset(log(person_time)), non_gaussian,
  family = poisson(), coords = ~ east + north,
  control = non_gaussian_control
)
prevalence_draws <- sp_predict_draws(prevalence_fit, nsim = 20, seed = 813,
                                     target = "latent", type = "response")
count_draws <- sp_predict_draws(count_fit, nsim = 20, seed = 814,
                                target = "latent", type = "response")
prevalence_table <- sp_prediction_table(prevalence_draws, threshold = 0.5)
count_table <- sp_prediction_table(count_draws, threshold = 2)
stopifnot(all(prevalence_draws$draws >= 0 & prevalence_draws$draws <= 1),
          all(count_draws$draws >= 0),
          all(prevalence_table$family == "binomial"),
          all(count_table$family == "poisson"),
          identical(names(prevalence_table), names(count_table)))

recovery <- sp_simulation_recovery(fit, nsim = 1, seed = 815)
recovery_table <- sp_validation_table(recovery)
stopifnot(all(recovery_table$assessment == "simulation_recovery"),
          all(recovery_table$scope == "parameter"),
          all(unique(recovery_table$id) %in% names(coef(fit))))
expect_error(sp_model_table(list()), "fitted")

cat("Structured reporting tests passed.\n")
