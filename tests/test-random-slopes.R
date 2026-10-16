if (file.exists(file.path("R", "random-effects.R"))) {
  for (file in list.files("R", full.names = TRUE)) source(file)
} else {
  library(SPmodels)
}

expect_error <- function(expression, pattern) {
  error <- tryCatch({ force(expression); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}

internal <- function(name) {
  if (exists(name, mode = "function", inherits = TRUE)) {
    get(name, mode = "function", inherits = TRUE)
  } else {
    getFromNamespace(name, "SPmodels")
  }
}
random_slope_training <- internal(".sp_random_slope_training")
random_slope_parameterization <- internal(
  ".sp_random_slope_parameterization"
)
random_slope_update <- internal(".sp_random_slope_update")
random_slope_parameters <- internal(".sp_random_slope_parameters")
random_slope_jacobian <- internal(".sp_random_slope_natural_jacobian")
estimated_slope_parameters <- internal(
  ".sp_random_slope_estimated_parameters"
)
parameter_specification <- internal(".sp_parameter_specification")
object_at_parameters <- internal(".sp_object_at_parameters")
prediction_design <- internal(".sp_prediction_design")
prediction_random_slope <- internal(".sp_prediction_random_slope")

dat <- data.frame(
  y = c(0.8, 1.3, 1.1, 2.0, 2.4, 2.8, 3.1, 3.7),
  x = c(-1, -.4, .1, .7, -.8, -.1, .6, 1.2),
  clinic = rep(c("A", "B"), each = 4L),
  east = c(0, .5, 1, 1.5, 2, 2.5, 3, 3.5),
  north = c(0, .3, .8, 1.2, 1.7, 2, 2.8, 3.1)
)
group_covariance <- matrix(c(.45, .12, .12, .25), 2L)
grouped <- sp_random_slope("clinic", "x", group_covariance)
pars <- c(spatial_variance = .6, scale = 1.4, nugget_variance = .15)
control <- sp_control(
  estimate = FALSE, smoothness = .5,
  spatial_variance = pars[[1L]], scale = pars[[2L]],
  nugget_variance = pars[[3L]], distance_units = "km"
)
fit <- sp_fit(y ~ x, dat, coords = ~ east + north,
              random_slope = grouped, control = control)

# Independent covariance, GLS, and likelihood calculation.
coordinates <- as.matrix(dat[c("east", "north")])
distance <- as.matrix(dist(coordinates))
w <- cbind(1, dat$x)
same_group <- outer(dat$clinic, dat$clinic, `==`)
group_part <- same_group * (w %*% group_covariance %*% t(w))
v <- pars[["spatial_variance"]] * exp(-distance / pars[["scale"]]) +
  group_part + diag(pars[["nugget_variance"]], nrow(dat))
x <- model.matrix(y ~ x, dat)
v_inverse <- solve(v)
beta <- solve(crossprod(x, v_inverse %*% x),
              crossprod(x, v_inverse %*% dat$y))
residual <- dat$y - as.vector(x %*% beta)
expected_log_likelihood <- -.5 * (nrow(dat) * log(2 * pi) +
  as.numeric(determinant(v, logarithm = TRUE)$modulus) +
  drop(crossprod(residual, v_inverse %*% residual)))
stopifnot(all.equal(unname(coef(fit)), as.vector(beta), tolerance = 1e-11))
stopifnot(all.equal(fit$fitted_state$v, v, tolerance = 1e-12))
stopifnot(all.equal(as.numeric(logLik(fit)), expected_log_likelihood,
                    tolerance = 1e-11))
stopifnot(identical(attr(logLik(fit), "df"), 2L))

# Known and new-level prediction covariance follows the declared group law.
new <- data.frame(x = c(.2, .9, -.3), clinic = c("A", "C", "C"),
                  east = c(.7, 4, 4.4), north = c(.4, 3.4, 3.7))
expect_error(predict(fit, new), "unseen levels")
known <- new[1L, , drop = FALSE]
cross_distance <- sqrt(rowSums((coordinates -
  matrix(as.numeric(known[c("east", "north")]), nrow(dat), 2L,
         byrow = TRUE))^2))
cross_group <- as.numeric(dat$clinic == "A") *
  as.vector(w %*% group_covariance %*% c(1, known$x))
cross_covariance <- pars[["spatial_variance"]] *
  exp(-cross_distance / pars[["scale"]]) + cross_group
expected_prediction <- as.vector(c(1, known$x) %*% beta +
  crossprod(cross_covariance, v_inverse %*% residual))
stopifnot(all.equal(predict(fit, known), expected_prediction,
                    tolerance = 1e-11))

draws <- sp_predict_draws(fit, new, nsim = 10L, seed = 12L,
                          unseen = "marginal")
new_w <- cbind(1, new$x)
new_group_part <- outer(new$clinic, new$clinic, `==`) *
  (new_w %*% group_covariance %*% t(new_w))
new_coordinates <- as.matrix(new[c("east", "north")])
cross_distance <- sqrt(outer(rowSums(coordinates^2),
  rowSums(new_coordinates^2), "+") - 2 * tcrossprod(coordinates,
                                                     new_coordinates))
cross_group <- outer(dat$clinic, new$clinic, `==`) *
  (w %*% group_covariance %*% t(new_w))
cross_covariance <- pars[["spatial_variance"]] *
  exp(-cross_distance / pars[["scale"]]) + cross_group
new_covariance <- pars[["spatial_variance"]] *
  exp(-as.matrix(dist(new_coordinates)) / pars[["scale"]]) + new_group_part
x_new <- model.matrix(~ x, new)
beta_vcov <- solve(crossprod(x, v_inverse %*% x))
mean_difference <- x_new - crossprod(cross_covariance, v_inverse %*% x)
expected_joint_covariance <- new_covariance -
  crossprod(cross_covariance, v_inverse %*% cross_covariance) +
  mean_difference %*% beta_vcov %*% t(mean_difference)
stopifnot(all.equal(draws$predictive_covariance, expected_joint_covariance,
                    tolerance = 1e-11))
stopifnot(new_group_part[2L, 3L] != 0,
          all(cross_group[, 2L:3L] == 0))
stopifnot(identical(dim(draws$draws), c(10L, 3L)))

# Laplace and MCML states use the same fixed grouped covariance.
binary <- dat
binary$y <- c(0, 0, 0, 1, 0, 1, 1, 1)
binary_fit <- sp_fit(y ~ x, binary, family = binomial(),
  coords = ~ east + north, random_slope = grouped, control = control)
expected_sigma <- pars[["spatial_variance"]] *
  exp(-distance / pars[["scale"]]) + group_part +
  diag(pars[["nugget_variance"]], nrow(dat))
stopifnot(all.equal(binary_fit$fitted_state$sigma, expected_sigma,
                    tolerance = 1e-12))
mcml <- sp_mcml(binary_fit,
  control = sp_mcml_control(nsim = 100L, maxit = 2L), seed = 31L)
stopifnot(all.equal(mcml$fitted_state$sigma, expected_sigma,
                    tolerance = 1e-12))

# Fixed covariance parameters are reported, but are not counted as estimated.
table <- sp_model_table(fit)
fixed_rows <- table$component == "random_slope"
stopifnot(sum(fixed_rows) == 3L, all(table$status[fixed_rows] == "fixed"),
          all(is.na(table$std_error[fixed_rows])))
stopifnot(length(summary(fit)$random_slope) == 3L)

# Additive components and a slope-only component retain their own group laws.
dat$survey <- rep(c("early", "late"), 4L)
slope_only <- sp_random_slope(
  "survey", "x", matrix(.09, 1L), intercept = FALSE,
  label = "survey-slope"
)
additive <- sp_fit(y ~ x, dat, coords = ~ east + north,
  random_slope = list(grouped, slope_only), control = control)
survey_part <- outer(dat$survey, dat$survey, `==`) *
  (.09 * outer(dat$x, dat$x))
stopifnot(all.equal(additive$fitted_state$v, v + survey_part,
                    tolerance = 1e-12),
          length(summary(additive)$random_slope) == 4L)

# Whole-group validation requires an explicit population-level prediction for
# levels deliberately absent from each training fit.
clinic_folds <- sp_group_folds(fit, group = "clinic", folds = 2L)
expect_error(sp_cross_validate(fit, clinic_folds), "unseen levels")
clinic_validation <- sp_cross_validate(
  fit, clinic_folds, unseen = "marginal"
)
stopifnot(all(is.finite(clinic_validation$predictions$predicted)),
          all(is.finite(clinic_validation$predictions$se.fit)),
          identical(clinic_validation$unseen, "marginal"))
binary_new <- transform(binary[1:2, ], clinic = "new-clinic")
binary_score <- sp_predictive_log_score(
  binary_fit, binary_new, observed = binary_new$y, unseen = "marginal"
)
stopifnot(all(is.finite(binary_score)))

# Constructor and data contracts fail explicitly.
expect_error(sp_random_slope("clinic", "x", diag(c(1, 0))),
             "positive definite")
expect_error(sp_random_slope("clinic", "x", matrix(c(1, .2, .1, 1), 2L)),
             "symmetric")
expect_error(sp_random_slope("clinic", "x", 1), "matrix")
expect_error(sp_fit(y ~ x, dat, coords = ~ east + north,
  random_slope = list(grouped, grouped), control = control), "unique")
expect_error(sp_fit(y ~ x, transform(dat, x = as.character(x)),
  coords = ~ east + north, random_slope = grouped, control = control),
  "finite numeric")

# Estimated covariance matrices use log-Cholesky diagonal parameters, so every
# optimizer evaluation and fitted result is positive definite. The natural-
# scale covariance and its delta-method uncertainty are reported explicitly.
estimated_data <- expand.grid(
  z = c(-1, 0, 1), group = paste0("g", seq_len(8L)),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
estimated_data$east <- seq_len(nrow(estimated_data))
estimated_data$north <- rep(0:2, 8L)
true_group_covariance <- matrix(c(.7, .18, .18, .25), 2L)
set.seed(601)
group_coefficients <- matrix(rnorm(16L), 8L, 2L) %*%
  chol(true_group_covariance)
group_index <- match(estimated_data$group, paste0("g", seq_len(8L)))
estimated_data$y <- 1 + .4 * estimated_data$z + rowSums(
  cbind(1, estimated_data$z) * group_coefficients[group_index, ]
) + rnorm(nrow(estimated_data), sd = .3)
estimated_specification <- sp_random_slope(
  "group", "z", true_group_covariance, estimate = TRUE
)
estimated_fit <- sp_fit(
  y ~ z, estimated_data, coords = ~ east + north,
  random_slope = estimated_specification,
  control = sp_control(distance_units = "km", maxit = 80L)
)
estimated_covariance <- estimated_fit$random_slope[[1L]]$covariance
stopifnot(isTRUE(estimated_fit$random_slope[[1L]]$estimate),
          min(eigen(estimated_covariance, symmetric = TRUE,
                    only.values = TRUE)$values) > 0,
          identical(dim(estimated_fit$random_slope_vcov), c(3L, 3L)),
          all(is.finite(estimated_fit$random_slope_vcov)),
          attr(logLik(estimated_fit), "df") == 8L)
estimated_table <- sp_model_table(estimated_fit)
estimated_rows <- estimated_table$component == "random_slope"
stopifnot(all(estimated_table$status[estimated_rows] == "estimated"),
          all(is.finite(estimated_table$std_error[estimated_rows])))

# A rebuilt parameter state propagates its Cholesky covariance into the
# prediction-side random-slope block, rather than retaining the covariance at
# the fitted centre. A small deterministic perturbation isolates this map from
# the separate question of whether an unconstrained Hessian draw is usable.
estimated_parameter_specification <- parameter_specification(estimated_fit)
estimated_design <- prediction_design(
  estimated_fit, estimated_data[1:4, ], "error"
)
perturbed_parameters <- estimated_parameter_specification$values
slope_start <- estimated_parameter_specification$beta_count +
  estimated_parameter_specification$covariance_count +
  estimated_parameter_specification$random_count
perturbed_parameters[[slope_start + 1L]] <-
  perturbed_parameters[[slope_start + 1L]] + .05
draw_object <- object_at_parameters(
  estimated_fit, perturbed_parameters, estimated_parameter_specification
)
prediction_slope <- prediction_random_slope(draw_object, estimated_design)
stopifnot(!isTRUE(all.equal(
            prediction_slope[[1L]]$covariance, estimated_covariance,
            tolerance = 0
          )),
          isTRUE(all.equal(
            prediction_slope[[1L]]$covariance,
            draw_object$random_slope_specification[[1L]]$covariance,
            tolerance = 0
          )))
estimated_recovery <- suppressWarnings(sp_simulation_recovery(
  estimated_fit, nsim = 1L, seed = 604L
))
stopifnot(sum(grepl("^random_slope::", names(estimated_recovery$truth))) == 3L,
          sum(grepl("^random_slope::",
                    colnames(estimated_recovery$estimates))) == 3L)
slope_truth <- estimated_slope_parameters(
  estimated_fit$random_slope_specification
)
grid_values <- as.list(unname(slope_truth))
names(grid_values) <- paste0("random_slope::", names(slope_truth))
slope_grid <- sp_simulation_grid(estimated_fit, grid_values,
                                 prefix = "slope")
stopifnot(length(slope_grid) == 1L,
          isTRUE(all.equal(
            slope_grid[[1L]]$random_slope[[1L]]$covariance,
            estimated_fit$random_slope[[1L]]$covariance,
            tolerance = 1e-14
          )))
covariance_name <- paste0("random_slope::",
  grep("::cov\\(", names(slope_truth), value = TRUE))
expect_error(sp_simulation_grid(
  estimated_fit, stats::setNames(list(1e6), covariance_name)
), "positive-definite")
estimated_interval <- confint(estimated_fit)
covariance_row <- grep("random_slope::.*::cov\\(",
                       row.names(estimated_interval))
stopifnot(length(covariance_row) == 1L,
          attr(estimated_interval, "transform")[covariance_row] == "identity")

# Direct reconstruction verifies the parameter map and its analytic Jacobian.
training_specification <- random_slope_training(
  list(estimated_specification), estimated_data
)
parameterization <- random_slope_parameterization(training_specification)
reconstructed <- random_slope_update(
  training_specification, parameterization$values, parameterization
)
stopifnot(isTRUE(all.equal(unname(reconstructed[[1L]]$covariance),
                          true_group_covariance, tolerance = 1e-14)))
jacobian <- random_slope_jacobian(reconstructed)
step <- 1e-6
numeric_jacobian <- sapply(seq_along(parameterization$values), function(j) {
  upper <- lower <- parameterization$values
  upper[j] <- upper[j] + step
  lower[j] <- lower[j] - step
  (random_slope_parameters(random_slope_update(
    training_specification, upper, parameterization
  )) - random_slope_parameters(random_slope_update(
    training_specification, lower, parameterization
  ))) / (2 * step)
})
stopifnot(isTRUE(all.equal(unname(jacobian), unname(numeric_jacobian),
                          tolerance = 1e-8)))

# The same constrained state is available to a non-Gaussian Laplace fit and
# its MCML correction.
count_data <- expand.grid(
  z = c(-1, 0, 1, 2), group = paste0("g", seq_len(8L)),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
count_data$east <- seq_len(nrow(count_data))
count_data$north <- rep(0:3, 8L)
count_covariance <- matrix(c(.15, .03, .03, .08), 2L)
set.seed(602)
count_coefficients <- matrix(rnorm(16L), 8L, 2L) %*%
  chol(count_covariance)
count_index <- match(count_data$group, paste0("g", seq_len(8L)))
count_eta <- .3 + .15 * count_data$z + rowSums(
  cbind(1, count_data$z) * count_coefficients[count_index, ]
)
count_data$y <- rpois(nrow(count_data), exp(count_eta))
count_fit <- suppressWarnings(sp_fit(
  y ~ z, count_data, family = poisson(), coords = ~ east + north,
  random_slope = sp_random_slope(
    "group", "z", count_covariance, estimate = TRUE
  ),
  control = sp_control(spatial_variance = .2, scale = 3,
    nugget_variance = .1, distance_units = "km", maxit = 30L)
))
count_mcml <- suppressWarnings(sp_mcml(
  count_fit,
  control = sp_mcml_control(nsim = 100L, maxit = 2L,
                            ess_warning = .001), seed = 603L
))
count_repeats <- suppressWarnings(sp_mcml_replicates(
  count_fit, nrep = 2L,
  control = sp_mcml_control(nsim = 100L, maxit = 1L,
                            ess_warning = .001), seed = 605L
))
stopifnot(min(eigen(count_fit$random_slope[[1L]]$covariance,
                    symmetric = TRUE, only.values = TRUE)$values) > 0,
          min(eigen(count_mcml$random_slope[[1L]]$covariance,
                    symmetric = TRUE, only.values = TRUE)$values) > 0,
          nrow(count_mcml$parameter_vcov) == 8L,
          sum(grepl("^random_slope::",
                    count_repeats$natural_summary$term)) == 3L,
          identical(dim(count_repeats$natural_sampling_covariance),
                    c(8L, 8L)))

# Parameter indexing remains aligned when grouped covariance and temporal
# scale are estimated together before and during MCML.
count_data$year <- rep(0:3, 8L) + rep(0:7, each = 4L) * .2
combined_fit <- suppressWarnings(sp_fit(
  y ~ z, count_data, family = poisson(), coords = ~ east + north,
  time = sp_temporal(~ year, 1, "years", estimate = TRUE),
  random_slope = sp_random_slope(
    "group", "z", count_covariance, estimate = TRUE
  ),
  control = sp_control(spatial_variance = .2, scale = 3,
    nugget_variance = .1, distance_units = "km", maxit = 30L)
))
combined_mcml <- suppressWarnings(sp_mcml(
  combined_fit,
  control = sp_mcml_control(nsim = 100L, maxit = 1L,
                            ess_warning = .001), seed = 606L
))
stopifnot(length(combined_fit$optimisation$par) == 9L,
          attr(logLik(combined_fit), "df") == 9L,
          length(combined_mcml$optimisation$par) == 9L,
          combined_mcml$temporal_scale > 0,
          min(eigen(combined_mcml$random_slope[[1L]]$covariance,
                    symmetric = TRUE, only.values = TRUE)$values) > 0)

# Identifiability checks reject designs that cannot distinguish the grouped
# coefficient covariance from another component.
expect_error(sp_random_slope("group", "z", true_group_covariance,
                             estimate = NA), "TRUE or FALSE")
expect_error(sp_fit(y ~ z, estimated_data, coords = ~ east + north,
  random_slope = estimated_specification,
  control = sp_control(spatial_variance = 1, scale = 1,
    nugget_variance = 1, estimate = FALSE, distance_units = "km")),
  "Estimated random-slope covariance requires")
expect_error(sp_fit(y ~ z, estimated_data, coords = ~ east + north,
  random = ~ group, random_slope = estimated_specification,
  control = sp_control(distance_units = "km")), "same grouping variable")
expect_error(sp_fit(y ~ z, estimated_data, coords = ~ east + north,
  random_slope = list(estimated_specification, sp_random_slope(
    "group", "z", true_group_covariance, label = "second", estimate = TRUE
  )), control = sp_control(distance_units = "km")), "Only one estimated")
too_few <- estimated_data[estimated_data$group %in% c("g1", "g2"), ]
expect_error(sp_fit(y ~ z, too_few, coords = ~ east + north,
  random_slope = estimated_specification,
  control = sp_control(distance_units = "km")), "at least 3 observed groups")
no_within_variation <- transform(
  estimated_data, z = match(group, unique(group))
)
expect_error(sp_fit(y ~ 1, no_within_variation, coords = ~ east + north,
  random_slope = estimated_specification,
  control = sp_control(distance_units = "km")), "distinct slope values")
one_per_group <- estimated_data[estimated_data$z == -1, ]
expect_error(sp_fit(y ~ 1, one_per_group, coords = ~ east + north,
  random_slope = sp_random_slope(
    "group", "z", matrix(.2, 1L), intercept = FALSE, estimate = TRUE
  ), control = sp_control(distance_units = "km")),
  "two non-zero slope observations")

cat("random-slope tests passed\n")
