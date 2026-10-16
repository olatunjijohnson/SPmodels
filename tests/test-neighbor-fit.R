if (file.exists(file.path("R", "neighbor-fit.R"))) {
  for (file in list.files("R", full.names = TRUE)) source(file)
} else {
  library(SPmodels)
}

set.seed(421)
n <- 12L
d <- data.frame(
  east = runif(n), north = runif(n), z = rnorm(n), exposure = runif(n, 1, 2)
)
distance <- as.matrix(stats::dist(as.matrix(d[c("east", "north")])))
truth <- c(spatial_variance = 1.1, scale = 0.6, nugget_variance = 0.15)
covariance <- truth[["spatial_variance"]] * exp(-distance / truth[["scale"]])
diag(covariance) <- diag(covariance) + truth[["nugget_variance"]]
d$y <- as.vector(0.5 - 0.35 * d$z + t(chol(covariance)) %*% rnorm(n))
control <- sp_control(
  smoothness = 0.5, spatial_variance = truth[["spatial_variance"]],
  scale = truth[["scale"]], nugget_variance = truth[["nugget_variance"]],
  estimate = FALSE, distance_units = "unit"
)

dense <- sp_fit(y ~ z, d, coords = ~ east + north, control = control)
neighbor <- sp_fit_neighbor(
  y ~ z, d, coords = ~ east + north, neighbors = n - 1L,
  control = control
)
stopifnot(
  inherits(neighbor, "spmodels_neighbor_fit"),
  identical(neighbor$likelihood_method, "Gaussian nearest-predecessor ML"),
  isTRUE(all.equal(as.numeric(logLik(neighbor)), dense$log_likelihood,
                   tolerance = 2e-10)),
  isTRUE(all.equal(coef(neighbor), dense$coefficients, tolerance = 2e-10)),
  isTRUE(all.equal(vcov(neighbor), dense$coefficient_vcov, tolerance = 2e-10)),
  length(neighbor$neighbor_graph) == n,
  identical(lengths(neighbor$neighbor_graph), 0:(n - 1L)),
  all(neighbor$fitted_state$conditional_variance > 0),
  length(fitted(neighbor)) == n,
  length(residuals(neighbor)) == n
)
fixed_table <- sp_model_table(neighbor)
stopifnot(
  identical(fixed_table$component,
            c("mean", "mean", rep("covariance", 4L))),
  all(fixed_table$status[fixed_table$component == "covariance"] == "fixed")
)

maximin <- sp_fit_neighbor(
  y ~ z, d, coords = ~ east + north, neighbors = n - 1L,
  ordering = "maximin", control = control
)
stopifnot(
  identical(maximin$observation_order[[1L]], 1L),
  isTRUE(all.equal(as.numeric(logLik(maximin)), dense$log_likelihood,
                   tolerance = 2e-10)),
  isTRUE(all.equal(coef(maximin), dense$coefficients, tolerance = 2e-10))
)

newdata <- d[1:3, ]
newdata$east <- newdata$east + 0.02
dense_prediction <- predict(dense, newdata, se.fit = TRUE,
                            target = "observation")
neighbor_prediction <- predict(neighbor, newdata, se.fit = TRUE,
                               target = "observation")
stopifnot(
  isTRUE(all.equal(neighbor_prediction$fit, dense_prediction$fit,
                   tolerance = 2e-10)),
  isTRUE(all.equal(neighbor_prediction$se.fit, dense_prediction$se.fit,
                   tolerance = 2e-10))
)

explicit_order <- rev(seq_len(n))
explicit <- sp_fit_neighbor(
  y ~ z, d, coords = ~ east + north, neighbors = 3,
  observation_order = explicit_order, control = control
)
stopifnot(
  identical(explicit$ordering, "user"),
  identical(explicit$observation_order, explicit_order),
  max(lengths(explicit$neighbor_graph)) == 3L,
  is.finite(as.numeric(logLik(explicit)))
)

prediction_data <- d[1:4, ]
prediction_data$east <- prediction_data$east + 0.03
prediction <- predict(explicit, prediction_data, se.fit = TRUE,
                      target = "observation")
new_distance <- sqrt(rowSums(sweep(
  as.matrix(d[c("east", "north")]), 2L,
  unname(unlist(prediction_data[1L, c("east", "north")])), "-"
)^2))
selected <- order(new_distance, seq_along(new_distance))[1:3]
neighbor_covariance <- truth[["spatial_variance"]] *
  matern_cor(as.matrix(stats::dist(as.matrix(
    d[selected, c("east", "north")]
  ))), truth[["scale"]], 0.5)
diag(neighbor_covariance) <- diag(neighbor_covariance) +
  truth[["nugget_variance"]]
cross <- truth[["spatial_variance"]] *
  matern_cor(new_distance[selected], truth[["scale"]], 0.5)
local_coefficient <- solve(neighbor_covariance, cross)
x_new <- model.matrix(~ z, prediction_data)[1L, ]
manual_mean <- sum(x_new * coef(explicit)) +
  sum(local_coefficient * explicit$residual_mean[selected])
mean_gradient_new <- x_new - as.vector(crossprod(
  local_coefficient, explicit$x[selected, , drop = FALSE]
))
manual_variance <- truth[["spatial_variance"]] +
  truth[["nugget_variance"]] - sum(cross * local_coefficient) +
  sum((mean_gradient_new %*% vcov(explicit)) * mean_gradient_new)
stopifnot(
  length(prediction$fit) == 4L,
  length(prediction$se.fit) == 4L,
  all(is.finite(prediction$fit)),
  all(is.finite(prediction$se.fit)),
  all(prediction$se.fit > 0),
  identical(prediction$approximation, "nearest-predecessor"),
  isTRUE(all.equal(prediction$fit[[1L]], manual_mean, tolerance = 1e-12)),
  isTRUE(all.equal(prediction$se.fit[[1L]]^2, manual_variance,
                   tolerance = 1e-12))
)
training_latent <- predict(explicit, se.fit = TRUE, target = "latent")
training_observation <- predict(explicit, se.fit = TRUE,
                                target = "observation")
mean_variance <- rowSums(
  (explicit$fitted_state$transformed_x %*% vcov(explicit)) *
    explicit$fitted_state$transformed_x
)
stopifnot(
  all(training_observation$se.fit >= training_latent$se.fit),
  isTRUE(all.equal(training_latent$se.fit^2,
    explicit$fitted_state$conditional_variance -
      explicit$covariance[["nugget_variance"]] + mean_variance,
    tolerance = 1e-12)),
  isTRUE(all.equal(training_observation$se.fit^2 -
    training_latent$se.fit^2,
    rep(explicit$covariance[["nugget_variance"]], n), tolerance = 1e-12))
)

offset_fit <- sp_fit_neighbor(
  y ~ z + offset(log(exposure)), d, coords = ~ east + north,
  neighbors = 3, control = control
)
stopifnot(length(predict(offset_fit, d[1:2, ])) == 2L)

set.seed(99)
seed_before <- .Random.seed
invisible(sp_fit_neighbor(y ~ z, d, coords = ~ east + north,
                          neighbors = 2, control = control))
stopifnot(identical(seed_before, .Random.seed))

estimated_control <- sp_control(
  smoothness = 0.5, spatial_variance = 1, scale = 0.5,
  nugget_variance = 0.2, estimate = TRUE, distance_units = "unit",
  maxit = 100
)
estimated <- suppressWarnings(sp_fit_neighbor(
  y ~ z, d, coords = ~ east + north, neighbors = 3,
  control = estimated_control
))
stopifnot(
  all(is.finite(estimated$covariance)),
  all(estimated$covariance > 0),
  is.finite(as.numeric(logLik(estimated))),
  identical(dim(estimated$covariance_vcov), c(3L, 3L))
)
estimated_table <- sp_model_table(estimated)
interval <- confint(estimated, parm = c("mean::(Intercept)",
                                        "covariance::spatial_variance"))
stopifnot(
  all(estimated_table$status[estimated_table$term %in%
    names(estimated$covariance)] == "estimated"),
  identical(dim(interval), c(2L, 2L)),
  identical(attr(interval, "transform"), c("identity", "log"))
)

invisible(capture.output(print(explicit)))
invisible(capture.output(print(summary(explicit))))

expect_error <- function(expression, pattern) {
  error <- tryCatch({ force(expression); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}
expect_error(sp_fit_neighbor(y ~ z, d, ~ east + north, neighbors = 0,
                             control = control), "positive integer")
expect_error(sp_fit_neighbor(y ~ z, d, ~ east + north, neighbors = 2,
  observation_order = rep(1, n), control = control), "permutation")
expect_error(sp_fit_neighbor(y ~ z, transform(d, east = 1, north = 1),
  ~ east + north, neighbors = 2, control = control), "distinct")
group_control <- control
group_control$random_variance <- 1
expect_error(sp_fit_neighbor(y ~ z, d, ~ east + north, neighbors = 2,
                             control = group_control), "grouped covariance")

message("Gaussian nearest-predecessor fitting tests passed.")
