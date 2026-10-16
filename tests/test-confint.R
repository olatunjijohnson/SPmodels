if (file.exists(file.path("R", "confint.R"))) {
  for (file in list.files("R", full.names = TRUE)) source(file)
} else {
  library(SPmodels)
}

close <- function(a, b, tolerance = 1e-11) {
  isTRUE(all.equal(a, b, tolerance = tolerance, check.attributes = FALSE))
}
expect_error <- function(expression, pattern) {
  error <- tryCatch({ force(expression); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}

dat <- data.frame(
  y = c(1.2, 0.7, 2.1, 1.8, 3.0, 2.6),
  x = c(-1, 0, 0.5, 1, 1.5, 2)
)
fit <- sp_fit(y ~ x, dat, spatial = FALSE)
interval <- confint(fit, level = 0.90)
z <- qnorm(0.95)
beta <- coef(fit)
beta_se <- sqrt(diag(vcov(fit)))
sigma2 <- fit$covariance[["nugget_variance"]]
sigma2_se <- sqrt(fit$covariance_vcov[
  "nugget_variance", "nugget_variance"
])
stopifnot(identical(row.names(interval),
                    c("mean::(Intercept)", "mean::x",
                      "covariance::nugget_variance")),
          close(interval[1:2, 1], beta - z * beta_se),
          close(interval[1:2, 2], beta + z * beta_se),
          close(interval[3, ], exp(log(sigma2) +
            c(-1, 1) * z * sigma2_se / sigma2)),
          identical(attr(interval, "transform"),
                    c("identity", "identity", "log")))

selected <- confint(fit, parm = c("covariance::nugget_variance", "mean::x"))
stopifnot(identical(row.names(selected),
                    c("covariance::nugget_variance", "mean::x")),
          identical(attr(selected, "transform"), c("log", "identity")))
expect_error(confint(fit, parm = "scale"), "Unknown")
expect_error(confint(fit, parm = 9), "valid")
expect_error(confint(fit, level = 1), "strictly between")

# Fixed covariance parameters are not mislabeled as estimated intervals.
spatial_dat <- transform(dat, east = c(0, 1, 0, 2, 1, 3),
                         north = c(0, 0, 1, 1, 2, 2))
fixed <- sp_fit(y ~ x, spatial_dat, coords = ~ east + north,
  control = sp_control(spatial_variance = 1, scale = 1.5,
    nugget_variance = 0.2, estimate = FALSE, distance_units = "km"))
stopifnot(identical(row.names(confint(fixed)),
                    c("mean::(Intercept)", "mean::x")))

# Probability and positive constrained components use logit and log scales.
constrained <- structure(list(
  coefficients = c("(Intercept)" = 0), coefficient_vcov = matrix(0.04, 1, 1),
  spatial = TRUE,
  covariance = c(spatial_variance = 1, scale = 2, nugget_variance = 0.2),
  covariance_vcov = matrix(NA_real_, 3, 3,
    dimnames = list(c("spatial_variance", "scale", "nugget_variance"),
                    c("spatial_variance", "scale", "nugget_variance"))),
  covariance_estimated = FALSE, smoothness = 0.5, distance_units = "km",
  random_variance = numeric(), intervention = list(
    alpha = 0.4, decay = 3, time_units = "years"
  ), intervention_estimated = TRUE,
  intervention_vcov = structure(diag(c(0.05^2, 0.4^2), 2),
    dimnames = list(c("alpha", "decay"), c("alpha", "decay"))),
  temporal = list(time_units = "years"), temporal_scale = 2,
  temporal_scale_estimated = TRUE, temporal_scale_se = 0.3
), class = "spmodels_fit")
bounded <- confint(constrained)
stopifnot(close(bounded["intervention::alpha", ],
  plogis(qlogis(0.4) + c(-1, 1) * qnorm(0.975) * 0.05 / (0.4 * 0.6))),
  all(bounded["intervention::alpha", ] > 0 &
      bounded["intervention::alpha", ] < 1),
  close(bounded["intervention::decay", ],
        exp(log(3) + c(-1, 1) * qnorm(0.975) * 0.4 / 3)),
  close(bounded["temporal::scale", ],
        exp(log(2) + c(-1, 1) * qnorm(0.975) * 0.3 / 2)))

cat("Fitted-parameter confidence interval tests passed.\n")
