if (file.exists(file.path("R", "intervention.R"))) {
  for (file in list.files("R", full.names = TRUE)) source(file)
} else {
  library(SPmodels)
}

expect_error <- function(expr, pattern = NULL) {
  value <- try(force(expr), silent = TRUE)
  stopifnot(inherits(value, "try-error"))
  if (!is.null(pattern)) stopifnot(grepl(pattern, as.character(value)))
}

history <- data.frame(
  area = c("A", "A", "B", "B"),
  year = c(1, 3, 2, 5),
  delivered = c(1, 1, 1, 0)
)
intervention <- sp_intervention(
  history, unit = "area", time = "year", alpha = 0.6, decay = 2,
  time_units = "years",
  treated = "delivered"
)
survey_rows <- data.frame(
  area = c("A", "A", "A", "B", "B", "C"),
  year = c(0, 1, 2, 2, 6, 4)
)
expected_multiplier <- c(
  1, 1, 1 - .6 * exp(-1 / 2), 1,
  1 - .6 * exp(-4 / 2), 1
)
stopifnot(all.equal(sp_intervention_multiplier(intervention, survey_rows),
                    expected_multiplier, tolerance = 1e-14),
          identical(intervention$time_rule, "strictly before survey"))

# The intervention-adjusted latent score and curvature equal centred finite
# differences of the independently evaluated binomial log posterior.
y_mode <- c(1, 0)
trials_mode <- c(3, 2)
linear_mode <- c(-.3, .2)
covariance_mode <- matrix(c(1, .2, .2, .8), 2)
multiplier_mode <- c(.7, .9)
mode <- sp_latent_mode(
  y_mode, linear_mode, covariance_mode, family = binomial(),
  trials = trials_mode, probability_multiplier = multiplier_mode
)
precision_mode <- solve(covariance_mode)
log_posterior <- function(s) {
  probability <- multiplier_mode * plogis(linear_mode + s)
  sum(dbinom(y_mode, trials_mode, probability, log = TRUE)) -
    .5 * drop(crossprod(s, precision_mode %*% s))
}
h <- 1e-5
numeric_gradient <- vapply(seq_along(mode$mode), function(j) {
  plus <- minus <- mode$mode
  plus[j] <- plus[j] + h
  minus[j] <- minus[j] - h
  (log_posterior(plus) - log_posterior(minus)) / (2 * h)
}, numeric(1))
numeric_hessian <- matrix(0, 2, 2)
for (j in 1:2) for (k in 1:2) {
  if (j == k) {
    plus <- minus <- mode$mode
    plus[j] <- plus[j] + h
    minus[j] <- minus[j] - h
    numeric_hessian[j, k] <-
      (log_posterior(plus) - 2 * log_posterior(mode$mode) +
         log_posterior(minus)) / h^2
  } else {
    pp <- pm <- mp <- mm <- mode$mode
    pp[j] <- pp[j] + h; pp[k] <- pp[k] + h
    pm[j] <- pm[j] + h; pm[k] <- pm[k] - h
    mp[j] <- mp[j] - h; mp[k] <- mp[k] + h
    mm[j] <- mm[j] - h; mm[k] <- mm[k] - h
    numeric_hessian[j, k] <-
      (log_posterior(pp) - log_posterior(pm) - log_posterior(mp) +
         log_posterior(mm)) / (4 * h^2)
  }
}
stopifnot(max(abs(numeric_gradient - mode$gradient)) < 1e-7,
          all.equal(mode$covariance, solve(-numeric_hessian),
                    tolerance = 2e-5))

dat <- data.frame(
  positive = c(0, 1, 0, 2, 1, 2, 1, 3),
  examined = c(4, 4, 5, 5, 4, 5, 5, 6),
  x = seq(-1, 1, length.out = 8),
  area = rep(c("A", "B"), each = 4),
  year = rep(2:5, 2),
  east = rep(0:3, 2), north = rep(0:1, each = 4)
)
control <- sp_control(
  spatial_variance = .3, scale = 1.5, nugget_variance = .7,
  estimate = FALSE, distance_units = "km"
)
fit <- sp_fit(
  cbind(positive, examined - positive) ~ x, dat,
  family = binomial(), coords = ~ east + north,
  intervention = intervention, control = control
)
training_multiplier <- sp_intervention_multiplier(intervention, dat)
stopifnot(identical(fit$intervention, intervention),
          all.equal(fit$intervention_multiplier, training_multiplier,
                    tolerance = 0),
          all.equal(fitted(fit), fit$trials * training_multiplier *
                      plogis(fit$latent_mode$linear_predictor),
                    tolerance = 1e-12))

new <- data.frame(
  x = c(-.4, .7), area = c("A", "B"), year = c(4, 6),
  east = c(.5, 2.5), north = c(.25, .8)
)
new_multiplier <- sp_intervention_multiplier(intervention, new)
link_prediction <- predict(fit, new, se.fit = TRUE, type = "link")
response_prediction <- predict(fit, new, se.fit = TRUE, type = "response")
moments <- vapply(seq_len(nrow(new)), function(i) {
  mu <- link_prediction$fit[i]
  standard_deviation <- link_prediction$se.fit[i]
  c(
    first = integrate(function(z) plogis(mu + standard_deviation * z) *
      dnorm(z), -Inf, Inf, rel.tol = 1e-9)$value,
    second = integrate(function(z) plogis(mu + standard_deviation * z)^2 *
      dnorm(z), -Inf, Inf, rel.tol = 1e-9)$value
  )
}, numeric(2))
expected_response <- new_multiplier * moments["first", ]
stopifnot(all.equal(response_prediction$fit, expected_response,
                    tolerance = 1e-10))

link_draws <- sp_predict_draws(fit, new, nsim = 20, seed = 710,
                               target = "latent", type = "link")
response_draws <- sp_predict_draws(fit, new, nsim = 20, seed = 710,
                                   target = "latent", type = "response")
stopifnot(all.equal(response_draws$draws,
                    sweep(plogis(link_draws$draws), 2, new_multiplier, "*"),
                    tolerance = 1e-14))

# Fitted simulation uses the same intervention multiplier before the binomial
# response draw.
set.seed(711)
z <- matrix(rnorm(nrow(dat) * 2), nrow(dat))
latent <- crossprod(fit$fitted_state$chol_sigma, z)
eta <- sweep(latent, 1, as.vector(fit$x %*% coef(fit) + fit$offset), "+")
probability <- sweep(plogis(eta), 1, training_multiplier, "*")
expected_simulation <- matrix(rbinom(
  length(probability), rep(fit$trials, 2), as.vector(probability)
), nrow(dat))
simulation <- simulate(fit, nsim = 2, seed = 711)
stopifnot(identical(unname(as.matrix(simulation)), expected_simulation))

model_table <- sp_model_table(fit)
stopifnot(identical(model_table$term[model_table$component == "intervention"],
                    c("alpha", "decay")),
          all(model_table$status[model_table$component == "intervention"] ==
                "fixed"))

# MCML weights and fold-wise refits retain the fixed history adjustment.
mc_fit <- suppressWarnings(sp_mcml(
  fit, control = sp_mcml_control(nsim = 100, maxit = 20,
                                 ess_warning = .01), seed = 712
))
stopifnot(identical(mc_fit$intervention, intervention),
          all(is.finite(mc_fit$mcml$log_weights)))
cv <- sp_cross_validate(
  fit, folds = c(1, 2, 1, 2, 1, 2, 1, 2)
)
stopifnot(nrow(cv$predictions) == nrow(dat),
          all(is.finite(cv$predictions$predicted)))

# Constrained Laplace estimation uses logit(alpha), log(decay), a normalized
# beta penalty, and the delta method for ordinary-scale uncertainty.
estimated_intervention <- sp_intervention(
  history, "area", "year", alpha = .5, decay = 2, time_units = "years",
  estimate = TRUE, alpha_penalty = c(2, 2), treated = "delivered"
)
estimated_fit <- suppressWarnings(sp_fit(
  cbind(positive, examined - positive) ~ x, dat,
  family = binomial(), coords = ~ east + north,
  intervention = estimated_intervention,
  control = sp_control(
    spatial_variance = .3, scale = 1.5, nugget_variance = .7,
    estimate = FALSE, distance_units = "km", maxit = 150
  )
))
stopifnot(estimated_fit$intervention$alpha > 0,
          estimated_fit$intervention$alpha < 1,
          estimated_fit$intervention$decay > 0,
          identical(dim(estimated_fit$intervention_vcov), c(2L, 2L)),
          identical(estimated_fit$likelihood_method, "Penalized Laplace"),
          attr(logLik(estimated_fit), "df") == length(coef(estimated_fit)) + 2L,
          all(sp_model_table(estimated_fit)$status[
            sp_model_table(estimated_fit)$component == "intervention"] ==
              "estimated"),
          all.equal(
            as.numeric(logLik(estimated_fit)),
            unname(estimated_fit$fitted_state$log_likelihood + dbeta(
              estimated_fit$intervention$alpha, 2, 2, log = TRUE
            )), tolerance = 1e-10
          ))
curve <- sp_intervention_curve(estimated_fit, lag = c(0, 1, 3), level = .9)
expected_reduction <- estimated_fit$intervention$alpha *
  exp(-curve$values$lag / estimated_fit$intervention$decay)
curve_gradient <- cbind(
  exp(-curve$values$lag / estimated_fit$intervention$decay),
  expected_reduction * curve$values$lag /
    estimated_fit$intervention$decay^2
)
expected_curve_se <- sqrt(pmax(0, rowSums(
  (curve_gradient %*% estimated_fit$intervention_vcov) * curve_gradient
)))
stopifnot(all.equal(curve$values$reduction, expected_reduction,
                    tolerance = 1e-13),
          all.equal(curve$values$multiplier, 1 - expected_reduction,
                    tolerance = 1e-13),
          all.equal(curve$values$std_error, expected_curve_se,
                    tolerance = 1e-12),
          all(curve$values$lower_reduction >= 0),
          all(curve$values$upper_reduction <= 1))
fixed_curve <- sp_intervention_curve(intervention, lag = 0:2)
stopifnot(all(is.na(fixed_curve$values$std_error)))
plot_file <- tempfile(fileext = ".pdf")
grDevices::pdf(plot_file)
plot(curve)
grDevices::dev.off()
stopifnot(file.info(plot_file)$size > 0)
unlink(plot_file)
estimated_mcml <- suppressWarnings(sp_mcml(
  estimated_fit,
  control = sp_mcml_control(nsim = 200, maxit = 60,
                            ess_warning = .001), seed = 713
))
maximum_log_weight <- max(estimated_mcml$mcml$log_weights)
mcml_log_mean_weight <- maximum_log_weight + log(mean(exp(
  estimated_mcml$mcml$log_weights - maximum_log_weight
)))
stopifnot(estimated_mcml$optimisation$convergence == 0L,
          estimated_mcml$intervention$alpha > 0,
          estimated_mcml$intervention$alpha < 1,
          estimated_mcml$intervention$decay > 0,
          identical(dim(estimated_mcml$intervention_vcov), c(2L, 2L)),
          all(is.finite(estimated_mcml$mcml$log_weights)),
          all.equal(as.numeric(logLik(estimated_mcml)),
                    mcml_log_mean_weight,
                    tolerance = 1e-12),
          all.equal(estimated_mcml$intervention_multiplier,
                    sp_intervention_multiplier(
                      estimated_mcml$intervention, dat
                    ), tolerance = 1e-14))
mc_samples <- estimated_mcml$mcml$latent_samples
mc_delta <- sweep(mc_samples, 2L, estimated_mcml$mcml$proposal_mean, "-")
mc_proposal_chol <- chol(estimated_mcml$mcml$proposal_covariance)
mc_standard <- t(forwardsolve(t(mc_proposal_chol), t(mc_delta)))
mc_log_q <- -.5 * (ncol(mc_samples) * log(2 * pi) +
  2 * sum(log(diag(mc_proposal_chol))) + rowSums(mc_standard^2))
mc_eta <- sweep(mc_samples, 2L,
  as.vector(estimated_mcml$x %*% coef(estimated_mcml) +
              estimated_mcml$offset), "+")
mc_probability <- sweep(plogis(mc_eta), 2L,
                         estimated_mcml$intervention_multiplier, "*")
mc_log_response <- rowSums(dbinom(
  matrix(estimated_mcml$y, nrow(mc_samples), ncol(mc_samples), byrow = TRUE),
  matrix(estimated_mcml$trials, nrow(mc_samples), ncol(mc_samples), byrow = TRUE),
  mc_probability, log = TRUE
)) + dbeta(estimated_mcml$intervention$alpha, 2, 2, log = TRUE)
mc_sigma_chol <- chol(estimated_mcml$fitted_state$sigma)
mc_precision_samples <- backsolve(
  mc_sigma_chol, forwardsolve(t(mc_sigma_chol), t(mc_samples))
)
mc_log_prior <- -.5 * (ncol(mc_samples) * log(2 * pi) +
  2 * sum(log(diag(mc_sigma_chol))) +
  rowSums(mc_samples * t(mc_precision_samples)))
stopifnot(all.equal(estimated_mcml$mcml$log_weights,
  mc_log_response + mc_log_prior - mc_log_q, tolerance = 1e-10))
estimated_parameter_draws <- sp_predict_draws(
  estimated_fit, new, nsim = 6, target = "latent", type = "response",
  parameter_uncertainty = TRUE, seed = 714
)
estimated_neighbor_draws <- sp_predict_draws_neighbor(
  estimated_fit, new, nsim = 6, neighbors = nrow(new),
  target = "latent", type = "response",
  parameter_uncertainty = TRUE, seed = 714
)
multiplier_expected <- t(vapply(seq_len(6), function(i) {
  current <- estimated_fit$intervention
  current$alpha <- plogis(estimated_parameter_draws$parameter_draws[i, 3])
  current$decay <- exp(estimated_parameter_draws$parameter_draws[i, 4])
  sp_intervention_multiplier(current, new)
}, numeric(nrow(new))))
stopifnot(identical(dim(estimated_parameter_draws$parameter_draws), c(6L, 4L)),
          all.equal(estimated_parameter_draws$intervention_multiplier_draws,
                    multiplier_expected, tolerance = 1e-14),
          all.equal(estimated_neighbor_draws$parameter_draws,
                    estimated_parameter_draws$parameter_draws,
                    tolerance = 0),
          all.equal(estimated_neighbor_draws$intervention_multiplier_draws,
                    multiplier_expected, tolerance = 1e-14),
          all.equal(estimated_neighbor_draws$draws,
                    estimated_parameter_draws$draws, tolerance = 1e-8),
          all(estimated_parameter_draws$draws >= 0),
          all(estimated_parameter_draws$draws <= multiplier_expected))

expect_error(sp_intervention(
  rbind(history, history[1, ]), "area", "year", .5, 1, "years"
), "at most one")
expect_error(sp_intervention(history, "area", "year", 1.1, 1, "years"),
             "between")
expect_error(sp_intervention(history, "area", "year", .5, 1, ""),
             "time_units")
expect_error(sp_intervention(history, "area", "year", 0, 1, "years",
  estimate = TRUE), "strictly between")
expect_error(sp_intervention(history, "area", "year", .5, 1, "years",
  estimate = TRUE, alpha_penalty = c(1, 0)), "positive")
expect_error(sp_intervention_multiplier(intervention,
  data.frame(area = "A")), "require columns")
expect_error(sp_fit(positive ~ x, dat, family = poisson(),
  coords = ~ east + north, intervention = intervention, control = control),
  "require a binomial")
expect_error(sp_latent_mode(0, 0, matrix(1), family = poisson(),
  probability_multiplier = 1), "only used for binomial")

cat("Decay-adjusted intervention-history tests passed.\n")
