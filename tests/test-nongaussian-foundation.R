if (file.exists(file.path("R", "latent-mode.R"))) {
  for (file in list.files("R", full.names = TRUE)) source(file)
} else {
  library(SPmodels)
}

close <- function(a, b, tolerance = 1e-7) {
  isTRUE(all.equal(a, b, tolerance = tolerance, check.attributes = FALSE))
}

# At smoothness one half, correlation is exponential and both derivatives have
# transparent closed forms.
d <- c(0, 0.2, 1, 3, 8)
rho <- 2
correlation <- exp(-d / rho)
first_expected <- correlation * d / rho^2
second_expected <- correlation * (d^2 / rho^4 - 2 * d / rho^3)
stopifnot(close(matern_scale_derivative(d, rho, 0.5, 1), first_expected,
                1e-12))
stopifnot(close(matern_scale_derivative(d, rho, 0.5, 2), second_expected,
                1e-12))

# General-smoothness derivatives agree with centred finite differences.
d_general <- c(0.1, 0.7, 2.5)
rho_general <- 1.4
nu_general <- 1.2
h1 <- 1e-5
first_fd <- (matern_cor(d_general, rho_general + h1, nu_general) -
               matern_cor(d_general, rho_general - h1, nu_general)) / (2 * h1)
h2 <- 1e-4
second_fd <- (matern_cor(d_general, rho_general + h2, nu_general) -
                2 * matern_cor(d_general, rho_general, nu_general) +
                matern_cor(d_general, rho_general - h2, nu_general)) / h2^2
stopifnot(close(matern_scale_derivative(d_general, rho_general, nu_general, 1),
                first_fd, 2e-8))
stopifnot(close(matern_scale_derivative(d_general, rho_general, nu_general, 2),
                second_fd, 2e-6))

# Binomial mode: direct score is zero and local covariance inverts curvature.
distance <- abs(outer(1:5, 1:5, "-"))
sigma <- 1.3 * exp(-distance / 1.7) + diag(0.2, 5)
y_binomial <- c(0, 2, 1, 3, 1)
trials <- c(2, 3, 2, 4, 2)
fixed <- c(-0.4, 0.2, 0.1, 0.7, -0.1)
binomial_mode <- sp_latent_mode(y_binomial, fixed, sigma,
                                family = binomial(), trials = trials)
stopifnot(binomial_mode$converged)
precision <- solve(sigma)
probability <- plogis(fixed + binomial_mode$mode)
score <- y_binomial - trials * probability -
  as.vector(precision %*% binomial_mode$mode)
curvature <- precision + diag(trials * probability * (1 - probability))
stopifnot(max(abs(score)) < 1e-9)
stopifnot(close(binomial_mode$covariance, solve(curvature), 1e-10))

# Poisson mode and row-permutation invariance.
y_poisson <- c(0, 3, 1, 5, 2)
poisson_mode <- sp_latent_mode(y_poisson, fixed, sigma, family = poisson())
stopifnot(poisson_mode$converged, max(abs(poisson_mode$gradient)) < 1e-9)
permutation <- c(5, 2, 4, 1, 3)
poisson_permuted <- sp_latent_mode(
  y_poisson[permutation], fixed[permutation],
  sigma[permutation, permutation], family = poisson()
)
stopifnot(close(poisson_permuted$mode,
                poisson_mode$mode[permutation], 1e-10))

expect_error <- function(expression, pattern) {
  error <- tryCatch({ force(expression); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}
expect_error(matern_scale_derivative(1, 1, 0.5, 3), "order")
expect_error(sp_latent_mode(c(0, 2), c(0, 0), diag(2),
                            family = binomial()), "y.*trials")
expect_error(sp_latent_mode(c(1, -1), c(0, 0), diag(2),
                            family = poisson()), "non-negative")
expect_error(sp_latent_mode(c(1, 1), c(0, 0), matrix(1, 2, 2),
                            family = poisson()), "positive definite")

# One-dimensional MALA moments agree with numerical integration within a
# predeclared stochastic margin, and a supplied seed is local/reproducible.
control <- sp_mcmc_control(iterations = 12000L, burnin = 2000L,
                           thin = 5L, step_size = 0.9)
set.seed(808)
seed_before <- .Random.seed
sample_one <- sp_latent_sample(1, 0, matrix(1, 1, 1), family = binomial(),
                               control = control, seed = 42)
stopifnot(identical(.Random.seed, seed_before))
sample_two <- sp_latent_sample(1, 0, matrix(1, 1, 1), family = binomial(),
                               control = control, seed = 42)
stopifnot(identical(sample_one$samples, sample_two$samples))
unnormalized <- function(s) {
  exp(dbinom(1, 1, plogis(s), log = TRUE) + dnorm(s, log = TRUE))
}
normalizer <- integrate(unnormalized, -Inf, Inf, rel.tol = 1e-11)$value
exact_mean <- integrate(function(s) s * unnormalized(s), -Inf, Inf,
                        rel.tol = 1e-11)$value / normalizer
exact_second <- integrate(function(s) s^2 * unnormalized(s), -Inf, Inf,
                          rel.tol = 1e-11)$value / normalizer
exact_sd <- sqrt(exact_second - exact_mean^2)
stopifnot(abs(mean(sample_one$samples) - exact_mean) < 0.04)
stopifnot(abs(sd(sample_one$samples) - exact_sd) < 0.05)
stopifnot(sample_one$acceptance_rate > 0.2,
          sample_one$acceptance_rate < 0.99)
sampler_diagnostics <- sp_sampler_diagnostics(sample_one, max_lag = 50)
stopifnot(inherits(sampler_diagnostics, "spmodels_sampler_diagnostics"))
stopifnot(nrow(sampler_diagnostics$summary) == 1L)
stopifnot(sampler_diagnostics$summary$ess > 1,
          sampler_diagnostics$summary$ess <= nrow(sample_one$samples))
stopifnot(length(sampler_diagnostics$autocorrelation) == 50L)
expect_error(sp_mcmc_control(iterations = 10, burnin = 10), "smaller")

# Multiple independently evolved chains retain their identity, support
# overdispersed starts, and produce rank-normalized split R-hat plus an
# autocorrelation-aware Monte Carlo standard error.
multi_control <- sp_mcmc_control(
  iterations = 6000L, burnin = 1000L, thin = 5L, step_size = 0.9,
  chains = 4L, initial_scale = 2
)
set.seed(809)
multi_seed_before <- .Random.seed
multi_sample <- sp_latent_sample(
  1, 0, matrix(1, 1, 1), family = binomial(),
  control = multi_control, seed = 91
)
stopifnot(identical(.Random.seed, multi_seed_before),
          identical(dim(multi_sample$chain_samples), c(1000L, 1L, 4L)),
          nrow(multi_sample$samples) == 4000L,
          identical(as.integer(table(multi_sample$chain_id)), rep(1000L, 4)),
          length(multi_sample$acceptance_rates) == 4L,
          all(multi_sample$acceptance_rates > 0.2 &
                multi_sample$acceptance_rates < 0.99))
multi_repeat <- sp_latent_sample(
  1, 0, matrix(1, 1, 1), family = binomial(),
  control = multi_control, seed = 91
)
stopifnot(identical(multi_sample$chain_samples,
                    multi_repeat$chain_samples))
multi_diagnostics <- sp_sampler_diagnostics(multi_sample, max_lag = 50)
stopifnot(multi_diagnostics$chains == 4L,
          multi_diagnostics$retained_per_chain == 1000L,
          isTRUE(multi_diagnostics$converged),
          is.finite(multi_diagnostics$summary$rhat),
          multi_diagnostics$summary$rhat < 1.01,
          multi_diagnostics$summary$ess > 1,
          multi_diagnostics$summary$ess <= 4000,
          close(multi_diagnostics$summary$mcse,
                multi_diagnostics$summary$sd /
                  sqrt(multi_diagnostics$summary$ess), 1e-12),
          identical(dim(multi_diagnostics$trace), c(1000L, 4L)))
selected_diagnostics <- sp_sampler_diagnostics(
  multi_sample, max_lag = 25, chains = c("chain_1", "chain_3")
)
stopifnot(selected_diagnostics$chains == 2L,
          identical(selected_diagnostics$selected_chains,
                    c("chain_1", "chain_3")),
          identical(dim(selected_diagnostics$trace), c(1000L, 2L)))

# Independently reconstruct rank-normalized split/folded R-hat for fixed
# chains, including the split order and Blom normal-score transformation.
known_chains <- cbind(
  c(-1.2, -0.7, -0.1, 0.4, -0.8, -0.2, 0.5, 1.0),
  c(-1.0, -0.5, 0.0, 0.6, -0.6, 0.1, 0.7, 1.2)
)
known_array <- array(known_chains, dim = c(8L, 1L, 2L),
  dimnames = list(NULL, "s1", c("chain_1", "chain_2")))
known_object <- structure(list(
  samples = rbind(known_chains[, 1, drop = FALSE],
                  known_chains[, 2, drop = FALSE]),
  chain_samples = known_array,
  acceptance_rate = 0.8,
  acceptance_rates = c(chain_1 = 0.75, chain_2 = 0.85)
), class = "spmodels_latent_sample")
known_diagnostics <- sp_sampler_diagnostics(known_object, max_lag = 3)
split <- cbind(known_chains[1:4, , drop = FALSE],
               known_chains[5:8, , drop = FALSE])
manual_normalize <- function(values) {
  ranks <- rank(values, ties.method = "average")
  qnorm((ranks - 3 / 8) / (length(values) + 1 / 4))
}
manual_rhat <- function(values) {
  within <- mean(apply(values, 2, var))
  between <- nrow(values) * var(colMeans(values))
  sqrt(((nrow(values) - 1) / nrow(values) * within +
          between / nrow(values)) / within)
}
ranked <- matrix(manual_normalize(as.vector(split)), nrow = 4)
folded <- abs(split - median(as.vector(split)))
folded_ranked <- matrix(manual_normalize(as.vector(folded)), nrow = 4)
known_rhat_expected <- max(manual_rhat(ranked), manual_rhat(folded_ranked))
stopifnot(close(known_diagnostics$summary$rhat,
                known_rhat_expected, 1e-12))
expect_error(sp_mcmc_control(chains = 0), "chains")
expect_error(sp_mcmc_control(initial_scale = -1), "initial_scale")
expect_error(sp_sampler_diagnostics(multi_sample, chains = c(1, 1)),
             "distinct")

cat("Non-Gaussian foundation tests passed.\n")
