if (file.exists(file.path("R", "matern.R"))) {
  for (file in list.files("R", full.names = TRUE)) source(file)
} else if (!exists("matern_cor", mode = "function")) {
  library(SPmodels)
}

internal_matern <- if (exists(".sp_matern_cor_unchecked", mode = "function")) {
  .sp_matern_cor_unchecked
} else {
  getFromNamespace(".sp_matern_cor_unchecked", "SPmodels")
}

close_enough <- function(actual, expected, tolerance = 1e-12) {
  isTRUE(all.equal(actual, expected, tolerance = tolerance,
                   check.attributes = TRUE))
}

expect_error <- function(expr) {
  inherits(try(force(expr), silent = TRUE), "try-error")
}

# Equation M-01 and its zero-distance limit (M-02).
d <- c(zero = 0, half = 0.5, one = 1, two = 2)
stopifnot(close_enough(matern_cor(d, scale = 1, smoothness = 0.5), exp(-d)))

# Closed forms M-05 and M-06.
x <- c(0, 0.01, 0.5, 1, 5, 20)
expected_3_2 <- (1 + x) * exp(-x)
expected_5_2 <- (1 + x + x^2 / 3) * exp(-x)
stopifnot(close_enough(matern_cor(x, 1, 1.5), expected_3_2))
stopifnot(close_enough(matern_cor(x, 1, 2.5), expected_5_2))

# Only the distance-to-scale ratio matters (M-03).
stopifnot(close_enough(
  matern_cor(d, scale = 2, smoothness = 0.8),
  matern_cor(10 * d, scale = 20, smoothness = 0.8)
))

# A covariance matrix constructed from valid Euclidean distances is symmetric
# and positive semidefinite up to floating-point error.
locations <- c(-1.0, -0.1, 0.0, 0.7, 2.2)
distance_matrix <- abs(outer(locations, locations, "-"))
correlation_matrix <- matern_cor(distance_matrix, 0.75, 1.2)
stopifnot(isSymmetric(correlation_matrix))
stopifnot(all(diag(correlation_matrix) == 1))
stopifnot(min(eigen(correlation_matrix, symmetric = TRUE,
                    only.values = TRUE)$values) > -1e-12)

# Shape and labels are preserved.
labelled <- matrix(c(0, 1, 1, 0), nrow = 2,
                   dimnames = list(c("a", "b"), c("a", "b")))
result <- matern_cor(labelled, 1, 0.5)
stopifnot(identical(dim(result), dim(labelled)))
stopifnot(identical(dimnames(result), dimnames(labelled)))

# The optimizer-only path receives already validated inputs. It must remain
# numerically identical to the public implementation for every specialized
# half-integer case and for the generic Bessel route, including attributes.
for (nu in c(0.5, 1.5, 2.5, 0.8, 3.1)) {
  public <- matern_cor(labelled, scale = 0.7, smoothness = nu)
  internal <- internal_matern(labelled, scale = 0.7, smoothness = nu)
  stopifnot(close_enough(internal, public, tolerance = 1e-14))
  stopifnot(identical(dim(internal), dim(labelled)))
  stopifnot(identical(dimnames(internal), dimnames(labelled)))
}

# Contract failures are explicit.
stopifnot(expect_error(matern_cor(-1, 1, 0.5)))
stopifnot(expect_error(matern_cor(NA_real_, 1, 0.5)))
stopifnot(expect_error(matern_cor(Inf, 1, 0.5)))
stopifnot(expect_error(matern_cor(1, 0, 0.5)))
stopifnot(expect_error(matern_cor(1, 1, 0)))
stopifnot(expect_error(matern_cor("1", 1, 0.5)))

