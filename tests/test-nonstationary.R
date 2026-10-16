if (file.exists(file.path("R", "nonstationary.R"))) {
  for (file in list.files("R", full.names = TRUE)) source(file)
} else {
  library(SPmodels)
}

expect_error <- function(expression, pattern) {
  message <- tryCatch({
    force(expression)
    NULL
  }, error = function(error) conditionMessage(error))
  if (is.null(message)) stop("Expected an error, but none was thrown.")
  if (!grepl(pattern, message, fixed = TRUE)) {
    stop(sprintf("Expected error matching '%s'; got '%s'.", pattern, message))
  }
  invisible(NULL)
}

coords <- cbind(x = c(0, 0.7, 1.4, 2), y = c(0, 0.2, -0.1, 0.5))
scale <- 1.3
local <- array(0, dim = c(nrow(coords), 2, 2))
for (i in seq_len(nrow(coords))) local[i, , ] <- diag(scale^2, 2)
actual <- sp_nonstationary_cor(coords, local, smoothness = 1.5)
expected <- matern_cor(as.matrix(stats::dist(coords)), scale, 1.5)
dimnames(actual) <- dimnames(expected) <- NULL
stopifnot(isTRUE(all.equal(actual, expected, tolerance = 1e-12)))

coords <- cbind(x = c(0, 1, 0.5, 2, 1.5), y = c(0, 0, 1, 1, 2))
local <- array(0, dim = c(nrow(coords), 2, 2))
angles <- c(0, pi / 6, pi / 3, pi / 2, pi / 4)
scales <- c(0.6, 1, 1.4, 0.8, 1.2)
for (i in seq_len(nrow(coords))) {
  rotation <- matrix(c(cos(angles[i]), -sin(angles[i]),
                       sin(angles[i]), cos(angles[i])), 2, 2)
  local[i, , ] <- rotation %*% diag(c(scales[i]^2, 1.5 * scales[i]^2)) %*%
    t(rotation)
}
correlation <- sp_nonstationary_cor(coords, local, smoothness = 0.5)
stopifnot(isTRUE(all.equal(correlation, t(correlation), tolerance = 1e-14)),
          isTRUE(all.equal(diag(correlation), rep(1, nrow(coords)),
                           tolerance = 0)),
          min(eigen(correlation, symmetric = TRUE, only.values = TRUE)$values) >=
            -1e-10,
          correlation[1, 2] < correlation[4, 5])

coords <- rbind(c(0, 0), c(1, 0))
good <- array(0, dim = c(2, 2, 2))
good[1, , ] <- diag(2)
good[2, , ] <- diag(2)
expect_error(sp_nonstationary_cor(c(0, 1), good), "coords")
expect_error(sp_nonstationary_cor(coords, matrix(1, 2, 2)), "n-by-2-by-2")
expect_error(sp_nonstationary_cor(coords, good, smoothness = 0), "smoothness")
asymmetric <- good
asymmetric[1, 1, 2] <- 0.3
expect_error(sp_nonstationary_cor(coords, asymmetric), "symmetric")
indefinite <- good
indefinite[2, , ] <- matrix(c(1, 2, 2, 1), 2)
expect_error(sp_nonstationary_cor(coords, indefinite), "positive definite")
