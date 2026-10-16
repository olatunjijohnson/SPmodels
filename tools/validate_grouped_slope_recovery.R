args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1L) args[[1L]] else
  "docs/validation/data/grouped_slope_recovery.csv"
output_path <- if (length(args) >= 2L) args[[2L]] else input_path
for (file in list.files("R", full.names = TRUE)) source(file)

data <- utils::read.csv(input_path, stringsAsFactors = FALSE,
                        check.names = FALSE)
cholesky_error <- rep(NA_real_, nrow(data))
jacobian_error <- rep(NA_real_, nrow(data))
for (row in seq_len(nrow(data))) {
  values <- c(
    data[["estimate__group:z::var((Intercept))"]][[row]],
    data[["estimate__group:z::cov(z,(Intercept))"]][[row]],
    data[["estimate__group:z::var(z)"]][[row]]
  )
  if (!all(is.finite(values))) next
  covariance <- matrix(c(values[[1L]], values[[2L]],
                         values[[2L]], values[[3L]]), 2L,
                       dimnames = list(c("(Intercept)", "z"),
                                      c("(Intercept)", "z")))
  factor <- tryCatch(chol(covariance), error = function(error) NULL)
  if (is.null(factor)) next
  cholesky_error[[row]] <- max(abs(t(factor) %*% factor - covariance))
  specification <- list(sp_random_slope(
    "group", "z", covariance, estimate = TRUE
  ))
  parameterization <- .sp_random_slope_parameterization(specification)
  analytic <- .sp_random_slope_natural_jacobian(specification)
  numeric <- analytic
  for (column in seq_along(parameterization$values)) {
    step <- 1e-6 * max(1, abs(parameterization$values[[column]]))
    plus <- parameterization$values
    minus <- parameterization$values
    plus[[column]] <- plus[[column]] + step
    minus[[column]] <- minus[[column]] - step
    natural_plus <- .sp_random_slope_parameters(.sp_random_slope_update(
      specification, plus, parameterization
    ))
    natural_minus <- .sp_random_slope_parameters(.sp_random_slope_update(
      specification, minus, parameterization
    ))
    numeric[, column] <- (natural_plus - natural_minus) / (2 * step)
  }
  jacobian_error[[row]] <- max(abs(analytic - numeric))
}
data$cholesky_reconstruction_error <- cholesky_error
data$jacobian_max_abs_error <- jacobian_error
data$mathematical_checks_pass <- is.finite(cholesky_error) &
  cholesky_error <= 1e-12 & is.finite(jacobian_error) &
  jacobian_error <= 1e-6
data$usable <- data$usable & data$mathematical_checks_pass
utils::write.csv(data, output_path, row.names = FALSE, na = "")
cat("Cholesky reconstruction passed:", sum(cholesky_error <= 1e-12,
                                              na.rm = TRUE), "/", nrow(data),
    "\n")
cat("Analytic Jacobian passed:", sum(jacobian_error <= 1e-6, na.rm = TRUE),
    "/", nrow(data), "\n")
cat("Mathematical checks passed:", sum(data$mathematical_checks_pass,
                                          na.rm = TRUE), "/", nrow(data),
    "\n")
