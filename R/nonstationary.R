#' Non-stationary Matérn correlation
#'
#' Construct a valid non-stationary Matérn correlation matrix from a local
#' positive-definite scale matrix at each location. The construction combines
#' local Mahalanobis distances with a determinant adjustment, so locally
#' varying anisotropy and range yield a valid covariance rather than an
#' arbitrary pairwise-distance kernel.
#'
#' Each local matrix has units of coordinate-distance squared. An isotropic
#' local matrix `rho^2 * diag(2)` therefore gives local range `rho`.
#'
#' @param coords A finite numeric matrix with two columns, one row per
#'   location, in common projected Euclidean units.
#' @param local_covariance A numeric `n x 2 x 2` array, where slice `i` is the
#'   symmetric positive-definite local scale matrix at `coords[i, ]`.
#' @param smoothness A finite positive Matérn smoothness parameter, held fixed
#'   across locations.
#' @examples
#' locations <- cbind(east = c(0, 1, 2), north = c(0, 0.5, 1))
#' local <- array(0, dim = c(nrow(locations), 2, 2))
#' for (i in seq_len(nrow(locations))) {
#'   local[i, , ] <- diag(c(1 + i / 4, 0.8 + i / 5)^2, 2)
#' }
#' sp_nonstationary_cor(locations, local, smoothness = 1.5)
#'
#' @return A symmetric `n x n` correlation matrix with unit diagonal. The
#'   determinant construction is valid for any supplied positive-definite
#'   local matrices; this function does not estimate them or fit a regression.
#'
#' @references Paciorek, C. J. and Schervish, M. J. (2006). Spatial modelling
#'   using a new class of nonstationary covariance functions.
#'   \doi{10.1002/env.785}.
#' @export
sp_nonstationary_cor <- function(coords, local_covariance, smoothness = 0.5) {
  if (!is.matrix(coords) || !is.numeric(coords) || ncol(coords) != 2L ||
      nrow(coords) < 1L || anyNA(coords) || any(!is.finite(coords))) {
    stop("`coords` must be a finite numeric matrix with two columns.",
         call. = FALSE)
  }
  n <- nrow(coords)
  if (!is.array(local_covariance) || !is.numeric(local_covariance) ||
      length(dim(local_covariance)) != 3L ||
      !identical(as.integer(dim(local_covariance)), c(as.integer(n), 2L, 2L)) ||
      anyNA(local_covariance) || any(!is.finite(local_covariance))) {
    stop("`local_covariance` must be a finite numeric n-by-2-by-2 array.",
         call. = FALSE)
  }
  if (!is.numeric(smoothness) || length(smoothness) != 1L ||
      is.na(smoothness) || !is.finite(smoothness) || smoothness <= 0) {
    stop("`smoothness` must be one finite positive number.", call. = FALSE)
  }

  local <- lapply(seq_len(n), function(i) {
    value <- matrix(local_covariance[i, , ], nrow = 2L, ncol = 2L)
    if (max(abs(value - t(value))) > 1e-10 * max(1, max(abs(value)))) {
      stop(sprintf("Local matrix %d must be symmetric.", i), call. = FALSE)
    }
    value <- (value + t(value)) / 2
    if (is.null(tryCatch(chol(value), error = function(e) NULL))) {
      stop(sprintf("Local matrix %d must be positive definite.", i),
           call. = FALSE)
    }
    value
  })
  local_log_determinant <- vapply(local, function(value) {
    as.numeric(determinant(value, logarithm = TRUE)$modulus)
  }, numeric(1))
  if (any(!is.finite(local_log_determinant))) {
    stop("Local matrices must have finite positive determinants.",
         call. = FALSE)
  }

  answer <- diag(1, n)
  if (n > 1L) {
    for (i in seq_len(n - 1L)) {
      for (j in seq.int(i + 1L, n)) {
        average <- (local[[i]] + local[[j]]) / 2
        chol_average <- tryCatch(chol(average), error = function(e) NULL)
        if (is.null(chol_average)) {
          stop("A local matrix average could not be factored.",
               call. = FALSE)
        }
        delta <- as.vector(coords[i, ] - coords[j, ])
        whitened <- forwardsolve(t(chol_average), delta)
        q <- sum(whitened * whitened)
        log_determinant_factor <- (local_log_determinant[i] +
          local_log_determinant[j]) / 4 -
          as.numeric(determinant(average, logarithm = TRUE)$modulus) / 2
        determinant_factor <- exp(log_determinant_factor)
        matern_value <- if (is.finite(q)) {
          .sp_matern_cor_unchecked(sqrt(q), 1, smoothness)
        } else {
          0
        }
        value <- determinant_factor * matern_value
        if (!is.finite(value) || value < 0 || value > 1 + 1e-12) {
          stop("Non-stationary correlation could not be evaluated stably.",
               call. = FALSE)
        }
        answer[i, j] <- answer[j, i] <- value
      }
    }
  }
  dimnames(answer) <- list(rownames(coords), rownames(coords))
  answer
}
