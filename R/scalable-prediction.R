#' Predict a large point set in bounded-memory chunks
#'
#' Apply the ordinary point-prediction method to consecutive row chunks and
#' combine the results in original order. This is an exact reorganization of
#' pointwise prediction, not a covariance approximation. It avoids holding the
#' complete training-by-prediction cross-covariance matrix in memory, but does
#' not return a joint covariance or joint draws across chunks.
#'
#' @param object A fitted `spmodels_fit` object.
#' @param newdata Non-empty prediction data containing every required
#'   predictor, coordinate, time and grouping column.
#' @param chunk_size Positive maximum number of prediction rows per chunk.
#' @param se.fit Return pointwise standard errors as well as fitted values?
#' @param target Predict the latent process or a future observation.
#' @param type Return response- or link-scale predictions.
#' @param unseen Policy for grouping levels absent during fitting.
#'
#' @return The same numeric vector or pointwise list schema as [predict()],
#'   with attributes `chunk_size` and `chunks` recording the computation.
#' @export
sp_predict_chunked <- function(object, newdata, chunk_size = 1000L,
                               se.fit = FALSE,
                               target = c("latent", "observation"),
                               type = c("response", "link"),
                               unseen = c("error", "marginal")) {
  if (!inherits(object, "spmodels_fit")) {
    stop("`object` must be a fitted `spmodels_fit` model.", call. = FALSE)
  }
  if (missing(newdata) || !is.data.frame(newdata) || !nrow(newdata)) {
    stop("`newdata` must be a non-empty data frame or `sf` object.",
         call. = FALSE)
  }
  if (!is.numeric(chunk_size) || length(chunk_size) != 1L ||
      is.na(chunk_size) || !is.finite(chunk_size) || chunk_size < 1 ||
      chunk_size != as.integer(chunk_size)) {
    stop("`chunk_size` must be one positive integer.", call. = FALSE)
  }
  if (!is.logical(se.fit) || length(se.fit) != 1L || is.na(se.fit)) {
    stop("`se.fit` must be TRUE or FALSE.", call. = FALSE)
  }
  target <- match.arg(target)
  type <- match.arg(type)
  unseen <- match.arg(unseen)
  chunk_size <- as.integer(chunk_size)
  chunk_id <- ceiling(seq_len(nrow(newdata)) / chunk_size)
  rows <- split(seq_len(nrow(newdata)), chunk_id)
  pieces <- lapply(rows, function(index) {
    stats::predict(object, newdata[index, , drop = FALSE], se.fit = se.fit,
                   target = target, type = type, unseen = unseen)
  })
  if (!se.fit) {
    answer <- unlist(pieces, use.names = FALSE)
    names(answer) <- row.names(newdata)
  } else {
    answer <- list(
      fit = unlist(lapply(pieces, `[[`, "fit"), use.names = FALSE),
      se.fit = unlist(lapply(pieces, `[[`, "se.fit"), use.names = FALSE),
      target = target
    )
    names(answer$fit) <- names(answer$se.fit) <- row.names(newdata)
    if (object$family$family != "gaussian") answer$type <- type
  }
  attr(answer, "chunk_size") <- chunk_size
  attr(answer, "chunks") <- length(rows)
  answer
}
