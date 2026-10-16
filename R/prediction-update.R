.sp_update_rows <- function(rows, data) {
  n <- nrow(data)
  ids <- row.names(data)
  if (is.null(rows)) return(seq_len(n))
  if (is.logical(rows)) {
    if (length(rows) != n || anyNA(rows)) {
      stop("Logical `rows` must contain one non-missing value per row.",
           call. = FALSE)
    }
    index <- which(rows)
  } else if (is.character(rows)) {
    if (anyNA(rows) || anyDuplicated(rows) || anyDuplicated(ids) ||
        any(!rows %in% ids)) {
      stop("Character `rows` must uniquely match prediction-data row names.",
           call. = FALSE)
    }
    index <- match(rows, ids)
  } else if (is.numeric(rows)) {
    if (anyNA(rows) || any(!is.finite(rows)) ||
        any(rows != as.integer(rows)) || any(rows < 1L | rows > n) ||
        anyDuplicated(rows)) {
      stop("Numeric `rows` must be unique valid row indices.", call. = FALSE)
    }
    index <- as.integer(rows)
  } else {
    stop("`rows` must be NULL, logical, numeric or character.", call. = FALSE)
  }
  if (!length(index)) stop("`rows` must select at least one row.", call. = FALSE)
  index
}

.sp_update_values <- function(value, index, ids, name, original) {
  if (!is.atomic(value) && !is.factor(value)) {
    stop(sprintf("Change `%s` must be an atomic vector or factor.", name),
         call. = FALSE)
  }
  if (length(value) == 1L) value <- rep(value, length(index))
  if (length(value) != length(index)) {
    stop(sprintf("Change `%s` must have length one or the number of selected rows.",
                 name), call. = FALSE)
  }
  if (!is.null(names(value))) {
    selected_ids <- ids[index]
    if (anyNA(names(value)) || anyDuplicated(names(value)) ||
        !setequal(names(value), selected_ids)) {
      stop(sprintf("Names for change `%s` must exactly match selected row names.",
                   name), call. = FALSE)
    }
    value <- value[match(selected_ids, names(value))]
  }
  if (is.factor(original)) {
    labels <- as.character(value)
    if (anyNA(labels) || any(!labels %in% levels(original))) {
      stop(sprintf("Change `%s` contains a new or missing factor level.", name),
           call. = FALSE)
    }
    value <- factor(labels, levels = levels(original), ordered = is.ordered(original))
  }
  value
}

#' Compare predictions under an updated predictor scenario
#'
#' Update selected predictor or in-formula offset variables in prediction data,
#' rebuild the fitted model's design and compare baseline with scenario
#' predictions. The fitted model, coordinates, geometry, row order and every
#' unmentioned data column are preserved.
#'
#' @param object A fitted `spmodels_fit` object.
#' @param newdata Prediction data accepted by [predict.spmodels_fit()].
#' @param changes A named list of replacement values for variables used on the
#'   right-hand side of the fitted formula, including variables inside an
#'   `offset()` term. Each value has length one or the number of selected rows.
#' @param rows Rows to update: `NULL` for every row, a logical vector, unique
#'   integer indices, or unique row names.
#' @param se.fit Include baseline and updated pointwise standard errors.
#' @param target Predict the latent process or a future observation.
#' @param type Return response- or link-scale predictions.
#'
#' @return An object of class `spmodels_prediction_update` containing updated
#'   data, baseline and updated design/offset values, row-level prediction
#'   differences and an explicit preservation audit.
#' @export
sp_prediction_update <- function(object, newdata, changes, rows = NULL,
                                 se.fit = FALSE,
                                 target = c("latent", "observation"),
                                 type = c("response", "link")) {
  if (!inherits(object, "spmodels_fit")) {
    stop("`object` must be a fitted `spmodels_fit` model.", call. = FALSE)
  }
  if (missing(newdata) || !is.data.frame(newdata) || !nrow(newdata)) {
    stop("`newdata` must be a non-empty data frame or `sf` object.",
         call. = FALSE)
  }
  if (isTRUE(object$external_offset)) {
    stop("Prediction updating requires offsets to be written in the formula.",
         call. = FALSE)
  }
  if (!is.list(changes) || !length(changes) || is.null(names(changes)) ||
      anyNA(names(changes)) || any(names(changes) == "") ||
      anyDuplicated(names(changes))) {
    stop("`changes` must be a non-empty uniquely named list.", call. = FALSE)
  }
  if (!is.logical(se.fit) || length(se.fit) != 1L || is.na(se.fit)) {
    stop("`se.fit` must be TRUE or FALSE.", call. = FALSE)
  }
  target <- match.arg(target)
  type <- match.arg(type)
  model_variables <- all.vars(stats::delete.response(object$terms))
  unknown <- setdiff(names(changes), model_variables)
  if (length(unknown)) {
    stop(sprintf("Changed variables are not in the fitted prediction design: %s.",
                 paste(unknown, collapse = ", ")), call. = FALSE)
  }
  coordinate_variables <- if (identical(object$coordinate_source, "columns"))
    all.vars(object$coords) else character()
  coordinate_changes <- intersect(names(changes), coordinate_variables)
  if (length(coordinate_changes)) {
    stop("Prediction scenarios cannot change coordinate variables.",
         call. = FALSE)
  }
  index <- .sp_update_rows(rows, newdata)
  ids <- row.names(newdata)
  baseline_design <- .sp_prediction_design(object, newdata)
  updated_data <- newdata
  for (name in names(changes)) {
    if (!name %in% names(updated_data)) {
      stop(sprintf("Prediction data do not contain changed variable `%s`.", name),
           call. = FALSE)
    }
    replacement <- .sp_update_values(changes[[name]], index, ids, name,
                                     updated_data[[name]])
    updated_data[[name]][index] <- replacement
  }
  updated_design <- .sp_prediction_design(object, updated_data)
  untouched <- setdiff(names(newdata), names(changes))
  columns_preserved <- identical(newdata[untouched], updated_data[untouched])
  coordinates_preserved <- identical(baseline_design$coordinates,
                                     updated_design$coordinates)
  geometry_preserved <- if (inherits(newdata, "sf")) {
    identical(sf::st_geometry(newdata), sf::st_geometry(updated_data))
  } else TRUE
  if (!columns_preserved || !coordinates_preserved || !geometry_preserved ||
      !identical(row.names(newdata), row.names(updated_data))) {
    stop("Internal preservation check failed while constructing the scenario.",
         call. = FALSE)
  }
  baseline <- stats::predict(object, newdata, se.fit = se.fit,
                             target = target, type = type)
  updated <- stats::predict(object, updated_data, se.fit = se.fit,
                            target = target, type = type)
  baseline_fit <- if (se.fit) baseline$fit else baseline
  updated_fit <- if (se.fit) updated$fit else updated
  comparison <- data.frame(
    row = ids, changed = seq_len(nrow(newdata)) %in% index,
    baseline = as.numeric(baseline_fit), updated = as.numeric(updated_fit),
    difference = as.numeric(updated_fit - baseline_fit),
    stringsAsFactors = FALSE
  )
  if (se.fit) {
    comparison$baseline_se <- as.numeric(baseline$se.fit)
    comparison$updated_se <- as.numeric(updated$se.fit)
  }
  structure(list(
    data = updated_data, comparison = comparison,
    baseline_design = baseline_design$x, updated_design = updated_design$x,
    baseline_offset = baseline_design$offset,
    updated_offset = updated_design$offset,
    changed_columns = names(changes), changed_rows = ids[index],
    preservation = c(row_order = TRUE, untouched_columns = columns_preserved,
                     coordinates = coordinates_preserved,
                     geometry = geometry_preserved,
                     fitted_parameters = TRUE),
    family = object$family$family, target = target, type = type,
    call = match.call()
  ), class = "spmodels_prediction_update")
}

#' @export
print.spmodels_prediction_update <- function(x, ...) {
  cat("SPmodels prediction scenario update\n")
  cat("Changed columns:", paste(x$changed_columns, collapse = ", "), "\n")
  cat("Changed rows:", length(x$changed_rows), "of", nrow(x$comparison), "\n")
  print(x$comparison, row.names = FALSE)
  invisible(x)
}
