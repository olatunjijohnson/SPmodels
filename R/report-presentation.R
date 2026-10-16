.sp_format_number <- function(value, digits, p_value = FALSE) {
  out <- rep(NA_character_, length(value))
  finite <- is.finite(value)
  if (any(finite)) {
    out[finite] <- formatC(value[finite], format = "f", digits = digits)
    if (p_value) {
      threshold <- 10^(-digits)
      small <- finite & value < threshold
      out[small] <- paste0("<", formatC(threshold, format = "f",
                                          digits = digits))
    }
  }
  out
}

#' Format a stable SPmodels table for presentation
#'
#' @param object A data frame, fitted model, prediction-draw/target result, or
#'   validation/simulation result supported by the stable table functions.
#' @param digits Number of decimal places for ordinary numeric columns.
#' @param p_digits Number of decimal places and lower display threshold for
#'   `p_value` columns.
#' @param missing Character string used for missing values.
#' @param detail Validation-table detail passed to [sp_validation_table()].
#'
#' @return A character data frame of class `spmodels_formatted_table`. The
#'   unrounded source table is retained in the `raw_table` attribute.
#' @export
sp_format_table <- function(object, digits = 3L, p_digits = 3L,
                            missing = "\u2014",
                            detail = c("summary", "fold")) {
  validate_digits <- function(value, name) {
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        !is.finite(value) || value < 0 || value > 15 ||
        value != as.integer(value)) {
      stop(sprintf("`%s` must be one integer between 0 and 15.", name),
           call. = FALSE)
    }
    as.integer(value)
  }
  digits <- validate_digits(digits, "digits")
  p_digits <- validate_digits(p_digits, "p_digits")
  if (!is.character(missing) || length(missing) != 1L || is.na(missing)) {
    stop("`missing` must be one non-missing character string.", call. = FALSE)
  }
  detail <- match.arg(detail)
  raw <- if (is.data.frame(object)) {
    object
  } else if (inherits(object, "spmodels_fit")) {
    sp_model_table(object)
  } else if (inherits(object, c("spmodels_prediction_draws",
                                "spmodels_prediction_target"))) {
    sp_prediction_table(object)
  } else if (inherits(object, c("spmodels_cross_validation",
                                "spmodels_simulation_assessment",
                                "spmodels_simulation_recovery"))) {
    sp_validation_table(object, detail = detail)
  } else {
    stop("`object` is not supported by the SPmodels table interface.",
         call. = FALSE)
  }
  formatted <- as.data.frame(lapply(names(raw), function(name) {
    column <- raw[[name]]
    if (is.numeric(column)) {
      value <- .sp_format_number(column,
                                 if (name == "p_value") p_digits else digits,
                                 p_value = name == "p_value")
    } else {
      value <- as.character(column)
    }
    value[is.na(value)] <- missing
    value
  }), stringsAsFactors = FALSE, optional = TRUE)
  names(formatted) <- names(raw)
  row.names(formatted) <- row.names(raw)
  attr(formatted, "raw_table") <- raw
  attr(formatted, "digits") <- digits
  attr(formatted, "p_digits") <- p_digits
  attr(formatted, "missing") <- missing
  class(formatted) <- c("spmodels_formatted_table", "data.frame")
  formatted
}

#' @export
print.spmodels_formatted_table <- function(x, ..., row.names = FALSE) {
  print.data.frame(x, ..., row.names = row.names, right = TRUE)
  invisible(x)
}

#' Assemble an SPmodels analysis report
#'
#' @param object A fitted `spmodels_fit` object.
#' @param prediction Optional prediction-draw or prediction-target object.
#' @param validation Optional cross-validation, predictive-assessment or
#'   simulation-recovery object.
#' @param title Non-empty report title.
#' @param digits,p_digits Presentation precision passed to [sp_format_table()].
#' @param validation_detail Overall or fold-level validation table.
#'
#' @return An object of class `spmodels_report` containing unrounded stable
#'   tables and presentation-only formatted views.
#' @export
sp_report <- function(object, prediction = NULL, validation = NULL,
                      title = "SPmodels analysis report", digits = 3L,
                      p_digits = 3L,
                      validation_detail = c("summary", "fold")) {
  if (!inherits(object, "spmodels_fit")) {
    stop("`object` must be a fitted `spmodels_fit` model.", call. = FALSE)
  }
  if (!is.character(title) || length(title) != 1L || is.na(title) ||
      !nzchar(title)) {
    stop("`title` must be one non-empty character string.", call. = FALSE)
  }
  validation_detail <- match.arg(validation_detail)
  tables <- list(model = sp_model_table(object))
  if (!is.null(prediction)) tables$prediction <- sp_prediction_table(prediction)
  if (!is.null(validation)) {
    tables$validation <- sp_validation_table(validation,
                                             detail = validation_detail)
  }
  formatted <- lapply(tables, sp_format_table, digits = digits,
                      p_digits = p_digits)
  likelihood <- object$likelihood_method
  if (is.null(likelihood)) likelihood <- "Gaussian maximum likelihood"
  structure(list(
    title = title,
    metadata = data.frame(
      family = object$family$family,
      likelihood = likelihood,
      observations = object$nobs,
      spatial = isTRUE(object$spatial),
      distance_units = if (is.null(object$distance_units)) NA_character_ else
        object$distance_units,
      stringsAsFactors = FALSE
    ),
    tables = tables, formatted = formatted,
    digits = as.integer(digits), p_digits = as.integer(p_digits),
    call = match.call()
  ), class = "spmodels_report")
}

#' @export
print.spmodels_report <- function(x, ...) {
  cat(x$title, "\n", sep = "")
  cat(strrep("=", nchar(x$title)), "\n\n", sep = "")
  cat("Analysis\n")
  print(x$metadata, row.names = FALSE)
  labels <- c(model = "Model estimates", prediction = "Predictions",
              validation = "Validation")
  for (name in names(x$formatted)) {
    cat("\n", labels[[name]], "\n", sep = "")
    print(x$formatted[[name]])
  }
  invisible(x)
}
