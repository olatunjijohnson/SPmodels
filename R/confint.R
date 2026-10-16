#' Confidence intervals for fitted SPmodels parameters
#'
#' Construct local Wald intervals on an appropriate unconstrained scale:
#' identity for mean coefficients and random intercept-slope covariance
#' entries, logarithmic for positive covariance and decay parameters, and
#' logit for intervention impact.
#'
#' @param object A fitted `spmodels_fit` object.
#' @param parm Optional numeric or character selection from the qualified
#'   parameter names returned in the row names.
#' @param level Confidence probability strictly between zero and one.
#' @param ... Reserved for method compatibility.
#'
#' @return A two-column matrix of lower and upper confidence limits. Missing
#'   fitted standard errors remain missing limits.
#' @export
confint.spmodels_fit <- function(object, parm, level = 0.95, ...) {
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      !is.finite(level) || level <= 0 || level >= 1) {
    stop("`level` must be one number strictly between zero and one.",
         call. = FALSE)
  }
  table <- sp_model_table(object, include_smoothness = FALSE)
  table <- table[table$status == "estimated", , drop = FALSE]
  qualified <- paste(table$component, table$term, sep = "::")
  if (anyDuplicated(qualified)) {
    stop("Fitted parameter identifiers are not unique.", call. = FALSE)
  }
  row.names(table) <- qualified
  if (!missing(parm)) {
    if (is.character(parm)) {
      unknown <- setdiff(parm, qualified)
      if (length(unknown)) {
        stop(sprintf("Unknown parameter name(s): %s.",
                     paste(unknown, collapse = ", ")), call. = FALSE)
      }
      table <- table[parm, , drop = FALSE]
    } else if (is.numeric(parm) && all(is.finite(parm)) &&
               all(parm == as.integer(parm)) &&
               all(parm >= 1L & parm <= nrow(table))) {
      table <- table[as.integer(parm), , drop = FALSE]
    } else {
      stop("`parm` must contain valid parameter names or positions.",
           call. = FALSE)
    }
  }
  estimate <- table$estimate
  standard_error <- table$std_error
  z <- stats::qnorm((1 + level) / 2)
  lower <- upper <- rep(NA_real_, nrow(table))
  usable <- is.finite(estimate) & is.finite(standard_error) &
    standard_error >= 0
  slope_covariance <- table$component == "random_slope" &
    grepl("::cov(", table$term, fixed = TRUE)
  identity <- usable & (table$component == "mean" | slope_covariance)
  lower[identity] <- estimate[identity] - z * standard_error[identity]
  upper[identity] <- estimate[identity] + z * standard_error[identity]
  probability <- usable & table$component == "intervention" &
    table$term == "alpha" & estimate > 0 & estimate < 1
  if (any(probability)) {
    logit_se <- standard_error[probability] /
      (estimate[probability] * (1 - estimate[probability]))
    centre <- stats::qlogis(estimate[probability])
    lower[probability] <- stats::plogis(centre - z * logit_se)
    upper[probability] <- stats::plogis(centre + z * logit_se)
  }
  positive <- usable & !identity & !probability & estimate > 0
  if (any(positive)) {
    log_se <- standard_error[positive] / estimate[positive]
    centre <- log(estimate[positive])
    lower[positive] <- exp(centre - z * log_se)
    upper[positive] <- exp(centre + z * log_se)
  }
  probability_label <- format(100 * c((1 - level) / 2, (1 + level) / 2),
                              trim = TRUE, scientific = FALSE)
  answer <- cbind(lower, upper)
  colnames(answer) <- paste0(probability_label, "%")
  row.names(answer) <- row.names(table)
  attr(answer, "level") <- level
  attr(answer, "transform") <- ifelse(
    table$component == "mean" | slope_covariance, "identity",
    ifelse(table$component == "intervention" & table$term == "alpha",
           "logit", "log")
  )
  answer
}
