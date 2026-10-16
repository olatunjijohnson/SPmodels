.sp_random_terms <- function(random) {
  if (is.null(random)) return(character())
  if (!inherits(random, "formula") || length(random) != 2L) {
    stop("`random` must be NULL or a one-sided additive formula such as `~ clinic + survey`.",
         call. = FALSE)
  }
  terms_object <- stats::terms(random)
  labels <- attr(terms_object, "term.labels")
  variables <- all.vars(random)
  if (!length(labels) || !identical(labels, variables)) {
    stop("`random` supports named grouping variables joined by `+`; interactions and transformations are not supported.",
         call. = FALSE)
  }
  if (anyDuplicated(variables)) {
    stop("Grouping variables in `random` must be unique.", call. = FALSE)
  }
  variables
}

.sp_random_training <- function(random, data) {
  variables <- .sp_random_terms(random)
  if (!length(variables)) return(NULL)
  missing <- setdiff(variables, names(data))
  if (length(missing)) {
    stop(sprintf("Grouping variables are missing from `data`: %s.",
                 paste(missing, collapse = ", ")), call. = FALSE)
  }
  groups <- lapply(variables, function(name) {
    value <- data[[name]]
    if (is.list(value) || is.matrix(value) || anyNA(value)) {
      stop(sprintf("Grouping variable `%s` must be an atomic vector without missing values.",
                   name), call. = FALSE)
    }
    factor(value, levels = unique(as.character(value)))
  })
  names(groups) <- variables
  levels <- lapply(groups, levels)
  z_parts <- lapply(seq_along(groups), function(j) {
    answer <- matrix(0, nrow(data), length(levels[[j]]),
                     dimnames = list(row.names(data),
                                     paste0(variables[j], "::", levels[[j]])))
    answer[cbind(seq_len(nrow(data)), as.integer(groups[[j]]))] <- 1
    answer
  })
  z <- do.call(cbind, z_parts)
  column_group <- rep(variables, lengths(levels))
  list(formula = random, variables = variables, levels = levels,
       groups = as.data.frame(groups, stringsAsFactors = FALSE),
       z = z, column_group = column_group)
}

.sp_random_new <- function(specification, data,
                           unseen = c("error", "marginal")) {
  if (is.null(specification)) return(NULL)
  unseen <- match.arg(unseen)
  missing <- setdiff(specification$variables, names(data))
  if (length(missing)) {
    stop(sprintf("Prediction data are missing grouping variables: %s.",
                 paste(missing, collapse = ", ")), call. = FALSE)
  }
  unseen_levels <- vector("list", length(specification$variables))
  names(unseen_levels) <- specification$variables
  z_parts <- lapply(seq_along(specification$variables), function(j) {
    name <- specification$variables[j]
    value <- as.character(data[[name]])
    if (anyNA(value)) {
      stop(sprintf("Grouping variable `%s` contains missing values.", name),
           call. = FALSE)
    }
    index <- match(value, specification$levels[[j]])
    if (anyNA(index) && unseen == "error") {
      unseen_values <- unique(value[is.na(index)])
      stop(sprintf("Grouping variable `%s` contains unseen levels: %s.",
                   name, paste(unseen_values, collapse = ", ")), call. = FALSE)
    }
    unseen_levels[[j]] <<- ifelse(is.na(index), value, NA_character_)
    answer <- matrix(0, nrow(data), length(specification$levels[[j]]),
                     dimnames = list(row.names(data),
                       paste0(name, "::", specification$levels[[j]])))
    known <- which(!is.na(index))
    answer[cbind(known, index[known])] <- 1
    answer
  })
  answer <- do.call(cbind, z_parts)
  attr(answer, "unseen_levels") <- unseen_levels
  answer
}

.sp_random_new_covariance <- function(z, specification, variance) {
  answer <- .sp_random_covariance(
    z, z, specification$column_group, variance
  )
  unseen_levels <- attr(z, "unseen_levels")
  if (is.null(unseen_levels)) return(answer)
  for (name in specification$variables) {
    value <- unseen_levels[[name]]
    is_same_unseen_level <- outer(value, value, function(left, right) {
      !is.na(left) & !is.na(right) & left == right
    })
    answer <- answer + unname(variance[[name]]) * is_same_unseen_level
  }
  answer
}

.sp_resolve_random_variance <- function(value, variables, estimate,
                                        default = 0.1) {
  if (!length(variables)) return(numeric())
  if (is.null(value)) {
    if (!estimate) {
      stop("Fixed grouped fitting requires `random_variance` in `sp_control()`.",
           call. = FALSE)
    }
    return(stats::setNames(rep(default, length(variables)), variables))
  }
  if (length(value) == 1L) value <- rep(value, length(variables))
  if (!is.null(names(value))) {
    if (anyDuplicated(names(value)) || !setequal(names(value), variables)) {
      stop("Named `random_variance` values must exactly match grouping variables.",
           call. = FALSE)
    }
    value <- value[variables]
  } else if (length(value) != length(variables)) {
    stop("`random_variance` must have length one or one value per grouping variable.",
         call. = FALSE)
  }
  stats::setNames(as.numeric(value), variables)
}

.sp_random_covariance <- function(z_left, z_right, column_group,
                                  variance) {
  if (is.null(z_left) || !length(variance)) {
    return(matrix(0, if (is.null(z_left)) 0L else nrow(z_left),
                  if (is.null(z_right)) 0L else nrow(z_right)))
  }
  weighted <- sweep(z_left, 2L, sqrt(variance[column_group]), "*")
  weighted_right <- sweep(z_right, 2L, sqrt(variance[column_group]), "*")
  tcrossprod(weighted, weighted_right)
}

#' Specify a grouped random slope
#'
#' Define a grouped Gaussian coefficient vector whose covariance is supplied
#' explicitly or used as the starting value for constrained estimation. Each
#' group receives an independent copy of that vector. With
#' `intercept = TRUE`, the vector contains a random intercept and a random
#' slope and `covariance` is therefore a 2 by 2 matrix. This explicit
#' covariance interface keeps assumptions visible in teaching; estimation is
#' opt-in rather than hidden behind a default.
#'
#' @param group Name of the grouping variable.
#' @param slope Name of the numeric random-slope variable.
#' @param covariance Positive-definite covariance matrix for the grouped
#'   coefficients.  Its order is intercept then slope when `intercept = TRUE`.
#' @param intercept Include a grouped random intercept as well as the slope?
#' @param label Optional unique label used in summaries.
#' @param estimate Estimate the covariance through a positive-definite
#'   Cholesky parameterization? The supplied matrix is then the starting value.
#'
#' @return An object for the `random_slope` argument of [sp_fit()].
#' @examples
#' fixed <- sp_random_slope(
#'   "clinic", "age_centered", matrix(c(0.3, 0.05, 0.05, 0.08), 2L)
#' )
#' estimated <- sp_random_slope(
#'   "clinic", "age_centered", matrix(c(0.3, 0.05, 0.05, 0.08), 2L),
#'   estimate = TRUE
#' )
#' @export
sp_random_slope <- function(group, slope, covariance, intercept = TRUE,
                            label = NULL, estimate = FALSE) {
  scalar_name <- function(value, argument) {
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
        !nzchar(value)) {
      stop(sprintf("`%s` must be one non-empty column name.", argument),
           call. = FALSE)
    }
    value
  }
  group <- scalar_name(group, "group")
  slope <- scalar_name(slope, "slope")
  if (identical(group, slope)) {
    stop("`group` and `slope` must name different columns.", call. = FALSE)
  }
  if (!is.logical(intercept) || length(intercept) != 1L || is.na(intercept)) {
    stop("`intercept` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(estimate) || length(estimate) != 1L || is.na(estimate)) {
    stop("`estimate` must be TRUE or FALSE.", call. = FALSE)
  }
  if (is.null(label)) label <- paste(group, slope, sep = ":")
  label <- scalar_name(label, "label")
  dimension <- if (intercept) 2L else 1L
  if (!is.matrix(covariance) || !is.numeric(covariance) ||
      !identical(dim(covariance), c(dimension, dimension)) ||
      anyNA(covariance) || any(!is.finite(covariance))) {
    stop(sprintf("`covariance` must be a finite numeric %d by %d matrix.",
                 dimension, dimension), call. = FALSE)
  }
  tolerance <- 100 * .Machine$double.eps * max(1, abs(covariance))
  if (max(abs(covariance - t(covariance))) > tolerance) {
    stop("`covariance` must be symmetric.", call. = FALSE)
  }
  covariance <- (covariance + t(covariance)) / 2
  if (is.null(tryCatch(chol(covariance), error = function(e) NULL))) {
    stop("`covariance` must be positive definite.", call. = FALSE)
  }
  coefficient_names <- if (intercept) c("(Intercept)", slope) else slope
  dimnames(covariance) <- list(coefficient_names, coefficient_names)
  structure(list(group = group, slope = slope, covariance = covariance,
                 intercept = intercept, label = label, estimate = estimate),
            class = "spmodels_random_slope")
}

#' @export
print.spmodels_random_slope <- function(x, ...) {
  cat("SPmodels grouped random slope",
      if (isTRUE(x$estimate)) "(estimated covariance)\n" else
        "(fixed covariance)\n")
  cat("Group:", x$group, " Slope:", x$slope,
      " Label:", x$label, "\n")
  cat(if (isTRUE(x$intercept))
    "Coefficient covariance (intercept, slope):\n" else
    "Slope variance:\n")
  print(x$covariance)
  invisible(x)
}

.sp_random_slope_list <- function(random_slope) {
  if (is.null(random_slope)) return(list())
  if (is.list(random_slope) && !length(random_slope)) return(list())
  if (inherits(random_slope, "spmodels_random_slope")) {
    random_slope <- list(random_slope)
  }
  if (!is.list(random_slope) || !length(random_slope) ||
      !all(vapply(random_slope, inherits, logical(1),
                  what = "spmodels_random_slope"))) {
    stop("`random_slope` must be created by `sp_random_slope()` or be a list of such objects.",
         call. = FALSE)
  }
  labels <- vapply(random_slope, `[[`, character(1), "label")
  if (anyDuplicated(labels)) {
    stop("Grouped random-slope labels must be unique.", call. = FALSE)
  }
  random_slope
}

.sp_random_slope_training <- function(random_slope, data) {
  specifications <- .sp_random_slope_list(random_slope)
  if (!length(specifications)) return(NULL)
  lapply(specifications, function(specification) {
    missing <- setdiff(c(specification$group, specification$slope), names(data))
    if (length(missing)) {
      stop(sprintf("Grouped random-slope columns are missing from `data`: %s.",
                   paste(missing, collapse = ", ")), call. = FALSE)
    }
    group <- data[[specification$group]]
    slope <- data[[specification$slope]]
    if (is.list(group) || is.matrix(group) || anyNA(group)) {
      stop(sprintf("Grouping variable `%s` must be an atomic vector without missing values.",
                   specification$group), call. = FALSE)
    }
    if (!is.numeric(slope) || is.matrix(slope) || anyNA(slope) ||
        any(!is.finite(slope))) {
      stop(sprintf("Random-slope variable `%s` must be finite numeric data without missing values.",
                   specification$slope), call. = FALSE)
    }
    group <- as.character(group)
    if (isTRUE(specification$estimate)) {
      group_count <- length(unique(group))
      minimum_groups <- if (specification$intercept) 3L else 2L
      if (group_count < minimum_groups) {
        stop(sprintf(
          "Estimated grouped random slope `%s` requires at least %d observed groups.",
          specification$label, minimum_groups), call. = FALSE)
      }
      if (specification$intercept) {
        variation <- tapply(slope, group, function(value) {
          length(unique(value)) >= 2L
        })
        if (any(!variation)) {
          stop(sprintf(
            "Estimated intercept-slope covariance `%s` requires at least two distinct slope values within every group.",
            specification$label), call. = FALSE)
        }
      } else {
        informative <- tapply(slope != 0, group, sum) >= 2L
        if (any(!informative)) {
          stop(sprintf(
            "Estimated slope-only covariance `%s` requires at least two non-zero slope observations within every group.",
            specification$label), call. = FALSE)
        }
      }
    }
    design <- if (specification$intercept) cbind(1, slope) else cbind(slope)
    colnames(design) <- colnames(specification$covariance)
    c(specification, list(group_values = group, levels = unique(group),
                          design = unname(design)))
  })
}

.sp_random_slope_identifiability <- function(specification,
                                              random_specification = NULL) {
  if (is.null(specification)) return(invisible(NULL))
  estimated <- specification[vapply(specification, function(component) {
    isTRUE(component$estimate)
  }, logical(1))]
  if (!length(estimated)) return(invisible(NULL))
  groups <- vapply(estimated, `[[`, character(1L), "group")
  if (anyDuplicated(groups)) {
    stop("Only one estimated random-slope covariance may use a given grouping variable.",
         call. = FALSE)
  }
  if (!is.null(random_specification) &&
      any(groups %in% random_specification$variables)) {
    overlap <- unique(groups[groups %in% random_specification$variables])
    stop(sprintf(
      "Estimated random-slope covariance cannot be combined with a separate random intercept for the same grouping variable: %s.",
      paste(overlap, collapse = ", ")), call. = FALSE)
  }
  invisible(NULL)
}

.sp_random_slope_new <- function(specification, data,
                                 unseen = c("error", "marginal")) {
  if (is.null(specification)) return(NULL)
  unseen <- match.arg(unseen)
  lapply(specification, function(reference) {
    missing <- setdiff(c(reference$group, reference$slope), names(data))
    if (length(missing)) {
      stop(sprintf("Prediction data are missing grouped random-slope columns: %s.",
                   paste(missing, collapse = ", ")), call. = FALSE)
    }
    group <- data[[reference$group]]
    slope <- data[[reference$slope]]
    if (is.list(group) || is.matrix(group) || anyNA(group)) {
      stop(sprintf("Grouping variable `%s` must be an atomic vector without missing values.",
                   reference$group), call. = FALSE)
    }
    if (!is.numeric(slope) || is.matrix(slope) || anyNA(slope) ||
        any(!is.finite(slope))) {
      stop(sprintf("Random-slope variable `%s` must be finite numeric data without missing values.",
                   reference$slope), call. = FALSE)
    }
    group <- as.character(group)
    new_levels <- setdiff(unique(group), reference$levels)
    if (length(new_levels) && unseen == "error") {
      stop(sprintf("Grouping variable `%s` contains unseen levels: %s.",
                   reference$group, paste(new_levels, collapse = ", ")),
           call. = FALSE)
    }
    design <- if (reference$intercept) cbind(1, slope) else cbind(slope)
    colnames(design) <- colnames(reference$covariance)
    reference$group_values <- group
    reference$design <- unname(design)
    reference
  })
}

.sp_random_slope_covariance <- function(left, right = left) {
  if (is.null(left) || is.null(right)) {
    n_left <- if (is.null(left)) 0L else length(left[[1L]]$group_values)
    n_right <- if (is.null(right)) 0L else length(right[[1L]]$group_values)
    return(matrix(0, n_left, n_right))
  }
  if (length(left) != length(right) ||
      !identical(vapply(left, `[[`, character(1), "label"),
                 vapply(right, `[[`, character(1), "label"))) {
    stop("Grouped random-slope specifications do not align.", call. = FALSE)
  }
  answer <- matrix(0, length(left[[1L]]$group_values),
                   length(right[[1L]]$group_values))
  for (j in seq_along(left)) {
    same_group <- outer(left[[j]]$group_values, right[[j]]$group_values,
                        `==`)
    coefficient_covariance <- left[[j]]$design %*%
      left[[j]]$covariance %*% t(right[[j]]$design)
    answer <- answer + same_group * coefficient_covariance
  }
  answer
}

.sp_random_slope_parameters <- function(specification) {
  if (is.null(specification)) return(numeric())
  values <- unlist(lapply(specification, function(component) {
    covariance <- component$covariance
    indices <- which(lower.tri(covariance, diag = TRUE), arr.ind = TRUE)
    terms <- ifelse(indices[, 1L] == indices[, 2L],
      paste0("var(", rownames(covariance)[indices[, 1L]], ")"),
      paste0("cov(", rownames(covariance)[indices[, 1L]], ",",
             colnames(covariance)[indices[, 2L]], ")"))
    stats::setNames(covariance[indices], paste(component$label, terms,
                                               sep = "::"))
  }), use.names = TRUE)
  stats::setNames(as.numeric(values), names(values))
}

.sp_random_slope_estimated_parameters <- function(specification) {
  if (is.null(specification)) return(numeric())
  values <- unlist(lapply(specification, function(component) {
    if (!isTRUE(component$estimate)) return(numeric())
    .sp_random_slope_parameters(list(component))
  }), use.names = TRUE)
  stats::setNames(as.numeric(values), names(values))
}

.sp_random_slope_set_parameters <- function(specification, values) {
  if (!length(values)) return(specification)
  if (is.null(names(values)) || anyDuplicated(names(values))) {
    stop("Grouped random-slope covariance updates must be uniquely named.",
         call. = FALSE)
  }
  available <- .sp_random_slope_parameters(specification)
  unknown <- setdiff(names(values), names(available))
  if (length(unknown)) {
    stop(sprintf("Unknown grouped random-slope covariance entries: %s.",
                 paste(unknown, collapse = ", ")), call. = FALSE)
  }
  for (j in seq_along(specification)) {
    component_names <- names(.sp_random_slope_parameters(
      list(specification[[j]])
    ))
    selected <- intersect(component_names, names(values))
    if (!length(selected)) next
    covariance <- specification[[j]]$covariance
    index <- which(lower.tri(covariance, diag = TRUE), arr.ind = TRUE)
    for (name in selected) {
      k <- match(name, component_names)
      row <- index[k, 1L]
      column <- index[k, 2L]
      covariance[row, column] <- covariance[column, row] <- values[[name]]
    }
    covariance <- (covariance + t(covariance)) / 2
    if (is.null(tryCatch(chol(covariance), error = function(e) NULL))) {
      stop(sprintf(
        "Grouped random-slope grid entries for `%s` do not form a positive-definite covariance matrix.",
        specification[[j]]$label), call. = FALSE)
    }
    specification[[j]]$covariance <- covariance
  }
  specification
}

.sp_random_slope_parameterization <- function(specification) {
  if (is.null(specification)) {
    return(list(values = numeric(), names = character(), components = list()))
  }
  estimated <- which(vapply(specification, function(component) {
    isTRUE(component$estimate)
  }, logical(1)))
  components <- lapply(estimated, function(j) {
    component <- specification[[j]]
    lower <- t(chol(component$covariance))
    index <- which(lower.tri(lower, diag = TRUE), arr.ind = TRUE)
    value <- vapply(seq_len(nrow(index)), function(k) {
      entry <- lower[index[k, 1L], index[k, 2L]]
      if (index[k, 1L] == index[k, 2L]) log(entry) else entry
    }, numeric(1L))
    coefficient <- rownames(component$covariance)
    term <- vapply(seq_len(nrow(index)), function(k) {
      if (index[k, 1L] == index[k, 2L]) {
        paste0("log_chol(", coefficient[index[k, 1L]], ")")
      } else {
        paste0("chol(", coefficient[index[k, 1L]], ",",
               coefficient[index[k, 2L]], ")")
      }
    }, character(1L))
    list(specification_index = j, index = index, value = value,
         names = paste("random_slope_chol", component$label, term, sep = "::"))
  })
  values <- unlist(lapply(components, `[[`, "value"), use.names = FALSE)
  names <- unlist(lapply(components, `[[`, "names"), use.names = FALSE)
  stats::setNames(values, names) -> values
  list(values = values, names = names, components = components)
}

.sp_random_slope_update <- function(specification, parameters,
                                    parameterization = NULL) {
  if (is.null(parameterization)) {
    parameterization <- .sp_random_slope_parameterization(specification)
  }
  if (!length(parameterization$values)) return(specification)
  if (length(parameters) != length(parameterization$values)) {
    stop("Grouped random-slope parameters do not align with their specification.",
         call. = FALSE)
  }
  cursor <- 0L
  for (details in parameterization$components) {
    count <- nrow(details$index)
    value <- parameters[cursor + seq_len(count)]
    cursor <- cursor + count
    dimension <- nrow(specification[[details$specification_index]]$covariance)
    lower <- matrix(0, dimension, dimension)
    for (k in seq_len(count)) {
      row <- details$index[k, 1L]
      column <- details$index[k, 2L]
      lower[row, column] <- if (row == column) exp(value[k]) else value[k]
    }
    covariance <- lower %*% t(lower)
    dimnames(covariance) <- dimnames(
      specification[[details$specification_index]]$covariance
    )
    specification[[details$specification_index]]$covariance <- covariance
  }
  specification
}

.sp_random_slope_natural_jacobian <- function(specification) {
  parameterization <- .sp_random_slope_parameterization(specification)
  parameter_count <- length(parameterization$values)
  natural <- .sp_random_slope_parameters(specification)
  jacobian <- matrix(0, length(natural), parameter_count,
                     dimnames = list(names(natural),
                                     names(parameterization$values)))
  if (!parameter_count) return(jacobian)
  cursor <- 0L
  for (details in parameterization$components) {
    component <- specification[[details$specification_index]]
    lower <- t(chol(component$covariance))
    natural_names <- names(.sp_random_slope_parameters(list(component)))
    if (nrow(lower) == 1L) {
      jacobian[natural_names, cursor + 1L] <- 2 * component$covariance[1L, 1L]
      cursor <- cursor + 1L
    } else {
      l11 <- lower[1L, 1L]
      l21 <- lower[2L, 1L]
      l22 <- lower[2L, 2L]
      columns <- cursor + seq_len(3L)
      jacobian[natural_names, columns] <- rbind(
        c(2 * l11^2, 0, 0),
        c(l11 * l21, l11, 0),
        c(0, 2 * l21, 2 * l22^2)
      )
      cursor <- cursor + 3L
    }
  }
  jacobian
}

.sp_random_slope_export <- function(specification) {
  lapply(specification, function(component) {
    structure(component[c("group", "slope", "covariance", "intercept",
                          "label", "estimate")],
              class = "spmodels_random_slope")
  })
}
