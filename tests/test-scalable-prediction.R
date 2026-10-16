if (file.exists(file.path("R", "scalable-prediction.R"))) {
  for (file in list.files("R", full.names = TRUE)) source(file)
} else {
  library(SPmodels)
}

close <- function(a, b, tolerance = 1e-11) {
  isTRUE(all.equal(a, b, tolerance = tolerance, check.attributes = FALSE))
}
expect_error <- function(expression, pattern) {
  error <- tryCatch({ force(expression); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}

dat <- data.frame(
  y = c(.8, 1.4, 1.7, 2.5, 2.2, 3.1, 3.4, 4),
  binary = c(0, 0, 0, 1, 0, 1, 1, 1),
  x = c(-1, -.5, 0, .5, 1, 1.5, 2, 2.5),
  clinic = rep(c("A", "B"), each = 4L),
  east = c(0, .5, 1, 1.5, 2, 2.5, 3, 3.5),
  north = c(0, .2, .8, 1, 1.7, 1.9, 2.7, 2.9)
)
control <- sp_control(
  spatial_variance = .7, scale = 1.3, nugget_variance = .15,
  estimate = FALSE, distance_units = "km"
)
random_slope <- sp_random_slope(
  "clinic", "x", matrix(c(.3, .05, .05, .12), 2L)
)
gaussian <- sp_fit(y ~ x, dat, coords = ~ east + north,
                   random_slope = random_slope, control = control)
new <- transform(dat[rep(seq_len(nrow(dat)), 3L), ],
                 east = east + rep(c(.05, .1, .15), each = nrow(dat)))
row.names(new) <- paste0("p", seq_len(nrow(new)))

ordinary <- predict(gaussian, new, se.fit = TRUE, target = "observation")
chunked <- sp_predict_chunked(
  gaussian, new, chunk_size = 5L, se.fit = TRUE,
  target = "observation"
)
stopifnot(close(unname(chunked$fit), ordinary$fit),
          close(unname(chunked$se.fit), ordinary$se.fit),
          identical(names(chunked$fit), row.names(new)),
          identical(attr(chunked, "chunks"), 5L),
          identical(attr(chunked, "chunk_size"), 5L))

ordinary_vector <- predict(gaussian, new)
chunked_vector <- sp_predict_chunked(gaussian, new, chunk_size = 7L)
stopifnot(close(unname(chunked_vector), ordinary_vector),
          identical(names(chunked_vector), row.names(new)))

binomial_fit <- sp_fit(binary ~ x, dat, family = binomial(),
  coords = ~ east + north, random_slope = random_slope, control = control)
ordinary_binary <- predict(binomial_fit, new, se.fit = TRUE,
                           type = "response")
chunked_binary <- sp_predict_chunked(
  binomial_fit, new, chunk_size = 3L, se.fit = TRUE, type = "response"
)
stopifnot(close(unname(chunked_binary$fit), ordinary_binary$fit),
          close(unname(chunked_binary$se.fit), ordinary_binary$se.fit),
          identical(chunked_binary$type, "response"))

# A shared unseen level may cross chunks: pointwise marginal results remain
# exact because cross-row covariance is not a returned target here.
unseen <- transform(new, clinic = "new-clinic")
ordinary_unseen <- predict(gaussian, unseen, se.fit = TRUE,
                           unseen = "marginal")
chunked_unseen <- sp_predict_chunked(
  gaussian, unseen, chunk_size = 4L, se.fit = TRUE, unseen = "marginal"
)
stopifnot(close(unname(chunked_unseen$fit), ordinary_unseen$fit),
          close(unname(chunked_unseen$se.fit), ordinary_unseen$se.fit))

expect_error(sp_predict_chunked(gaussian, new[0, ]), "non-empty")
expect_error(sp_predict_chunked(gaussian, new, chunk_size = 0), "positive")
expect_error(sp_predict_chunked(list(), new), "fitted")

cat("scalable point-prediction tests passed\n")
