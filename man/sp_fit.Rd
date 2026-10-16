\name{sp_fit}
\alias{sp_fit}
\alias{coef.spmodels_fit}
\alias{fitted.spmodels_fit}
\alias{logLik.spmodels_fit}
\alias{print.spmodels_fit}
\alias{residuals.spmodels_fit}
\alias{summary.spmodels_fit}
\alias{print.summary.spmodels_fit}
\alias{vcov.spmodels_fit}
\title{Fit a Geostatistical Regression}
\usage{
sp_fit(
  formula,
  data,
  family = stats::gaussian(),
  coords = NULL,
  spatial = TRUE,
  time = NULL,
  random = NULL,
  intervention = NULL,
  weights = NULL,
  subset = NULL,
  na.action = stats::na.fail,
  offset = NULL,
  control = sp_control(),
  verbose = FALSE,
  random_slope = NULL
)
}
\arguments{
  \item{formula, data}{A model formula and either a data frame or projected
  point \code{sf} object.}
  \item{family}{A Gaussian identity, binomial logit, or Poisson log family.}
  \item{coords}{For a data frame, a one-sided formula selecting exactly two
  numeric Euclidean coordinate columns. Do not supply this for \code{sf} data.}
  \item{spatial}{Fit the spatial model? \code{FALSE} fits the independent
  Gaussian limiting model and does not require coordinates.}
  \item{random}{For any supported spatial family, an optional one-sided
  additive formula of independent latent Gaussian random-intercept grouping
  variables, such as \code{~ clinic + survey}.}
  \item{random_slope}{An optional grouped random-slope specification made by
  \code{sp_random_slope()}, or a list of additive specifications. Each
  covariance may be fixed or estimated through constrained Cholesky
  parameters.}
  \item{time}{Optional separable temporal specification created by
  \code{sp_temporal()}. Every supported spatial family may use a fixed or
  jointly estimated temporal scale.}
  \item{intervention}{Optional decay-adjusted intervention history created by
  \code{sp_intervention()} for a binomial-logit fit. Its impact parameters may
  be fixed or estimated by the labelled Laplace route.}
  \item{weights}{Integer binomial trial counts for a proportion response.
  Gaussian and Poisson weights are not implemented.}
  \item{subset}{Optional row subset following \code{model.frame} conventions.}
  \item{na.action}{Missing-value action; the default fails explicitly.}
  \item{offset}{Optional known linear-predictor offset.}
  \item{control}{An object returned by \code{sp_control()}.}
  \item{verbose}{Whether to trace optimisation.}
}
\value{An object of class \code{spmodels_fit}. Its standard methods include
\code{coef}, \code{vcov}, \code{fitted}, \code{residuals}, \code{logLik},
\code{summary}, \code{print}, and \code{predict}.}
\description{Fits Gaussian identity-link spatial regression by maximum
likelihood, or binomial-logit and Poisson-log latent Gaussian spatial
regression by a Laplace approximation. All models use a linear mean and
stationary isotropic Matern covariance with fixed smoothness, optionally
multiplied by exponential temporal correlation.}
\details{Data frames use two numeric Euclidean coordinate columns and require
explicit common units. Point \code{sf} data must have a known projected CRS;
longitude/latitude geometry is rejected rather than silently projected. New
\code{sf} prediction data must use the fitted CRS. Binomial responses may be
Bernoulli values, a two-column successes/failures matrix, or proportions with
integer trial weights. Poisson exposure should be supplied through a log
offset. Non-Gaussian fits support point and joint prediction; those predictive
calculations use a local Gaussian latent approximation at the fitted parameter
values. Latent Gaussian random intercepts add one independent effect per
observed level of each grouping variable. Prediction is conditional at known
levels by default; \code{unseen = "marginal"} in prediction methods integrates
a new zero-mean effect with its fitted variance. Positive-definite grouped
random-slope covariance matrices may additionally be supplied or estimated
with \code{sp_random_slope()}; correlated intercept-slope covariance is
carried through fitting, prediction, simulation, diagnostics and MCML.
Grouped effects are not
available for non-spatial fits. A
decay-adjusted intervention history may be supplied to a binomial-logit
spatial fit with fixed or constrained estimated impact parameters. A
\code{sp_temporal()} specification supplies one finite numeric time per row.
The separable correlation is the product of the Matern spatial correlation
and \code{exp(-abs(time_i-time_j)/temporal_scale)}. New prediction data must
contain the declared time variable. Temporal scale may be estimated only when
covariance estimation is enabled; Gaussian ML and non-Gaussian Laplace/MCML
fits carry it through prediction, simulation and validation.}
\examples{
d <- data.frame(y = c(1, 2, 1.5, 3), x = c(0, 1, 0, 1),
                east = c(0, 1, 0, 1), north = c(0, 0, 1, 1))
fit <- sp_fit(y ~ x, d, coords = ~ east + north,
              control = sp_control(spatial_variance = 1, scale = 1,
                                   nugget_variance = 0.2, estimate = FALSE,
                                   distance_units = "km"))
coef(fit)

survey <- transform(d, positive = c(0, 1, 0, 1), examined = 1)
binomial_fit <- sp_fit(
  cbind(positive, examined - positive) ~ x, survey,
  family = binomial(), coords = ~ east + north,
  control = sp_control(spatial_variance = 1, scale = 1,
                       nugget_variance = 0.2, estimate = FALSE,
                       distance_units = "km")
)
fitted(binomial_fit)
}
