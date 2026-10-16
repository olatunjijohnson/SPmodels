# Binomial geostatistical teaching example: Loa loa prevalence

## Scientific question

How do elevation and vegetation relate to village-level *Loa loa* prevalence,
after allowing for residual spatial dependence?

This is a worked example of the implemented binomial-logit geostatistical
workflow. It illustrates model construction and diagnostics; it is not a
causal analysis or an operational disease-risk assessment.

## Data provenance and fixed input identity

- Package: `lsda` version 0.2.0.
- Object: `loaloa`.
- Package source inspected: 15 October 2026.
- Input-file SHA-256:
  `0acc9426c9842f93c264914a86e62641c160bdcdf134391cc57be7facb3c4811`.
- Support: 197 village survey locations in Cameroon and southern Nigeria.
- Response: `no_inf` positives among `no_exam` examined people.
- Predictors used here: elevation and mean 1999--2001 NDVI.
- Original coordinates: `longitude` and `latitude` in decimal degrees.
- Primary study citation recorded by lsda: Diggle et al. (2007), *Spatial
  modelling and prediction of Loa loa risk: decision making under
  uncertainty*.
- The dataset remains in lsda. It is neither copied into SPmodels nor made a
  runtime dependency.

The file hash makes the checked input unambiguous even if a later lsda release
changes the object.

## Prepare a projected analysis dataset

Longitude and latitude are angular coordinates, whereas the fitted covariance
uses Euclidean distances. The example therefore uses a Lambert azimuthal
equal-area projection centred on the study region and expressed in kilometres.
The projection is declared in full instead of relying on a machine-specific
default.

```r
library(SPmodels)
library(sf)

data("loaloa", package = "lsda")

loaloa$elevation_100m <- loaloa$elevation / 100
loaloa$mean_ndvi <- loaloa$mean9901

loaloa_sf <- st_as_sf(
  loaloa,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)

loaloa_laea <- st_transform(
  loaloa_sf,
  "+proj=laea +lat_0=7 +lon_0=12 +datum=WGS84 +units=km +no_defs"
)
```

Keeping `longitude` and `latitude` as ordinary columns preserves provenance;
the projected point geometry is what SPmodels uses for distances.

## Fit the prevalence model

For village \(i\), let \(Y_i\) be the number infected and \(m_i\) the number
examined. The fitted model is

\[
Y_i\mid S_i \sim \operatorname{Binomial}(m_i,p_i),\qquad
\operatorname{logit}(p_i)=
\beta_0+\beta_1\,\text{elevation}_{i}/100+
\beta_2\,\text{NDVI}_i+S_i.
\]

The latent field uses an exponential Matérn correlation because smoothness is
fixed at 0.5. The three supplied covariance values are starting values, not
fixed answers, because `estimate` retains its default `TRUE` value.

```r
fit <- sp_fit(
  cbind(no_inf, no_exam - no_inf) ~ elevation_100m + mean_ndvi,
  data = loaloa_laea,
  family = binomial(),
  control = sp_control(
    smoothness = 0.5,
    spatial_variance = 0.5,
    scale = 150,
    nugget_variance = 0.5,
    distance_units = "km",
    maxit = 120
  )
)

summary(fit)
confint(fit)
```

In the checked Linux run under R 4.6.1, optimization ended with convergence
code 0. The estimates were approximately -7.222 for the intercept, -0.127 per
100 metres of elevation and 11.557 per unit of mean NDVI. The estimated latent
spatial variance, Matérn scale and independent latent variance were 1.630,
79.479 km and 0.087. These numbers document reproducibility for the fixed input
above; they are not universal constants and do not substitute for scientific
interpretation or reference-environment validation.

## Diagnose and validate

Randomized quantile residuals account for the discreteness of binomial counts.
Their randomization is made reproducible without changing the caller's random
number stream.

```r
diagnostics <- sp_diagnostics(fit, seed = 20261015)
plot(diagnostics)

folds <- sp_spatial_folds(fit, folds = 3, seed = 20261015)
cv <- sp_cross_validate(fit, folds)
sp_validation_table(cv)
plot(cv)
```

The spatial folds are formed from unique locations, so observations at an
identical location cannot leak between training and test sets. For a stronger
separation check, use a scientifically chosen distance buffer:

```r
buffered <- sp_buffered_folds(
  fit,
  folds = 3,
  buffer = 25,
  seed = 20261015
)
buffered_cv <- sp_cross_validate(fit, buffered)
buffered_cv$training_size
```

The held-out log score integrates binomial probability over the approximate
Gaussian predictive link distribution. It is not a plug-in score evaluated
only at the mean prevalence.

## Predict prevalence jointly

This example uses observed locations to demonstrate the prediction object. A
substantive mapping analysis should construct and document an externally
validated prediction grid and covariate surfaces.

```r
draws <- sp_predict_draws(
  fit,
  newdata = loaloa_laea,
  nsim = 500,
  target = "latent",
  type = "response",
  parameter_uncertainty = TRUE,
  seed = 20261016
)

prevalence_summary <- sp_prediction_table(
  draws,
  probs = c(0.025, 0.5, 0.975),
  threshold = 0.20
)
head(prevalence_summary)
```

The draws retain spatial dependence across locations. With
`parameter_uncertainty = TRUE`, the local observed-Hessian approximation also
varies fitted mean and covariance parameters. This remains an approximation;
it is not a full Bayesian posterior.

## Scope and limitations

- the likelihood is a Laplace approximation unless `sp_mcml()` is run;
- Matérn smoothness is fixed rather than estimated;
- randomized residuals, spatial folds and predictive scores are model checks,
  not proof that the model is scientifically adequate;
- the projection is suitable for this regional teaching analysis, but every
  production study should justify its own distance representation;
- the checked run establishes internal execution and reproducibility, not
  numerical agreement with the separately managed validation environment;
- mapping requires reviewed prediction covariates and support, which are not
  supplied by this example.
