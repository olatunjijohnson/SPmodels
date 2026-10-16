# Gaussian geostatistical teaching example: Galicia moss lead

## Scientific question

How is log-transformed lead concentration in moss associated with co-located
PM10, and does residual spatial dependence remain after including PM10?

This example demonstrates the currently implemented Gaussian workflow. It is
not presented as a completed substantive environmental analysis.

## Data provenance and version

- Package: `lsda` version 0.2.0.
- Object: `galicia`.
- Package source inspected: 13 October 2026.
- Support: 132 point sampling stations across Galicia, north-west Spain.
- Response: `loglead`, documented by inspection as the base-10 logarithm of
  `lead` (for example, `log10(1.9)` equals the first recorded value to the
  displayed precision).
- Predictor: co-located `pm10`, particulate matter below 10 micrometres.
- Original coordinates: `long` and `lat` in decimal degrees.
- The dataset remains in lsda. It is not copied into SPmodels and is not a
  runtime dependency of SPmodels.

The lsda documentation does not presently give a bibliographic source beyond
the study description. That missing upstream citation must be resolved before
the example is used in a release claiming complete provenance.

## Prepare projected point data

Spatial distances must not be computed directly from decimal degrees. Galicia
falls within UTM zone 29N for these locations, so the example projects WGS84
longitude/latitude to ETRS89 / UTM zone 29N (EPSG:25829), whose coordinate unit
is metres.

```r
library(SPmodels)
library(sf)

data("galicia", package = "lsda")

galicia_sf <- st_as_sf(
  galicia,
  coords = c("long", "lat"),
  crs = 4326,
  remove = FALSE
)

galicia_utm <- st_transform(galicia_sf, 25829)
```

The transformation is explicit and inspectable. SPmodels will reject the
unprojected longitude/latitude object.

## Compare independent and spatial models

```r
independent_fit <- sp_fit(
  loglead ~ pm10,
  data = galicia,
  spatial = FALSE
)

spatial_fit <- sp_fit(
  loglead ~ pm10,
  data = galicia_utm,
  control = sp_control(smoothness = 0.5)
)

summary(independent_fit)
summary(spatial_fit)
```

The coefficient describes the expected change in base-10 log lead for a
one-unit increase in PM10. It must not be interpreted causally without a study
design and confounding argument.

The spatial fit decomposes residual variation into a Matérn spatial component
and an independent nugget. Its scale is reported in metres because that is the
CRS unit. Smoothness is fixed at 0.5 in the current implementation.

## Check residual dependence

```r
vg <- sp_variogram(spatial_fit, bins = 10)
vg$values
plot(vg)

diagnostics <- sp_diagnostics(spatial_fit)
plot(diagnostics)
```

The empirical variogram uses residuals after the fitted regression mean by
default. Each bin's pair count must be considered: a value based on few pairs
is unstable. The fitted curve is a diagnostic overlay, not a second fit.

The decorrelated-residual quantile plot checks whether the fitted Gaussian mean
and covariance have removed major distributional structure. It is not proof of
model adequacy.

## Reproducible model simulation

```r
simulated <- simulate(spatial_fit, nsim = 20, seed = 20261013)
dim(simulated)
```

These are unconditional fitted-model response replicates at the 132 observed
sites. They condition on the estimated parameters and therefore do not include
parameter uncertainty.

## Current limitations and available extensions

- covariance smoothness is fixed;
- covariance standard errors are local Hessian/delta-method approximations;
- ordinary point standard errors condition on fitted covariance parameters;
  `sp_predict_draws(..., parameter_uncertainty = TRUE)` supplies an optional
  local observed-Hessian approximation when joint draws are needed;
- formal variogram envelopes, spatial or buffered cross-validation, integrated
  predictive scoring, calibration and recovery tools are available, but
  should be added only with analysis-specific choices such as the buffer and
  simulation count;
- `galicia_pred` and `galicia_grid` do not expose a documented exact coordinate
  alignment in the inspected lsda metadata. They must not be joined by row
  position or column-name guesswork; grid prediction remains blocked pending
  a verified alignment record;
- reference-environment comparison has not yet been run.
