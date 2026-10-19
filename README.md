# SPmodels

SPmodels is an independent R package for spatial modelling in research and
teaching. It is designed around a familiar formula interface, explicit model
assumptions, reproducible examples, and full mathematical derivations for
implemented methods.

> **Development status:** active development. The package has passed
> its current Linux `R CMD check`, but it is not yet a CRAN release. Some
> capabilities are implemented while reference comparisons, cross-platform
> checks, consumer migrations, and release work remain open.

SPmodels is released under the [MIT License](LICENSE-MIT).

## What it can do

- Fit Gaussian, binomial-logit, and Poisson-log spatial regression models.
- Simulate fitted responses; make point and joint predictions; build spatial
  grids and weighted polygon targets.
- Model grouped random effects, intervention histories, and separable
  space–time dependence within documented limits.
- Run diagnostics, cross-validation, simulation recovery, and reporting
  workflows.
- Build a valid non-stationary Matérn correlation matrix from supplied local
  range/anisotropy matrices. This new covariance primitive is implemented, but
  local covariance estimation and a fitted non-stationary regression workflow
  are not yet available.

All model claims are bounded by the tested configurations and documented
limitations. In particular, no universal nearest-neighbour budget or MCML
proposal setting is recommended. Study results and limitations are linked in
the guide and validation reports below.

## Quick start

The following is the basic stationary Gaussian workflow. Coordinates must be
projected Euclidean coordinates in common units.

```r
fit <- sp_fit(
  outcome ~ rainfall + elevation,
  data = survey_points,
  coords = ~ east_km + north_km,
  control = sp_control(distance_units = "km")
)

summary(fit)
predict(fit, newdata = prediction_points, se.fit = TRUE)
```

From a shell in a local source checkout, install the package with:

```sh
R CMD INSTALL .
```

For browser-based use, SPmodels is built as a WebAssembly package repository.
Once the first deployment is complete, the WebR install URL and JavaScript
example are documented in [WebR setup](docs/WEBR.md).

## Non-stationary covariance

The first non-stationary capability is `sp_nonstationary_cor()`. It accepts
one symmetric positive-definite 2 × 2 local scale matrix per projected point
and returns a valid Matérn correlation matrix. The local matrices control
range and directional anisotropy. See the [full derivation and current
boundary](docs/derivations/21_NONSTATIONARY_MATERN_COVARIANCE.md) before using
it. This function does not fit the local matrices or a regression model.

## Documentation

- [Implemented guide (PDF)](docs/SPmodels_Implemented_Guide.pdf) — how the
  current package works, its assumptions, and what remains.
- [Feature inventory](docs/02_FEATURE_INVENTORY.md) — the agreed functional
  scope.
- [Implementation register](docs/11_IMPLEMENTATION_REGISTER.md) — methods,
  evidence, and limitations.
- [Completion checklist](docs/13_PROJECT_COMPLETION_TASKS.md) — open research,
  validation, and release work.
- [All validation reports](docs/validation/README.md) — retained study designs
  and outcomes.
- [WebR / WebAssembly setup](docs/WEBR.md) — build and smoke-check status,
  hosting boundary, installation instructions, and browser limitations.
- [Derivations](docs/derivations/) — equation-level chapters, including the
  new [non-stationary Matérn derivation](docs/derivations/21_NONSTATIONARY_MATERN_COVARIANCE.md).
- [Technical derivations (PDF)](docs/SPmodels_Technical_Derivations.pdf) and
  [derivation collection (PDF)](docs/SPmodels_Derivations.pdf).

GitHub Pages will host the WebR package binaries after its first successful
deployment. A separate searchable
teaching and documentation website remains future work.
