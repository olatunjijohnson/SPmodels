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

## Installation

SPmodels is not yet on CRAN. To install the development version in desktop R,
run:

```r
install.packages("remotes")
remotes::install_github("olatunjijohnson/SPmodels")
library(SPmodels)
```

For a reproducible course, install the tagged course snapshot rather than
following the moving `main` branch:

```r
remotes::install_github(
  "olatunjijohnson/SPmodels",
  ref = "v0.1.0"
)
```

For browser-based R, install the prebuilt WebAssembly package after WebR has
been initialized. The published WebR repository currently tracks the latest
successful build; record the SPmodels commit and WebR version alongside course
materials. If your tutorial has a WebR instance named `webR`, run:

```js
await webR.installPackages("SPmodels", {
  repos: "https://olatunjijohnson.github.io/SPmodels/"
});
await webR.evalRVoid("library(SPmodels)");
```

See [WebR setup](docs/WEBR.md) for the complete example, package limitations,
and the live package repository. The WebR repository URL is for package
installation; the source code and user documentation are on GitHub.

Use `sp_predict_chunked(..., se.fit = TRUE)` for pointwise means and standard
errors over large grids. Reserve `sp_predict_draws()` for smaller sets where
joint spatial draws or draw-based target summaries are needed; its cost and
memory use rise quickly with grid size and number of draws.

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

From a shell in a local source checkout, you can also install the package with:

```sh
R CMD INSTALL .
```

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

GitHub Pages hosts the WebR package binaries. A separate searchable
teaching and documentation website remains future work.
