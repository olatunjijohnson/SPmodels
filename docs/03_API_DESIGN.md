# SPmodels proposed interface

## Design intent

Make a routine spatial fit feel like `glm`, with the compact workflow already used by SDALGCP2: a formula and data, optional controls, then `summary`, `predict` and `plot`. The installed SDALGCP2 0.1.1 help exposes `sdalgcp(formula, data, time, rasters, covariates, popden, control, verbose)`. Its public README demonstrates the same workflow. [SDALGCP2 README](https://github.com/olatunjijohnson/SDALGCP2)

All signatures and examples below are design proposals, not executable deliverables or settled API commitments. The provisional entry point is `sp_fit()` with `sp_control()` and an owned `spmodels_fit` class. These names remain open for review.

The package will support an expanding model collection for research and teaching. Keep common conventions for formulas, data, families, controls and result methods. Distinguish the response family from the spatial model class: `family = poisson()` alone cannot identify whether the model is a point-referenced Gaussian-process model, an areal model or a point-process model. Add explicit model selection or dedicated model-specific front ends when those models are designed, with their own required inputs and clear error messages.

## Main entry point

```r
sp_fit(
  formula, data, family = gaussian(),
  coords = NULL, spatial = TRUE, time = NULL,
  random = NULL, intervention = NULL,
  weights = NULL, subset = NULL, na.action = na.fail,
  offset = NULL, control = sp_control(), verbose = FALSE
)
```

`formula` describes the mean, including interactions, transformations and `offset(...)`. `data` accepts a point `sf` object, or a data frame with coordinates specified by a two-variable formula and a CRS supplied through controls. Point geometry supplies coordinates automatically. Polygon aggregation is a prediction target; polygon input must not silently be treated as point support.

Accept standard family objects, family functions and recognised character names. Start with Gaussian identity, binomial logit and Poisson log links. Inventory all reference custom inverse-link behaviours; map supported links to complete, validated family specifications or an advanced extension interface. Reject unsupported links explicitly. A custom inverse link alone does not supply all derivatives and distribution assumptions needed for fitting.

Use R's model-frame conventions for factors, contrasts, transformations, formula environments, subsetting and offsets. Preserve row identifiers when missing data are omitted or excluded. The proposed explicit `na.fail` default is a deliberate safety choice to revisit at kickoff, not a claim that it equals the usual glm default. [R glm documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/glm.html)

## Response and exposure conventions

| Response | Proposed input | Required interpretation |
| --- | --- | --- |
| Gaussian | `value ~ predictors` | Measurement scale and variance conventions recorded |
| Binomial | `cbind(successes, failures) ~ predictors` | Trials are the sum of the two columns |
| Binomial proportions | `proportion ~ predictors`, `weights = trials` | Equivalent likelihood to the two-column form when valid |
| Bernoulli | A zero or one response | One trial per observation |
| Poisson | `count ~ predictors + offset(log(exposure))` | Exposure enters once on the log scale |

Define valid weights for each family before implementation. Never treat generic weights, binomial denominators and Poisson exposure as interchangeable. Reject duplicate denominator specifications, negative counts, successes above trials, invalid exposure and inconsistent lengths. Decide and test zero-trial and zero-exposure behaviour explicitly.

## Example workflows

```r
# Point sf data with an explicit CRS
continuous <- sp_fit(value ~ altitude, data = samples,
                     family = gaussian())

prevalence <- sp_fit(
  cbind(positive, examined - positive) ~ rainfall + elevation,
  data = surveys, family = binomial()
)

counts <- sp_fit(
  cases ~ deprivation + offset(log(person_years)),
  data = surveillance, family = poisson()
)

summary(prevalence)
coef(prevalence)
predict(prevalence, newdata = grid_points,
        type = "response", se.fit = TRUE)
```

Keep a `spatial = FALSE` route for independent-model limiting checks. This does not authorise replacing the spatial fit with `glm` during ordinary use.

For Gaussian point data that cannot use a dense training covariance,
`sp_fit_neighbor()` is the explicit approximate front end. It uses the same
formula, coordinate and `sp_control()` conventions but requires a declared
neighbor budget and ordering. It owns a separate `spmodels_neighbor_fit` class
so an approximate likelihood cannot be mistaken for dense ML. The initial
route supports only Matérn-plus-nugget Gaussian models; unsupported grouped,
temporal and non-Gaussian structures fail explicitly.

## Spatial and advanced controls

`sp_control()` should expose named sections for covariance, optimisation, simulation, geometry and output retention. Record all resolved defaults in the fitted object. Retain controls for Matérn smoothness and scale, fixed or estimated variance components, initial parameters, convergence tolerances, latent draws, burn-in, thinning and random-number handling where reference parity requires them.

Require explicit distance units and a documented CRS policy. Do not guess a projected CRS or silently convert degrees to kilometres. Record transformations and the parameter conversion used to compare spatial scales.

Use `random = ~ cluster + survey` for additive independent intercepts and
`sp_random_slope()` for an explicit fixed positive-definite grouped slope or
correlated intercept-slope covariance. Both have tested known-level and
explicit marginal unseen-level policies. Do not imply arbitrary mixed-model
formula syntax or estimation of grouped slope covariance matrices.

Use `time` for the GP time index where supported. An `intervention` specification for DAST must separately carry survey time, MDA times, coverage/history, power and fixed or estimated impact parameters. The proposed design must allow survey time and the GP time index to differ and validate their separate roles against the mathematical specification.

Whether DAST uses this entry point or a dedicated `sp_fit_dast()` front end is an open design decision. Both must use the same fitting-object conventions.

## Fitted objects and methods

Use an owned class such as `spmodels_fit`, with `call`, `formula`, `terms`, family/link, model frame, contrasts, factor levels, row mapping, CRS/units, coefficient estimates, covariance parameters, uncertainty, log likelihood, convergence diagnostics, resolved controls and optional draws.

`coef()` should return fixed effects as a named vector, with a separate accessor for covariance and random-effect parameters. Document explicit mappings for existing consumer outputs. Supply `summary`, `predict`, `fitted`, `residuals`, `vcov`, `logLik`, `simulate`, `plot` and `update` only with documented statistical meanings. Do not inherit from `glm` solely to obtain its methods.

For `predict`, distinguish link-scale and response-scale means, latent process predictions, observation predictions, marginal and joint simulation, conditional and marginal random effects, standard errors and predictive standard deviations. State whether fitted parameter uncertainty is included. Nonlinear response means require integration, not merely applying the inverse link to a latent mean.

Plan a separate `predict_target()` for grid and polygon functionals, including weights, thresholds, quantiles and exceedance probabilities. Regional uncertainty requires joint spatial draws. Preserve the target definition and component inclusion flags.

## Teaching conventions

Give each model a short worked example that starts with the scientific question and data support, shows the formula and assumptions, fits the model, interprets estimates and uncertainty, checks assumptions and predicts. Provide small deterministic examples or bounded simulations for class use, with seeds and runtime guidance. Store resolved defaults so students can inspect them. Advanced controls should remain discoverable without making a first fit unnecessarily complicated. Link each documented method to its full MSc-accessible derivation, using the same parameter names and conventions in equations, help and summaries. Use appropriate lsda datasets for the worked examples and explain every transformation from the original data to the model inputs.

Use the same generic methods across model classes only where their statistical meanings are clear. Make distinctions such as latent versus observed outcomes, conditional versus marginal predictions, confidence versus credible intervals and point versus area support visible in help and summaries.

## Consumer integration

Use the owned SPmodels API and classes consistently. Update MBGapp and lsda calls to preserve the statistical meaning of inputs and outputs. If a temporary argument or result mapping is necessary, prefer locating it in consumer code. All computation must use the independent SPmodels implementation. Review any public integration helper against the distribution plan before including it in the submitted package.
