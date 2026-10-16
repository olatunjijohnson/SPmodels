# Initial consumer audit

## Scope and evidence boundary

This is a read-only audit of the owner's downstream projects. It records the
scientific inputs and outputs that SPmodels must support; it is not an
authorization to modify those projects. No consumer source, Git state or
generated artefact was changed.

The source snapshots inspected on 14 October 2026 were:

| Consumer | Working source | Branch | Commit/version |
| --- | --- | --- | --- |
| lsda and course | `Documents/Manchester Job/Teaching/Courses/Longitudinal and Spatial Data Analysis/LSDA_2026` | `lda-slides-2026` | `e0541c11139244e190ecb2179efebb813f0ddc76`; lsda 0.2.0 |
| MBGapp | `Documents/Lancaster Job/projects/MBGApp/package/MBGapp` | `main` | `3243f40b1226e1b0d702614bc7ae82ccd3353210`; MBGapp 0.1.0 |
| SDALGCP2 | `Documents/Lancaster/Lancaster PhD Work/Lancaster PhD Work/RPackage/SDALGCP2` | `main` | `203bb67cc40767517d62f36b7152bf8f497f92f8`; SDALGCP2 0.1.1 |

MBGapp already contained modified and untracked files before inspection. Those
changes belong to the owner and were left untouched. Generated checks and
installed-library copies were excluded from the audit.

## MBGapp contract

### Fit inputs

The app constructs all three initial response families from a data frame and
two user-selected coordinate columns. Geographic input is transformed to one
projected UTM CRS chosen from the first input location; already projected input
is used directly. The user selects fixed Matérn smoothness, whether the nugget
is fixed at zero or estimated, optional covariates and nonlinear formula terms.

The required scientific schemas are:

| Outcome | Response | Exposure/trials | Mean model | SPmodels mapping |
| --- | --- | --- | --- | --- |
| Continuous | numeric response | none | intercept plus selected terms | Gaussian `sp_fit()` |
| Prevalence | infected count | examined count | intercept plus selected terms | `cbind(infected, examined - infected)` with binomial family |
| Prevalence teaching approximation | empirical logit | none | intercept plus selected terms | Gaussian `sp_fit()` after an explicit documented transformation |
| Count | non-negative integer count | positive exposure | intercept plus selected terms | Poisson `sp_fit()` with `offset(log(exposure))` in the formula |

The app currently accepts MCMC iteration, burn-in and thinning controls for its
stochastic route. SPmodels now has multiple-chain conditional latent sampling,
chain-specific acceptance, pooled ESS/MCSE and rank-normalized split R-hat.
Connecting these controls to an eventual app fitting adapter still requires an
authorized integration design; they must not be presented as parameter-chain
diagnostics when the sampler is conditional on supplied parameters.

### Downstream fields and results

The app directly prints `summary(fit)`. Its estimate table currently assumes
four implementation-specific summary components: regression coefficients,
spatial parameters, random-effect variance and measurement-error variance.
SPmodels must not reproduce those internal names. `sp_model_table()` now
provides its own stable schema with component, term, estimate, standard
error, statistic, p-value, estimation status and distance units.

Prediction uses an uploaded or regular point grid, optionally accompanied by a
same-row predictor data frame. It requires joint draws on the link scale and
then applies identity, inverse-logit or exponential transformations. The
resulting application table is two coordinate columns followed by one column
per predictive draw. Maps display:

- draw mean;
- draw standard deviation;
- the fraction of draws above a threshold; and
- a selected draw quantile.

SPmodels already supplies joint response or link draws with
`sp_predict_draws()`, stable target summaries with `sp_prediction_table()`, and
weighted grouped summaries with `sp_predict_target()`. A future MBGapp adapter
can construct its raster input as `cbind(coordinates, t(draws$draws))` without
reading fit internals. The adapter must also carry the prediction CRS and
geographic-display CRS explicitly; SPmodels does not infer them from app state.

### Smallest migration tests

1. Continuous: fit a projected six-point fixture, predict on a four-point grid,
   and verify table schema, coordinate order, mean, SD, quantile and exceedance
   values against direct draw calculations.
2. Prevalence: fit successes/failures, verify probability draws lie in
   `[0,1]`, and verify the same grid-summary calculations.
3. Count: fit an in-formula log exposure offset, verify non-negative intensity
   draws and preservation of the prediction-grid exposure.
4. App smoke test: for each family, reach fit summary, table, prediction and
   raster construction without implementation-specific object access.

## lsda and course contract

The package requires R 4.1 or later. The audited week 8 lesson currently uses
projected point `sf` data and teaches a formula/data workflow, ordinary
`summary()`, grid prediction, predictive mean/SD and threshold-exceedance maps.
The lesson sequence is Gaussian fitting first, followed by prediction,
binomial fitting and a further continuous worked example.

Relevant data schemas are:

| Object | Support and fields | Use |
| --- | --- | --- |
| `galicia` | 132 points; `loglead`, longitude, latitude and `pm10` | Gaussian spatial regression |
| `galicia_grid` | 1,508 projected grid coordinates | prediction locations |
| `galicia_pred` | 1,508 aligned rows with `pm10`, longitude and latitude | prediction covariate |
| `loaloa` | 197 villages; infected, examined, coordinates, elevation and vegetation summaries | binomial spatial regression |

The week 8 Gaussian exercises compare no covariate/no nugget, no
covariate/estimated nugget and covariate/no nugget models with smoothness fixed
at 1.5. Prediction requires the grid and predictor rows to agree and asks
students to distinguish predictive mean from predictive uncertainty. The
binomial lesson uses infected and examined counts and stresses seed sensitivity
in stochastic estimation. The exceedance exercise calculates the proportion of
predictive samples above a scientifically chosen threshold.

Week 7 additionally requires an empirical residual variogram, pair summaries,
a no-spatial-correlation permutation envelope and a global permutation
p-value. SPmodels now supplies all four numerical components through
`sp_distance_summary()`, `sp_variogram()` and `sp_variogram_envelope()`; an
authorized course migration and reviewed lesson wording remain outstanding.

The smallest course migration test is a non-interactive script that loads lsda
data explicitly, transforms geographic points to the declared projected CRS,
fits the three Gaussian models, predicts on the aligned grid, computes the four
draw summaries, and fits one binomial model. Expected teaching values must be
regenerated and reviewed rather than copied from another implementation.
Student installation, PDF rendering and runtime on the supported classroom
platforms remain release gates.

## Interface-style evidence

SDALGCP2 confirms the owner's preferred style: one compact formula/data/control
entry point, ordinary `print()`, `summary()`, `predict()` and `confint()`
methods, and spatial prediction results that can be joined directly to `sf`
geometry. SPmodels already follows the formula and ordinary-method parts and
now provides `confint.spmodels_fit`. It does not return an `sf` prediction
object directly; instead, `sp_target_geometry()` performs an explicit
identifier-safe join to point or polygon geometry. Broader polygon/time model
classes remain future requirements. These are design requirements, not
permission to copy implementation details.

## SPmodels-owned contract added from this audit

The consumer contract is now executable in `tests/test-reporting.R`:

- `sp_model_table()` exposes estimates without internal-field coupling;
- `sp_prediction_table()` exposes point or grouped draw summaries;
- `sp_validation_table()` exposes overall/fold cross-validation, predictive
  assessment and simulation-recovery metrics in one long schema; and
- the tests recompute reported values directly from owned SPmodels objects.

Remaining consumer work is deliberately separated from package correctness.
Integration edits require owner authorization; reference comparisons run in a
separate pinned environment; and consumer projects must remain independently
installable.
