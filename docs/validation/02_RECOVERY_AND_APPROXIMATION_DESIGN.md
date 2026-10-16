# Predeclared recovery and approximation study design

## Status and purpose

**Status: NNL statistical accuracy/stability, per-fit timing and whole-run peak
RSS recorded; GRS results recorded with mixed/inconclusive precision outcomes;
TSC recovery and MCP proposal-policy results are recorded, with no tested
MCML policy operational; the NNP accuracy study is recorded;
per-configuration external RSS and cross-platform NNP measurements remain
outstanding.**

This design fixes the local simulation and approximation studies for the
newly estimated grouped covariance, temporal scale, adaptive MCML and
nearest-predecessor prediction and Gaussian training-likelihood routes before
their final results are examined.
It complements the separately managed reference-comparison matrix; it neither
freezes that external environment nor substitutes internal recovery for
reference validation.

The NNL accuracy/budget-stability, GRS grouped-slope, TSC temporal-scale,
synthetic MCP proposal-policy and NNP prediction studies in this design were
executed on 16 October 2026. NNL results are
reported in `docs/benchmarks/04_NEIGHBOR_LIKELIHOOD.md`; GRS results and the
case-specific interpretation are in `03_GROUPED_SLOPE_RECOVERY_RESULTS.md`;
TSC results are in `04_TEMPORAL_SCALE_RECOVERY_RESULTS.md`; MCP results are in
`05_MCML_PROPOSAL_RESULTS.md`; NNP results are in
`06_NEIGHBOR_PREDICTION_RESULTS.md`.

Every run must record the SPmodels archive SHA-256, R/platform/BLAS details,
full configuration, seeds, warnings, convergence codes, raw replicate values,
elapsed time and peak memory. Failed fits remain in the denominator and record.
The first NNL execution retained its accuracy metrics but did not capture
performance metadata. A second complete run records elapsed time for each fit
and externally measured peak RSS for the study process as a whole. It does not
provide per-configuration externally measured RSS; see the benchmark report.

## Common rules

- Synthetic response seeds are `410001 + scenario_index * 10000 + replicate`.
- Any MCML refit seed is the response seed plus 5,000.
- Scenario order is the table order below; no scenario is dropped after seeing
  results.
- Recovery begins with 100 successful refits, then uses the already implemented
  geometric schedule 200 and 400. It stops only when every declared bias-MCSE
  target passes or 400 attempted replicates have been evaluated.
- A convergence warning is not silently removed. Summaries report attempted,
  successful and failed counts separately.
- Wald coverage is descriptive near variance boundaries. It is not used to
  validate boundary-corrected inference, which is not implemented.
- Exact analytic/small-matrix checks remain primary when available.

For parameter \(j\), mean bias and its simulation standard error are

\[
 \widehat b_j=R_s^{-1}\sum_{r\in\mathcal S}
 (\widehat\theta_{rj}-\theta_j),\qquad
 \operatorname{MCSE}(\widehat b_j)=
 \operatorname{sd}(\widehat\theta_{rj}-\theta_j)/\sqrt{R_s},
\]

where \(\mathcal S\) contains the finite successful fits and \(R_s\) is its
size. A study is inconclusive if fewer than 90% of attempted fits are usable,
even when the successful subset meets its bias-MCSE target.

## GRS: estimated grouped random-slope covariance

Each scenario has groups with centered slope values \(-1,0,1\), an intercept
and slope fixed mean, exponential spatial correlation, and independent
Gaussian grouped coefficient vectors. Coordinates and latent/response draws
are generated once per declared response seed.

| ID | Groups | True grouped covariance \(G\) | Spatial variance / scale / nugget | Purpose |
|---|---:|---|---|---|
| GRS-01 | 8 | `diag(0.50, 0.20)` | 0.30 / 3 / 0.15 | zero covariance interior check |
| GRS-02 | 8 | `matrix(c(.70,.18,.18,.25),2)` | 0.30 / 3 / 0.15 | moderate positive correlation |
| GRS-03 | 20 | same as GRS-02 | 0.30 / 3 / 0.15 | increased-group information |
| GRS-04 | 20 | `matrix(c(.30,-.12,-.12,.25),2)` | 0.30 / 3 / 0.15 | negative covariance |

All spatial and grouped covariance parameters are estimated from the same
declared starts: `(0.4, 2.5, 0.2)` and `diag(0.4, 0.2)`. Targets are the three
natural entries of \(G\), fixed effects, spatial components, positive-
definiteness, convergence and 95% local-Wald coverage.

To make the simulation reproducible, each group contributes one observation
at each centered slope value `-1, 0, 1`. Group centres occupy the first `G`
cells of a regular grid with 2-unit spacing; each observation receives
independent coordinate jitter uniform on `[-0.1, 0.1]` in each axis. The fixed
mean is `0.3 + 0.4 * z`. For each response seed, draw one independent Gaussian
intercept/slope vector per group with covariance \(G\), then draw the spatial
field and nugget jointly as a Gaussian vector with covariance
`0.30 * exp(-distance / 3) + 0.15 * I`; add these components to the fixed
mean. This supplies within-group slope replication while retaining the
spatial component.

For this GRS table, `scenario_index` is the one-based row number (1--4) and
replicate indices start at 1. Thus response seeds follow the common formula
`410001 + 10000 * scenario_index + replicate`. At each stage, attempt seeds
until the target of 100, 200 or 400 successful fits is reached, or 400 total
attempts have been used. A fit is usable for point recovery only if it returns
finite target estimates, positive-definite \(G\), zero optimizer convergence
code and no warning. Record every warning and failure. Hessian-based
uncertainty availability and local-Wald coverage are reported separately;
missing Hessian uncertainty does not remove finite converged estimates from
the bias-MCSE calculation. Evaluate bias-MCSE targets after each
successful-fit target; stop when all pass. Fewer than 90% usable attempts is
inconclusive.

Bias-MCSE targets are 0.03 for grouped variances, 0.02 for grouped covariance,
0.03 for fixed effects and 0.05 times the true value (with floor 0.02) for
positive spatial parameters. The primary mathematical checks reconstruct
\(LL^{\mathsf T}\) and the analytic natural-scale Jacobian for every finite fit.
The execution results are in `03_GROUPED_SLOPE_RECOVERY_RESULTS.md`: GRS-01/02
did not meet the scale bias-MCSE target at 400 fits; GRS-03/04 met the
bias-MCSE precision thresholds, but their mean spatial-scale estimates were
still far from the generating value. No scenario is described as full model
recovery, and local-Wald coverage was poor where intervals were available.

## TSC: estimated non-Gaussian temporal scale

Spatial locations form a 4 by 4 grid observed at times 0, 1, 2 and 4. Temporal
correlation is exponential and separable from an exponential Matérn spatial
correlation. The same covariate pattern is reused across response families.

| ID | Family | True temporal scale | Spatial variance / scale / nugget | Mean specification |
|---|---|---:|---|---|
| TSC-01 | binomial, 10 trials | 0.75 | 0.45 / 2 / 0.20 | intercept -0.5, slope 0.35 |
| TSC-02 | binomial, 10 trials | 2.00 | 0.45 / 2 / 0.20 | same |
| TSC-03 | Poisson log | 0.75 | 0.30 / 2 / 0.15 | intercept 0.2, slope 0.20 |
| TSC-04 | Poisson log | 2.00 | 0.30 / 2 / 0.15 | same |

Laplace fitting estimates spatial components and temporal scale. Starts are
1.25 times the true temporal scale and `(0.4, 1.5, 0.2)` for spatial values.
Targets are log- and natural-scale temporal estimates, fixed effects,
convergence, covariance positive definiteness and local-Wald coverage.
Temporal bias MCSE must be at most `max(0.05, 0.05 * true_scale)`.

For the first 25 successful replicates of each scenario, normalized MCML is
also run at 1,000 draws with three independent seeds per response. Laplace--
MCML movement, between-MCML covariance, ESS fraction and relative likelihood
MCSE are retained; no MCML estimate is called stable unless every repeat has
ESS fraction at least 0.10 and relative MCSE at most 0.05.

The completed TSC results are in `04_TEMPORAL_SCALE_RECOVERY_RESULTS.md`.
Every case reached 400 attempts but had fewer than 90% usable fits, temporal-
scale bias-MCSE far above target, frequent optimizer convergence warnings and
poor MCML repeat-gate rates. No case is considered recovered or validated.

The exact shared design uses the Cartesian product of east/north coordinates
`0:3` and times `0, 1, 2, 4`, in that order. The deterministic covariate is
the standardized value of `east + 0.5 * north + 0.25 * time`, reused unchanged
in all four cases. For each response seed, draw a zero-mean latent Gaussian
vector with covariance
`spatial_variance * exp(-distance / spatial_scale) * exp(-time_lag / temporal_scale) + nugget_variance * I`.
Add it to the declared fixed linear predictor. Draw binomial successes from
`Binomial(10, plogis(eta))` or Poisson counts from `Poisson(exp(eta))` as
appropriate. The spatial smoothness is 0.5. Scenario indices are one-based
table row numbers and replicates start at 1, using the common response-seed
formula.

Use the registered spatial starts `(0.4, 1.5, 0.2)` and an initial temporal
scale 1.25 times its truth. Laplace refits use the ordinary covariance
optimizer defaults. For MCML, call the independent-repeat route with 1,000
draws, default optimizer settings, three consecutive seeds beginning at
`response_seed + 5,000`, and no proposal adaptation. Capture each repeat's
optimizer status and Monte Carlo diagnostics as well as the between-repeat
natural-parameter covariance.

## MCP: MCML proposal policies

The binomial fixture has 24 sites on a 6-by-4 integer grid; the Poisson fixture
has 48 sites on an 8-by-6 integer grid. In each, `z` is the standardized value
of `east + 0.5 * north`. Generate a Gaussian spatial field and independent
nugget jointly with covariance `0.40 * exp(-distance / 2) + 0.20 * I`, using
smoothness 0.5. The binomial mean is `-0.4 + 0.3*z` with 8 trials; the Poisson
mean is `0.2 + 0.2*z`. Response seeds are 510001 and 510002 in case order.
Fit one common starting model per fixture using spatial starts `(0.4, 1.5,
0.2)`, `maxit=200`, and jointly estimated covariance. A failed or warned base
fit makes that fixture inconclusive and is retained as such.

The authorised Loa loa teaching fit is an additional diagnostic case, not a
known-truth recovery case; its exact source version, file hash, projection and
model call are recorded in the checked example. It is not used to select a
policy or claim recovery. For reproducibility, it uses `lsda` 0.2.0's
`loaloa.rda` (SHA-256
`0acc9426c9842f93c264914a86e62641c160bdcdf134391cc57be7facb3c4811`), the
projection and binomial formula in the checked example, and the same 0.5
smoothness / `(0.5, 150, 0.5)` starting values. For this descriptive case
only, run five independent attempts per policy at 2,000 main draws, with seeds
`730001 + 10000*policy_index + attempt`. Report convergence and Monte Carlo
diagnostics without claiming operational policy selection or known-truth
recovery. Each synthetic-fixture policy uses five independent attempts at main
sample sizes 500, 1,000 and 2,000:

1. scale 1 without adaptation;
2. scales 0.75 and 1.5 without adaptation;
3. scale 1 with 500-draw pilot, shrinkage 0.5 and ESS gate 0.10;
4. scales 0.75 and 1.5 with the same pilot rule.

Use common response data across policies and sample sizes. MCML attempt seeds
are `520001 + 100000*case_index + 10000*policy_index +
100*size_index + attempt`, where all indices are one-based in the listed case,
policy, size and attempt order. The shared-parameter objective is recomputed
for each retained main sample at the common base Laplace parameter vector,
using that attempt's complete normalized proposal density. It is reported as a
noisy importance-sampling estimate, not as an exact likelihood.

Primary outputs are ESS fraction, maximum normalized weight, relative
likelihood MCSE, objective value at a shared parameter vector, optimized
parameters, adaptation-used status and optimizer convergence. Proposal policy
is not selected from one best seed. A policy is considered operational for a
case only if at least four of five attempts meet ESS fraction 0.10 and relative
MCSE 0.05 at 2,000 main draws. Otherwise the case remains inconclusive and
requires proposal redesign; increasing draws alone is not called a remedy.

The completed synthetic and Loa loa diagnostic results are in
`05_MCML_PROPOSAL_RESULTS.md`. No policy met the 4-of-5 operational criterion
for either synthetic fixture. All four Loa policies had warnings on every
attempt and are descriptive only. Proposal redesign and a separately
predeclared follow-up study remain necessary; no policy is selected.

## NNP: nearest-predecessor prediction

Exact ordinary joint prediction is the internal mathematical baseline on
manageable grids of 25, 100 and 225 prediction rows. The fitted models cover
Gaussian, binomial-logit and Poisson-log responses, one space-time case and one
grouped random-slope case. For each model, neighbor budgets are 5, 10, 20 and
40 (capped at \(M-1\)); input and maximin orders are both run with 20,000
draws and seeds 430001--430005.

The study records:

- relative Frobenius covariance error on the link scale;
- largest absolute marginal-variance error;
- absolute error in grouped target mean, SD and threshold probability;
- combined Monte Carlo standard error for draw summaries;
- elapsed time, full point-draw size, streamed target size and peak active
  point vectors; and
- sensitivity between input and maximin orders.

For a candidate neighbor budget to be recommended in a worked example, every
predeclared grouped target must have relative SD error at most 0.05 and absolute
threshold-probability error at most 0.02, with the latter difference's 95%
Monte Carlo interval wholly inside \([-0.02,0.02]\). Failure does not invalidate
the algorithm; it means that budget/order is not adequate for that case.
The all-predecessor same-seed equality test remains the exact boundary.

### Frozen fixture and target details

These details complete the NNP preregistration before its execution. For the
three point-family cases, generate 36 training sites on the Cartesian grid
`0:5` by `0:5`, and use prediction grids of side 5, 10 and 15 over `[0,5]^2`
(25, 100 and 225 rows, east varying fastest). Let `z` be the centered,
unit-standard-deviation version of `east + 0.5*north` on the training set;
apply the training mean and standard deviation unchanged to prediction rows.
Use fixed covariance `(spatial variance, scale, nugget) = (0.8, 2, 0.2)`,
exponential spatial correlation, and fixed-effect truth `(0.3, 0.4)`. For
each of the five response seeds, generate one Gaussian latent vector from
`0.8*exp(-distance/2) + 0.2*I`, add the fixed mean, then draw Gaussian
observations with the nugget represented in that covariance, binomial
successes with 10 trials and inverse-logit mean, or Poisson counts with
exponential mean. Fit with those covariance values fixed; the response-seed
order is 430001 through 430005. The fitted object and response are reused for
all prediction-grid sizes and approximation settings for that seed.

The space-time case has 16 locations on `0:3` by `0:3`, observed at times
`0,1,2`; the prediction grids use the same three side lengths and spatial
domain, with time cycling `0,1,2` in row order. Its covariance is
`0.8*exp(-distance/2)*exp(-time_lag/1.5) + 0.2*I`, with time measured in
years, fixed effects `(0.3, 0.4)`, and the same standardized spatial
covariate. The grouped-slope case uses the 36-site grid split into nine
contiguous 2-by-2 groups, a centered standardized slope `z`, and fixed
grouped covariance `matrix(c(0.3,0.05,0.05,0.08), 2)`. Its total training
covariance adds that grouped component to the same fixed spatial-plus-nugget
covariance. Prediction rows retain the group of their nearest training site
(distance ties resolved by first training row). Both additional cases use
Gaussian responses and the same five response seeds and fixed-effect truth.

For every case/grid/seed, the grouped targets are the equally weighted
overall prediction mean and four fixed spatial quadrants (west/east crossed
with south/north; points on the midlines go east/north). The threshold is
zero on the Gaussian link scale, 0.5 on the binomial response scale, and 1
on the Poisson response scale. Gaussian targets use link-scale draws;
binomial and Poisson targets use response-scale draws. Report each target's
mean and SD difference and threshold-probability difference between the
approximate and ordinary draws. The Monte Carlo interval is the normal 95%
interval for the difference, using the independent-seed replicate estimates;
also retain within-run draw MCSE. These fixed thresholds are diagnostic
cut-points, not domain-specific clinical claims.

For covariance comparisons, the denominator is the exact fitted predictive
link covariance Frobenius norm. The candidate criteria above apply separately
to all four quadrant targets and the overall target, and both orderings must
pass for a budget to be called adequate. The Monte Carlo decision is
inconclusive if the 95% interval is not wholly inside the stated margin.
Report unsuccessful fits or predictions without dropping their configurations.

## NNL: Gaussian nearest-predecessor likelihood

Generate independent Gaussian point designs at `n = 100, 225, 400` under
short-, medium- and long-range exponential correlation. Each design has an
intercept 0.3, one standard-normal covariate with coefficient 0.4, spatial
variance 0.8 and nugget variance 0.2; smoothness is 0.5. Training and 50 held-
out prediction coordinates are independently uniform on the unit square.
Scales are respectively 0.05, 0.20 and 0.50 times the square's diameter.
Generate the latent Gaussian field jointly over training and held-out sites,
then add independent nugget noise. Use response seeds
`460001 + 10000 * scenario_index + replicate`, with size-major, range-minor
scenario order and zero-based scenario indices; use replicate indices 1--25.
For each size/range setting, fit the declared dense baseline once and evaluate
its 50 held-out conditional means.

For input and maximin orderings, fit budgets 5, 10, 20 and 40, capped at
`n - 1`, plus the all-predecessor boundary. Covariance parameters stay fixed
at their known generating values to isolate the likelihood approximation.
Record regression estimates, approximate log likelihood, held-out prediction
RMSE, elapsed time, retained object size and externally measured peak resident
memory. Dense ML is the internal manageable-size baseline.

The all-predecessor fit must match dense fixed-parameter likelihood,
coefficients, conditional coefficient covariance and held-out means/standard
errors within `2e-10`. A candidate budget is called stable for a scenario only when
at least 90% of usable replicates have every regression difference below 0.05
standard errors, every positive covariance difference below 10% on the log
scale, and out-of-sample point prediction RMSE within 2% of dense ML. The
smallest passing budget may be recommended only if both orders pass. Fewer than
90% usable fits is inconclusive. Runtime or memory results from one platform
are reported, never generalized.

## Execution and reporting order

1. Build and checksum one source archive; do not change code during a run.
2. Run GRS and TSC analytic/simulation cases sequentially first.
3. Run MCP attempts with all seeds retained, without selecting a policy mid-run.
4. Run NNP and NNL exact baselines before approximate configurations.
5. Produce machine-readable raw tables and a readable discrepancy report.
6. Mark every case `pass`, `fail` or `inconclusive` against the rules above.
7. If code changes are needed, assign a new archive checksum and rerun every
   affected case; do not splice favorable results across versions.

