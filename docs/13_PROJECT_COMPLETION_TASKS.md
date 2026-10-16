# SPmodels project-completion task list

## How to read this list

This is the working completion checklist for the full initial feature inventory.
`Done` means implemented with an independent derivation or algorithm, local
tests and a clean package check recorded in the implementation register.
`Active` is the current development block. `Pending` remains required.
`Owner/external decision` cannot be completed truthfully without information or
authority outside this package directory.

## Verified foundation

- [x] Matérn correlation and first/second scale derivatives.
- [x] Gaussian spatial and independent maximum-likelihood fitting.
- [x] Projected point geometry, prediction and regular grids.
- [x] Gaussian simulation, variograms and residual diagnostics.
- [x] Binomial-logit and Poisson-log Laplace fitting.
- [x] Latent modes, conditional MALA and basic sampler diagnostics.
- [x] Importance-sampling MCML correction and likelihood diagnostics.
- [x] Point and joint prediction, parameter draws and response transformations.
- [x] Grid/polygon-centre targets, weighting and exceedance summaries.
- [x] Fixed-fold cross-validation and predictive-draw assessment.

## Completed current block

- [x] F11/F29: fitted non-Gaussian response simulation.
- [x] M07: fitted non-Gaussian residual diagnostics.
- [x] Empirical Pearson, deviance or seeded randomized-quantile residual
  variograms for binomial and Poisson fits, without an invented fitted curve.
- [x] Non-Gaussian residual permutation envelopes and global curve tests with
  explicit exchangeability and fitted-parameter limitations.
- [x] Non-Gaussian fitted-model parametric residual envelopes that repeat the
  declared fitting route and retain convergence/warning evidence.
- [x] F02/F05: geometry-driven folds that keep duplicate locations together.
- [x] Buffered spatial validation folds with fold-specific near-test training
  exclusion and recorded remaining sample sizes.
- [x] Regular spatial-block folds with explicit widths/origin, whole-block
  allocation, deterministic or locally seeded assignment and retained audit
  metadata.
- [x] Environmental-space folds with optional standardization, duplicate-
  location consistency, local seed handling and retained centres/scaling.
- [x] Complete derivations, tests and help for this block.
- [x] Complete final installation and PDF refresh
  for this block.

## Remaining statistical and numerical development

- [x] F32 foundation: implement and independently derive a positive-definite
  non-stationary Matérn correlation primitive with local range/anisotropy
  matrices and stationary-limit checks. The full estimation and fitted-model
  workflow remains pending.
- [ ] Specify local-matrix parameterization/regularization, identifiability,
  and validation designs for non-stationary estimation before implementing
  fitted-model support.

- [x] F01/F04/F28: multiple-chain controls, between-chain diagnostics and
  autocorrelation-aware Monte Carlo standard errors.
- [x] F03/M05/M08: repeated simulation/refit orchestration and parameter-
  recovery summaries for one or multiple named fitted scenarios.
- [x] Optional deterministic parallel dispatch across named recovery scenarios
  on non-Windows systems, with sequential equivalence tests.
- [x] Optional deterministic parallel dispatch across recovery replicates on
  non-Windows systems, with explicit MCML seeds and exact equivalence tests.
- [x] Precision-directed simulation recovery with bias-MCSE targets, minimum
  successful refits, geometric escalation and an exact auditable cap.
- [x] Cartesian parameter-grid generation over recoverable fitted parameters,
  with rebuilt Gaussian/non-Gaussian generative states and retained design.
- [x] F17/F20/F21: predictive interval calibration with finite-location bounds,
  plus validation/simulation display methods over tested underlying values.
- [x] M02/M03: point/grid and polygon target geometry joins and display method
  built from tested unrounded summaries.
- [x] F09: explicit distance-summary utility with duplicate-location policy.
- [x] F25: projected CRS/zone helper with boundary and hemisphere contracts.
- [x] F24: exact fractional cell-area polygon weights and explicit unmatched/
  empty-region policies.
- [x] F30/M04/M06 foundation: stable unrounded model, prediction and validation
  tables with direct schema/value tests.
- [x] F30/M04/M06 presentation: fixed-precision formatting, plots and assembled
  dependency-free reports built as views over retained unrounded tables.
- [x] Local transformed Wald confidence intervals for estimated coefficients,
  variances, spatial/temporal scales and intervention parameters.
- [x] F31: fitted-model prediction/design updating with explicit preservation
  rules, exact row matching and in-formula offset reconstruction.
- [x] F26: additive Gaussian random-intercept grouping, multiple variance
  components, known-level prediction, simulation and reporting.
- [x] Additive non-Gaussian latent random intercepts with fixed or estimated
  variances and known-level prediction.
- [x] Explicit marginal prediction policy for unseen grouping levels, with
  shared-new-level covariance and a safe error default.
- [x] Fixed-covariance grouped random slopes and correlated intercept-slope
  effects, including known/unseen-level prediction and all fitted workflows.
- [x] Optional constrained estimation of grouped random-slope covariance
  matrices using a positive-definite Cholesky parameterization, minimum-design
  checks, natural-scale uncertainty and complete Gaussian/Laplace/MCML flow.
- [x] F08/W01: validated decay-adjusted intervention
  histories propagated through binomial fitting, MCML, prediction, simulation,
  validation and reporting.
- [x] F08/F18: constrained Laplace estimation, optional normalized alpha
  penalty, delta-method uncertainty and intervention-effect display.
- [x] Intervention extension: MCML re-estimation of impact parameters with
  transformed curvature and optional alpha penalty.
- [x] Spatio-temporal foundation: separable Matérn-exponential correlation
  primitive with limiting-case and positive-semidefinite checks.
- [x] Wire temporal scale and survey times through fitting, prediction,
  simulation, diagnostics and uncertainty.
- [x] Joint positive log-scale temporal estimation for Gaussian ML and
  binomial/Poisson Laplace fitting, including MCML re-estimation.
- [x] Robust MCML foundation: independent replicates, proposal refresh between
  attempts, precision targets and automatic sample-size escalation.
- [x] Normalized multiple-scale deterministic-mixture MCML proposals with
  balanced allocation, exact mixture density and stratified likelihood MCSE.
- [x] Parameter-level Monte Carlo covariance from independent MCML replicates,
  kept separate from and optionally combined with fitted-Hessian covariance.
- [x] Opt-in within-attempt MCML proposal adaptation with an independent
  discarded pilot, regularized weighted moments, exact normalized retained
  proposal and separate pilot/main diagnostics.
- [x] Group-constrained folds with whole-unit retention, greedy row balancing,
  randomized allocation and retained audit metadata; arbitrary supplied fold
  identifiers remain supported directly.
- [x] Exact bounded-memory point prediction for large grids with results
  checked against ordinary Gaussian and non-Gaussian prediction.
- [x] Bounded-memory nearest-predecessor approximate joint-grid draws with an
  exact all-predecessor boundary, explicit order dependence and retained local
  graph/variance metadata; point chunks are not treated as joint covariance.
- [x] Direct bounded-memory grouped-target streaming over the predecessor
  graph, including sparse polygon allocations and recorded peak active storage.
- [x] Analytic implied-covariance diagnostics for selecting neighbor budgets
  on manageable subsets without Monte Carlo error.
- [x] Explicit opt-in Gaussian nearest-predecessor training-likelihood
  foundation with bounded covariance blocks, local prediction, auditable
  graph/order metadata and an exact all-predecessor boundary.
- [x] Execute and archive the predeclared multi-scenario neighbor-likelihood
  statistical accuracy/stability study and exact boundary results; no general
  production budget is recommended because medium/long-range cases at n =
  225/400 require every predecessor to meet the registered criteria.
- [x] Record NNL per-fit elapsed time and externally measured whole-run peak
  RSS in a full seeded rerun; verify all deterministic accuracy outcomes match
  the archived study.
- [ ] Measure externally observed peak RSS by NNL configuration and on
  additional supported platforms before drawing comparative memory conclusions.
- [x] Execute and archive the preregistered NNP joint-prediction study across
  all five model cases, grid sizes, orderings, budgets and seeds; retain failed
  configurations and make no budget recommendation without the registered
  grouped-target criteria passing. Results: `validation/06_NEIGHBOR_PREDICTION_RESULTS.md`;
  no universal budget passed because the large space-time case was inconclusive.
- [ ] Optimize the pure-R graph/local-factor route only after broader profiling
  on supported platforms identifies a defensible target.
- [ ] Derive and validate separate scalable likelihoods for non-Gaussian,
  grouped or temporal models before exposing any such extension.

## Validation and adoption gates

- [x] Expand the consumer audit into exact call schemas and minimal migration
  tests for MBGapp, lsda and the course; keep those projects read-only until the
  owner requests integration edits.
- [x] Predeclare configuration-specific deterministic and stochastic criteria
  and the initial neutral case matrix.
- [x] Predeclare recovery and approximation studies for fitted grouped
  covariance, non-Gaussian temporal scale, MCML proposals and bounded-memory
  joint prediction, including fixed seeds, escalation and failure rules.
- [x] Execute the registered GRS grouped-slope study and report every case,
  including inconclusive bias-MCSE results, systematic scale bias, Hessian
  uncertainty availability and local-Wald coverage limitations.
- [x] Diagnose spatial-range identifiability and interval uncertainty under
  the registered GRS scenarios without relaxing their criteria; existing
  Hessian uncertainty is broad and poorly calibrated, and this is not a formal
  profile-likelihood analysis.
- [x] Execute the registered TSC temporal-scale and MCML study through its
  attempt caps, preserve all fit/repeat rows, validate seeds and gates, and
  report inconclusive usability, temporal-scale recovery and MCML diagnostics.
- [x] Execute all predeclared MCML proposal policies at all registered main
  sample sizes on both synthetic cases, plus the five-attempt-per-policy Loa
  loa diagnostic; preserve common-parameter objectives and weight diagnostics.
  No synthetic policy met the 4-of-5 operational gate.
- [ ] Redesign MCML proposals based on retained ESS, maximum-weight and MCSE
  evidence, then predeclare a new independent validation study before fitting.
- [ ] Freeze the separately managed reference-validation environment.
- [ ] Run and record every required reference comparison; unresolved or
  inconclusive cases do not count as passes.
- [x] Add a checked binomial teaching example using suitable lsda data with
  version, file hash, provenance, projection and a recorded successful fit.
- [ ] Add checked count, grid-target and intervention teaching examples using
  suitable lsda data. The audited lsda count dataset is polygon-supported and
  the Galicia grid metadata does not yet establish an exact coordinate join,
  so neither is being forced into an invalid point workflow.
- [ ] Run supported cross-platform checks, including macOS, and record platform
  limitations.
- [x] Record a reproducible first small/moderate dense-workflow timing baseline
  before optimization.
- [x] Profile the measured moderate-fit bottleneck and optimize the exact
  half-integer Matérn cases, with a controlled same-process kernel A/B check.
- [ ] Repeat broader controlled cross-platform timings before selecting a
  scalable approximation or compiled kernel.

## Release gates requiring owner or external decisions

- [ ] Owner/external decision: final public package name.
- [ ] Owner/external decision: distribution licence.
- [ ] After identity/licence decisions, resolve CRAN release notes by setting
  a release version and separating package-shipped user documentation from
  top-level project-management records. The current `--as-cran` run also noted
  that HTML manual validation was skipped because `tidy` is unavailable.
- [ ] Owner/external decision: authorize edits and migration tests in consumer
  projects.
- [x] Create a private, versioned GitHub source repository after owner
  authorization; public publication remains pending.
- [ ] Prepare release metadata, citation, website and CRAN/R-universe route.
- [ ] Optional release presentation: journal-specific typesetting and
  HTML/office-document exporters.
- [ ] Submit to external services only when explicitly authorized.

## Completion rule

Do not describe SPmodels as complete until every applicable item above and each
row of `02_FEATURE_INVENTORY.md` has recorded implementation, derivation,
validation and consumer evidence. A clean `R CMD check` is necessary but is not
by itself project completion.
