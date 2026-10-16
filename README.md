# SPmodels

SPmodels is Olatunji Johnson's independently developed R package for spatial
modelling in research and teaching. It is growing into a broad model collection
with a consistent formula-based interface, independently implemented methods
and full mathematical derivations accessible to MSc students. Suitable
datasets from the owner's lsda package support checked worked examples.

The first model collection covers Gaussian, binomial and Poisson geostatistical regression, simulation, prediction, diagnostics and intervention models. Future additions will follow research and curriculum needs. MBGapp, lsda and the course are the initial adoption targets; SDALGCP2 provides the preferred interface style.

**Status: implementation in progress.** The implemented subset includes dense
Gaussian, binomial-logit and Poisson-log point geostatistical regression;
formula processing; standard R methods; projected point `sf` input; conditional
prediction; fitted response simulation; variograms; and family-appropriate
residual diagnostics. Important validation, scalability and release gates
remain open. Deterministic Matérn scale derivatives,
binomial/Poisson latent-mode calculation, conditional MALA sampling and basic
multiple-chain diagnostics with pooled ESS, MCSE and rank-normalized split
R-hat are implemented. `sp_fit()` also fits binomial-logit
and Poisson-log spatial models by a clearly labelled Laplace approximation;
point prediction integrates through their nonlinear links. `sp_mcml()` can
refine those fits using a normalized importance-sampling estimate of the
marginal likelihood with explicit Monte Carlo diagnostics. Independent MCML
replicates, parameter-level Monte Carlo covariance, normalized multiple-scale
proposals, optional discarded-pilot proposal adaptation and precision-directed
sample-size escalation are available through `sp_mcml_control()`,
`sp_mcml_replicates()` and `sp_mcml_precision()`.

The registered grouped-covariance and non-Gaussian temporal-scale recovery
studies are now archived. They exposed weak spatial/temporal range recovery,
poor local interval behavior and unstable MCML diagnostics in the tested
designs; these studies do not constitute reference agreement. See the
case-specific results in `docs/validation/` before using those fitted
extensions scientifically.
The registered MCML proposal study likewise found no operational policy in
its synthetic cases; no proposal setting is recommended.

Use this `Research/SPmodels` directory as the working master on each machine.
Machine-specific absolute paths may differ between Linux and macOS; synced copies
should not be treated as separate project identities. The Codex outputs copy is
only a delivery snapshot.

## Read the planning documents

| Document | Purpose |
| --- | --- |
| [Project charter](docs/01_PROJECT_CHARTER.md) | Research and teaching aims, scope and acceptance criteria |
| [Feature inventory](docs/02_FEATURE_INVENTORY.md) | Planned modelling and supporting capabilities |
| [Proposed interface](docs/03_API_DESIGN.md) | Formula, data, family and standard R methods |
| [Validation protocol](docs/04_VALIDATION_PROTOCOL.md) | Independent checks and numerical reference comparisons |
| [Adoption plan](docs/05_MIGRATION_PLAN.md) | MBGapp, lsda and course transition |
| [Roadmap and decisions](docs/06_ROADMAP_AND_DECISIONS.md) | Milestones and a restart prompt |
| [Reference validation preparation](docs/07_REFERENCE_BASELINE.md) | Reproducibility requirements for future comparisons |
| [Recovery/approximation design](docs/validation/02_RECOVERY_AND_APPROXIMATION_DESIGN.md) | Predeclared grouped, temporal, MCML and neighbor studies |
| [Grouped-slope recovery results](docs/validation/03_GROUPED_SLOPE_RECOVERY_RESULTS.md) | Synthetic covariance recovery, mathematical checks and spatial-range diagnostics |
| [Temporal-scale recovery results](docs/validation/04_TEMPORAL_SCALE_RECOVERY_RESULTS.md) | Non-Gaussian temporal recovery, boundary behavior and MCML repeat diagnostics |
| [MCML proposal-policy results](docs/validation/05_MCML_PROPOSAL_RESULTS.md) | Independent-repeat, ESS, weight and adaptation outcomes; no policy passed |
| [Teaching and future models](docs/08_TEACHING_AND_FUTURE_MODELS.md) | Curriculum needs and model expansion |
| [Distribution and CRAN](docs/09_DISTRIBUTION_AND_CRAN.md) | Release, installation and maintenance |
| [Derivations and teaching data](docs/10_DERIVATIONS_AND_TEACHING_DATA.md) | Complete MSc-level derivations for every method and selected lsda datasets |
| [Implemented package guide](docs/SPmodels_Implemented_Guide.pdf) | Reader-focused guide to everything implemented, its mathematics, limitations and next steps |
| [Technical derivations PDF](docs/SPmodels_Technical_Derivations.pdf) | Compact compilation of the equation-level derivation records |
| [Implementation register](docs/11_IMPLEMENTATION_REGISTER.md) | Implemented slices, evidence and unresolved release metadata |
| [Initial consumer audit](docs/12_CONSUMER_AUDIT.md) | Available-source audit and retained downstream contract |
| [Galicia Gaussian example](docs/examples/01_GALICIA_GAUSSIAN.md) | Checked lsda teaching workflow with explicit projection and provenance boundary |
| [Loa loa binomial example](docs/examples/02_LOALOA_BINOMIAL.md) | Checked prevalence workflow with projected distances, diagnostics and spatial validation |
| [lsda data suitability audit](docs/examples/00_LSDA_DATA_AUDIT.md) | Read-only support/provenance decisions for the remaining teaching examples |
| [Dense-workflow benchmark](docs/benchmarks/01_DENSE_WORKFLOWS.md) | Reproducible first timing baseline and measured scaling boundary |
| [MCML proposal benchmark](docs/benchmarks/02_MCML_PROPOSALS.md) | Real-data proposal diagnostics and limits of scale-mixture tuning |
| [Neighbor-prediction diagnostic](docs/benchmarks/03_NEIGHBOR_PREDICTION.md) | Synthetic joint-covariance error and ordering sensitivity |
| [Neighbor-likelihood diagnostic](docs/benchmarks/04_NEIGHBOR_LIKELIHOOD.md) | Gaussian training-likelihood exact boundary, approximation differences and timing limits |
| [Instructions for future work](AGENTS.md) | Planning boundary and project conventions |

## Intended user experience

The first implemented fitting path uses ordinary data frames with two numeric
Euclidean coordinate columns and explicit common distance units, or projected
point `sf` data:

```r
fit <- sp_fit(
  outcome ~ rainfall + elevation,
  data = survey_points,
  coords = ~ east_km + north_km,
  control = sp_control(distance_units = "km")
)
summary(fit)
predict(fit, newdata = prediction_points, se.fit = TRUE,
        target = "observation")
```

Gaussian identity, binomial logit and Poisson log responses support point and
joint prediction, fitted response simulation and diagnostics. Geographic
longitude/latitude geometry is rejected until explicitly projected. Random
intercepts, fixed or constrained-estimated grouped random slopes, and separable space-time
dependence are available within their documented boundaries; see the guide for
unsupported structures.

## Distribution

Plan for CRAN submission and straightforward installation, with R-universe as
a possible complementary or alternative route. SPmodels is a working name;
the public identity must be reviewed before submission.

## Current implementation

This folder contains a development package with Matérn operations, Gaussian
fitting and prediction, and binomial/Poisson Laplace and MCML fitting.
Gaussian models also have an explicit opt-in nearest-predecessor training
likelihood through `sp_fit_neighbor()`. It profiles the linear mean from local
conditional densities, retains its order and graph, and becomes exactly the
dense Gaussian likelihood when all predecessors are kept. The initial route
does not support non-Gaussian, grouped or temporal covariance models and is a
bounded-covariance-memory foundation, not a speed guarantee.
Local analytic tests, projected geometry, Gaussian simulation, variograms,
diagnostics (including seeded randomized quantile residuals for discrete
responses and their empirical residual variograms), conditional latent
sampling, joint predictive draws, projected
grids, polygon-region allocation and draw-based exceedance or grouped targets
are implemented and tested. Predictive draws can optionally include a local
Hessian approximation to fitted-parameter uncertainty. Fixed-fold model
refitting and draw-based recovery, coverage and CRPS assessment are also
implemented. Cross-validation repeats MCML inside every training fold when the
supplied model was MCML-corrected. Geometry-driven coordinate folds keep
repeated locations together.
Large point sets can use `sp_predict_chunked()` for exact bounded-memory
pointwise means and marginal standard errors; joint grid dependence still
uses either the ordinary exact joint-draw route or the explicitly approximate
bounded-memory nearest-predecessor route in `sp_predict_draws_neighbor()`.
Grouped grid targets can be accumulated directly with
`sp_predict_target_neighbor()`, which releases point-draw vectors after their
last graph use instead of automatically retaining every point draw.
For manageable validation subsets, `sp_neighbor_covariance()` reconstructs
the approximation analytically and reports exact fitted-covariance errors.
Buffered folds additionally exclude training rows within a declared fitted-unit
distance of each held-out location and record each fold's remaining sample size.
Regular block folds with declared axis widths and origin keep whole occupied
blocks together and support deterministic or locally seeded balanced allocation.
Environmental folds cluster unique locations by declared numeric covariates,
with optional standardization and repeated-coordinate consistency checks.
Group-constrained folds keep every level of a declared fitted-data unit wholly
within one fold while greedily balancing held-out row counts.
Predictive interval calibration curves include pointwise Wilson bounds, and
cross-validation, predictive assessment and parameter recovery have display
methods over their stored unrounded values.
Parameter-recovery designs can be generated as explicit Cartesian grids over
recoverable fitted parameters and dispatched sequentially or across scenario
or within-scenario replicate workers with deterministic seeds.
Held-out binomial and Poisson log scores integrate response probability over
the approximate Gaussian predictive link distribution rather than substituting
only its marginal response mean.
Matched prediction scenarios can update selected predictor or in-formula
offset values while explicitly preserving row order, locations, geometry,
unmentioned columns and the fitted parameter state.
Repeated fitted-model simulation and refitting reports parameter bias, RMSE,
empirical spread, reported standard errors and Wald interval coverage.
Named multi-scenario orchestration retains every individual recovery result,
uses reproducible scenario seeds and produces an unpooled combined summary.
Named scenarios can optionally be dispatched in parallel on non-Windows
systems without changing their preassigned random streams.
`sp_simulation_precision()` can instead increase a single scenario's complete
recovery experiment until every bias-MCSE target passes or an exact cap is
recorded.
Gaussian, binomial and Poisson spatial fits support one or more additive latent
random-intercept grouping variables, with fixed or fitted variance components,
known-level conditional prediction, opt-in marginal prediction for unseen
levels, and joint simulation; Gaussian fits also provide group-aware fitted
variograms.
Grouped slope-only or correlated intercept-slope components may likewise use
fixed covariance matrices or opt-in positive-definite Cholesky estimation,
with natural-scale uncertainty and explicit minimum-design checks.
Binomial spatial fits can use a decay-adjusted intervention history with an
explicit same-time survey-before-rollout rule; its fixed or constrained fitted parameters and multipliers
are carried through Laplace fitting, MCML, prediction, simulation and
validation. A delta-method intervention-effect curve and display are available.
Both impact and decay can be re-estimated by the normalized MCML route with
their optional alpha penalty retained, and joint predictive parameter draws
rebuild the intervention multiplier for every draw.
Formal permutation or fitted-model residual-variogram envelopes, auditable
pair-distance summaries, explicit UTM CRS selection, exact fractional-area
polygon targets and spatial target geometry joins are also implemented.
Permutation envelopes accept supported discrete-family residuals; fitted-model
parametric envelopes for discrete families repeat the declared fitting route
and retain every refit's convergence code and warnings.
Stable model, prediction and validation tables are available through
`sp_model_table()`, `sp_prediction_table()` and `sp_validation_table()` for
downstream applications and teaching materials.
`confint()` supplies local identity-, log- or logit-scale Wald intervals that
respect parameter support and retain missing limits for unusable curvature.
`sp_format_table()` and `sp_report()` add readable presentation views and an
assembled analysis object while retaining the exact unrounded source tables.
Separable Matérn-exponential space-time models are available through
`spacetime_cor()` and `sp_temporal()`. Survey times and a fixed or jointly
estimated temporal scale are carried through Gaussian, binomial and Poisson
fitting, prediction, simulation, variograms and validation; non-Gaussian MCML
also re-estimates the declared scale.
No numerical reference comparison or consumer migration has been completed.
A private, versioned GitHub source repository now holds the working master;
there has been no public release or external submission. Future numerical
claims require the evidence specified in the validation protocol. The
implementation register records the latest package-check evidence. No release
identity or numerical comparison result has been invented.
