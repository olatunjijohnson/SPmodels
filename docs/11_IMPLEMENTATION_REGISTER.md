# SPmodels implementation register

## Status

Implementation began on 12 October 2026 after the owner's explicit instruction.
This register distinguishes implemented evidence from planned scope. The full
feature inventory in `02_FEATURE_INVENTORY.md` remains in scope.

## First implementation slice

| Field | F15 record |
| --- | --- |
| Capability | Matérn correlation |
| Public API | `matern_cor(distance, scale, smoothness)` |
| Parameter convention | `x = distance / scale`; scale is not labelled a practical range |
| Supported inputs | Finite non-negative numeric distances; one positive scale and smoothness |
| Output | Correlations with input dimensions and names preserved |
| Failure conditions | Non-numeric, missing, infinite or negative distances; invalid scalar parameters; unrecoverable numerical evaluation |
| Owning module | `R/matern.R` |
| Derivation | `docs/derivations/01_MATERN_CORRELATION.md`, equations M-01--M-09 |
| MSc prerequisites | Introductory probability, calculus and linear algebra |
| Synthetic checks | Zero limit, half-integer closed forms, scale invariance, PSD matrix, invalid inputs |
| Provisional deterministic criterion | `1e-12` for analytic well-conditioned cases; PSD eigenvalue floor `-1e-12` |
| Reference comparison | Planned in a separate pinned environment; not yet run |
| Consumers | Read-only call-site audit completed in `docs/12_CONSUMER_AUDIT.md`; authorized integration edits remain outstanding |
| Status | Implemented locally; independent deterministic tests and package checks pass under R 4.6.1 on Linux |

## Gaussian geostatistical fitting slice

| Field | F10/F12/F22 record |
| --- | --- |
| Capability | Gaussian identity-link geostatistical regression, covariance fitting and conditional prediction |
| Public API | `sp_control()`, `sp_fit()`, `predict()` and standard methods |
| Mean model | R formula/model-frame design matrix, optional subset and offset, full-rank requirement; independent Gaussian ML route through `spatial = FALSE` |
| Spatial input | Data frame plus two finite numeric Euclidean coordinates and explicit units, or non-empty projected point `sf` geometry with a known CRS |
| Covariance | Stationary isotropic Matérn spatial covariance plus independent nugget; fixed smoothness |
| Estimation | Dense ordinary maximum likelihood with regression coefficients profiled by generalized least squares; log-scale bounded optimisation |
| Prediction | Universal-kriging conditional mean; latent or future-observation standard error |
| Uncertainty boundary | Point-prediction standard errors include mean-estimation correction conditional on covariance parameters; covariance standard errors use observed profile curvature and the delta method; joint predictive draws optionally add the local Hessian approximation recorded below |
| Owned class | `spmodels_fit` with `coef`, `vcov`, `fitted`, `residuals`, `logLik`, `summary`, `print` and `predict` methods |
| Owning modules | `R/control.R`, `R/fit.R`, `R/methods.R` |
| Derivation | `docs/derivations/02_GAUSSIAN_GEOSTATISTICAL_REGRESSION.md`, equations G-01--G-17 and G-12a--G-12b |
| Synthetic checks | Direct GLS, likelihood and conditional Gaussian equations; row permutation; latent/observation variance; contract failures |
| Provisional deterministic criterion | `1e-11` for the small well-conditioned analytic fixture |
| Reference comparison | Planned in a separate pinned environment; not yet run |
| Consumer audit | Exact read-only call/schema audit completed; see `docs/12_CONSUMER_AUDIT.md`; no external project was changed |
| Status | Implemented locally; package-level verification recorded below |

## Gaussian nearest-predecessor likelihood slice

| Field | F10/F12 scalable extension record |
| --- | --- |
| Capability | Explicit order-dependent Gaussian nearest-predecessor likelihood with profiled mean and local point prediction |
| Public API | `sp_fit_neighbor()` with `coef`, `vcov`, `logLik`, `fitted`, `residuals`, `summary`, `print` and `predict` methods |
| Approximation | Product of Gaussian conditionals using at most the declared number of earlier nearest observations; input, deterministic maximin or supplied order is retained |
| Exact boundary | Retaining every predecessor reproduces the dense Gaussian log likelihood, regression estimate and conditional coefficient covariance to floating-point precision for any tested valid order |
| Storage boundary | Graph storage is `O(nm)` and local covariance workspace is `O(m^2)`; exact neighbor discovery remains `O(n^2)` time and the current R local factorizations are not claimed faster than dense BLAS |
| Prediction | Each new point conditions on its nearest training observations; latent/future-observation marginal uncertainty is returned without inventing a joint covariance |
| Supported model | Gaussian identity-link stationary isotropic Matérn plus nugget, with fixed smoothness; fixed or jointly estimated covariance parameters |
| Explicit exclusions | Non-Gaussian, grouped, random-slope and temporal covariance models; production-budget recommendation pending the predeclared NNL study |
| Owning module | `R/neighbor-fit.R` |
| Derivation | `docs/derivations/20_GAUSSIAN_NEAREST_PREDECESSOR_LIKELIHOOD.md`, equations NL-01--NL-15 |
| Synthetic checks | Dense all-predecessor equality in input/maximin orders, graph bound, explicit order, fixed/estimated fits, local prediction, offset reconstruction, RNG preservation and failure contracts |
| Diagnostic | `docs/benchmarks/04_NEIGHBOR_LIKELIHOOD.md`; at n=225 small budgets retained 188--232 kB versus 1,007 kB dense, but took 0.113--0.125 s versus 0.009 s dense in the warmed fixed-parameter run |
| Process memory | Separate fixed-parameter runs at n=2,000 used 130,496 KiB peak RSS and a 1.59 MB retained object for neighbor 20, versus 300,132 KiB and 65.22 MB dense; neighbor 20 also completed n=5,000 at 134,112 KiB peak RSS in 6.25 s |
| Status | Implemented and locally tested; broader recovery, peak-memory and cross-platform studies remain planned |

## Gaussian simulation and assessment slice

| Field | F11/F27 supporting record |
| --- | --- |
| Capability | Unconditional fitted-model simulation, binned residual variogram, fitted variogram, formal variogram envelopes and decorrelated residual diagnostics |
| Public API | `simulate()`, `sp_variogram()`, `sp_variogram_envelope()`, `sp_diagnostics()` and plot methods |
| Simulation target | Response vectors at observed sites, conditional on fitted parameters |
| Variogram | Equal-width distance bins, explicit pair counts, mean-model or conditional residuals |
| Diagnostics | Mean-model, conditional and Cholesky-decorrelated residuals |
| Owning module | `R/assessment.R` |
| Derivation | `docs/derivations/03_SIMULATION_VARIOGRAM_AND_DIAGNOSTICS.md`, equations A-01--A-09, and `docs/derivations/11_SPATIAL_AUDIT_ENVELOPES_AND_AREA_TARGETS.md`, equations D-05--D-10 |
| Synthetic checks | Exact seeded Cholesky transform, complete pair accounting, direct whitening solve, projected-geometry contracts, and direct reconstruction of every seeded permutation curve and the global Monte Carlo p-value |
| Not yet included | Additional formal residual tests beyond the implemented variogram envelope |
| Status | Implemented locally; clean package-check evidence recorded below |

## Non-Gaussian deterministic foundation

| Field | F13/F14/F16 record |
| --- | --- |
| Capability | First and second Matérn scale derivatives; binomial-logit and Poisson-log latent Gaussian modes and conditional single/multiple-chain sampling |
| Public API | `matern_scale_derivative()`, `sp_latent_mode()`, `sp_mcmc_control()`, `sp_latent_sample()`, `sp_sampler_diagnostics()` and its plot method |
| Estimation boundary | Supplied fixed predictor and covariance only; no non-Gaussian parameter estimator is claimed |
| Algorithm | Analytic Bessel derivatives; damped Newton solve; fixed-preconditioner MALA with asymmetric proposal correction; dispersed starts; rank-normalized split/folded R-hat; pooled initial-positive-pair ESS and MCSE |
| Derivation | `docs/derivations/04_MATERN_DERIVATIVES_AND_LATENT_MODE.md`, equations N-01--N-36 |
| Synthetic checks | Exponential closed forms, centred finite differences, direct score/curvature, permutation invariance, seeded one/four-chain reproducibility, one-dimensional posterior moments, chain identity/selection, ESS/MCSE/R-hat structure and invalid inputs |
| Status | Implemented locally; deterministic tests and package check recorded below |

## Non-Gaussian Laplace fitting slice

| Field | F13/F14 record |
| --- | --- |
| Capability | Point-supported binomial-logit and Poisson-log latent Gaussian geostatistical regression |
| Public API | `sp_fit()` plus `coef()`, `vcov()`, `fitted()`, `residuals()`, `logLik()`, `summary()` and `print()` |
| Response conventions | Bernoulli values; successes/failures matrix; proportions with integer trial weights; Poisson integer counts with exposure in a log offset |
| Covariance | Matérn spatial covariance plus independent latent variance; fixed smoothness |
| Estimation | Joint maximisation of a Laplace-approximated marginal log likelihood; fixed-covariance fitting is also supported |
| Approximation boundary | This is explicitly not MCML; fitted values are conditional response means at the latent mode |
| Owning module | `R/nongaussian-fit.R` with shared mode and covariance operations |
| Derivation | `docs/derivations/05_NONGAUSSIAN_LAPLACE_FITTING.md`, equations L-01--L-16 |
| Synthetic checks | Direct objective reconstruction, equivalent binomial encodings, row permutation, Poisson offset, response ranges and invalid inputs |
| Prediction | Approximate link distribution at new points; Poisson lognormal moments; binomial logistic-normal quadrature; latent and future-observation targets |
| Not yet included | An independent non-Gaussian limiting route; fitted simulation, residual diagnostics, MCML correction and optional Hessian parameter draws are recorded separately below |
| Status | Implemented locally; package-level verification recorded below |

## Monte Carlo maximum-likelihood correction slice

| Field | F13/F14 supporting record |
| --- | --- |
| Capability | Importance-sampling Monte Carlo marginal-likelihood refinement of binomial-logit and Poisson-log spatial fits |
| Public API | `sp_mcml()`, `sp_mcml_control()`, `sp_mcml_replicates()`, `sp_mcml_precision()`, ordinary fitted-model methods and downstream prediction/validation methods |
| Proposal | One normalized multivariate Gaussian, or a declared multiple-scale deterministic mixture with balanced allocation and exact density; optional adaptation uses an independent discarded pilot and regularized importance-weighted centre/covariance before a new main sample |
| Objective | Absolute marginal likelihood from normalized joint-density/proposal weights; fixed common random numbers during optimization |
| Diagnostics | Effective sample size and fraction, maximum normalized weight, single-proposal or stratified-mixture relative likelihood Monte Carlo error and ESS warning threshold |
| Parameterization | Regression coefficients on their ordinary scale; estimated covariance parameters on bounded log scales |
| Prediction boundary | Parameters are MCML-corrected; downstream conditional latent prediction remains a local Gaussian approximation |
| Cross-validation | Preliminary Laplace fitting and MCML correction are both repeated inside every training fold |
| Owning module | `R/mcml.R`, with existing latent-mode and covariance primitives |
| Robustness workflow | Independent fitted replicates retain every result, empirical parameter-level MC covariance, average fitted-Hessian sampling covariance and their explicit sum; geometric escalation stops only at declared relative-MCSE and ESS targets or the exact cap, with optional proposal refresh |
| Derivation | `docs/derivations/08_MONTE_CARLO_MAXIMUM_LIKELIHOOD.md`, equations C-01--C-16, including C-05a--C-05g, C-12a and C-13a--C-13b |
| Synthetic checks | Direct scalar integration, exact seeded reproduction, normalized original/adapted mixture densities and allocation, pilot/main separation, weighted-moment regularization, seed isolation, likelihood reconstruction, binomial/Poisson paths, prediction, fold-wise correction, direct MC covariance and first/capped precision attempts |
| Not yet included | Bridge sampling and sparse large-data methods |
| Status | Implemented, checked and installed; package-level verification recorded below |

## Joint prediction, grid and target slice

| Field | F06/F07/F22/F23/F24 supporting record |
| --- | --- |
| Capability | Exact and bounded-memory approximate joint point-prediction draws, bounded-memory exact point prediction, projected regular grids, point summaries, exceedance probabilities, grouped weighted targets and exact polygon-cell area weights |
| Public API | `sp_predict_draws()`, `sp_predict_draws_neighbor()`, `sp_predict_target_neighbor()`, `sp_neighbor_covariance()`, `sp_predict_chunked()`, `sp_grid()`, `sp_polygon_weights()`, `sp_polygon_area_weights()`, `sp_predict_target()`, `sp_target_geometry()` and print/plot methods |
| Joint uncertainty | Exact Gaussian conditional covariance or an explicit order-dependent nearest-predecessor factorization; Laplace Gaussian link covariance for binomial/Poisson fits |
| Response draws | Gaussian values, binomial probabilities or Bernoulli outcomes, Poisson intensities or counts |
| Grid policy | Known projected CRS only; polygon union, explicit boundary or point convex hull; stable top-to-bottom/left-to-right identifiers |
| Target policy | Aggregation within every joint draw; non-negative weights normalised within group; draw quantiles and exceedance frequencies |
| Parameter uncertainty | Optional local observed-Hessian draws with fitted-state reconstruction; Gaussian regression/log-covariance block independence is explicit |
| Polygon allocation | Either projected point-within-polygon centre allocation or exact projected cell-region intersection areas; optional density weights and explicit unmatched/empty-region policies |
| Derivation | `docs/derivations/06_JOINT_PREDICTION_GRIDS_AND_TARGETS.md`, equations J-01--J-20; `docs/derivations/18_CHUNKED_POINT_PREDICTION.md`, equations CH-01--CH-03; `docs/derivations/19_NEAREST_PREDECESSOR_PREDICTION.md`, equations NP-01--NP-14; and area-target equations D-11--D-14 |
| Synthetic checks | Direct Gaussian joint covariance and seeded transform, exact chunked/one-shot Gaussian and binomial equality, same-seed full-predecessor factorization, reordered moment agreement, grouped unseen-level preservation, PSD duplicate handling, fitted-parameter state propagation, streamed numeric/repeated-polygon targets, non-Gaussian marginal consistency, response support, density weighting, identifier joins, CRS and boundary failures |
| Not yet included | Exact cross-chunk covariance, non-Gaussian/grouped/temporal scalable training likelihoods and consumer-specific target schemas; the separate Gaussian nearest-predecessor likelihood is recorded above |
| Status | Implemented, checked and installed; package-level verification recorded below |

## Spatial audit and target-display slice

| Field | F09/F25/M02/M03 supporting record |
| --- | --- |
| Capability | Pair-distance audit, explicit UTM CRS selection, and exact joins of unrounded target summaries to point/grid or polygon geometry |
| Public API | `sp_distance_summary()`, `sp_utm_crs()`, `sp_target_geometry()` and its plot method |
| Distance policy | Euclidean projected distances; exact repeated-coordinate detection; include, exclude or reject duplicates; type-8 quantiles |
| CRS policy | Standard WGS84 UTM zone/hemisphere choice, including Norway/Svalbard exceptions; explicit errors or centre policies for multi-zone and equator-spanning inputs |
| Geometry join | Exact identifier match, or declared same-row order; the numerical summary remains unrounded in the returned `sf` object |
| Derivation | `docs/derivations/11_SPATIAL_AUDIT_ENVELOPES_AND_AREA_TARGETS.md`, equations D-01--D-04 and the target-display algorithm |
| Synthetic checks | Direct pair enumeration and duplicate policies; northern/southern, zone-boundary, longitude-180, Norway, Svalbard, polar, equator and antimeridian cases; reordered identifier joins and mismatch failures |
| Not yet included | Automatic cartographic styling or a general projection recommender outside the stated UTM domain |
| Status | Implemented locally; final package verification is recorded below |

## Cross-validation and predictive-assessment slice

| Field | F02/F03/M08/M09 supporting record |
| --- | --- |
| Capability | Fixed-fold refitting with held-out predictions, draw-based recovery assessment, predictive interval calibration, named multi-scenario recovery and validation/simulation displays |
| Public API | `sp_cross_validate()`, `sp_group_folds()`, `sp_simulation_assess()`, `sp_simulation_grid()`, `sp_simulation_scenarios()`, `sp_calibration()` and print/plot methods |
| Fold policy | Preserve supplied identifiers, keep declared groups whole, create balanced locally seeded random folds, or exclude a declared distance around geometry-driven test folds; full refit within every resulting training set |
| Metrics | RMSE, MAE and negative log score overall/by fold; non-Gaussian response mass is integrated over the approximate Gaussian link prediction; draw bias, coverage, interval width and empirical CRPS |
| Leakage boundary | Held-out observations are excluded before mean and covariance fitting |
| Derivation | `docs/derivations/07_CROSS_VALIDATION_AND_SIMULATION_ASSESSMENT.md`, equations V-01--V-16 and explicit display algorithms |
| Synthetic checks | Manual independent fold refit, local/reproducible balanced folds, Gaussian and binomial predictions, brute-force CRPS, direct exceedance frequencies, direct multi-level coverage/Wilson bounds and rendering of every display view |
| Not yet included | Further domain-specific fold constructors without a declared consumer use case |
| Status | Implemented, checked and installed; package-level verification recorded below |

## Non-Gaussian simulation, diagnostics and spatial-fold slice

| Field | F02/F05/F11/F29/M07 supporting record |
| --- | --- |
| Capability | Fitted binomial/Poisson response simulation, conditional residual diagnostics, coordinate-clustered validation folds and repeated parameter recovery |
| Public API | `simulate()`, `sp_diagnostics()`, `sp_spatial_folds()`, `sp_simulation_recovery()`, `sp_simulation_grid()` and print/plot methods |
| Simulation | Unconditional fitted latent Gaussian draw followed by fitted-trial binomial or Poisson response generation |
| Residuals | Response, Pearson, signed deviance and caller-seeded randomized quantile residuals evaluated at the final latent mode |
| Fold policy | K-means on unique coordinates; duplicate locations remain together; labels canonicalized top-to-bottom then left-to-right; optional buffers exclude near-test training rows per fold |
| Recovery | Repeat the declared fitting route, including MCML when applicable; report successful count, bias, RMSE, empirical spread, mean reported SE and Wald coverage |
| Derivation | `docs/derivations/09_NONGAUSSIAN_SIMULATION_DIAGNOSTICS_AND_SPATIAL_FOLDS.md`, equations S-01--S-24 |
| Synthetic checks | Exact seeded latent/response transform, family support, direct Pearson/deviance and randomized-CDF definitions, seed isolation, duplicate-location protection, fold refitting, independent first recovery refit and direct bias/RMSE summaries |
| Not yet included | Further domain-specific folds without a declared use case |
| Status | Implemented, checked and installed; package-level verification recorded below |

## Structured reporting and consumer-contract slice

| Field | F30/M04/M06 record |
| --- | --- |
| Capability | Stable unrounded model-estimate, prediction-summary and validation-metric tables, plus formatted views and an assembled analysis report |
| Public API | `sp_model_table()`, `sp_prediction_table()`, `sp_validation_table()`, `sp_format_table()` and `sp_report()` |
| Model schema | Component, term, estimate, standard error, Wald statistic/p-value where appropriate, estimation status and distance units |
| Prediction schema | Target identifier, draw mean/SD, requested type-8 quantiles, optional exceedance probability, family, target and scale |
| Validation schema | Assessment type, scope, identifier, metric and value for cross-validation, predictive assessment and simulation recovery |
| Consumer evidence | Read-only exact call/schema audit at immutable MBGapp, lsda/course and SDALGCP2 revisions in `docs/12_CONSUMER_AUDIT.md` |
| Presentation contract | Fixed-decimal display, explicit small-p-value threshold and missing marker; every formatted object retains its exact source table |
| Report contract | Model, optional prediction and optional validation sections are assembled only from the public stable tables |
| Derivation | `docs/derivations/10_REPORTING_AND_CONSUMER_CONTRACTS.md`, including R-01 and explicit assembly algorithm; statistical definitions continue to derive from chapters 06, 07 and 09 |
| Synthetic checks | Fixed schema, direct coefficient/prediction calculations, grouped target, overall/fold validation, exact formatted strings, raw-table identity, report sections and failure contracts |
| Not yet included | Journal-specific typesetting, HTML/office export and authorized consumer integration edits |
| Status | Computational and dependency-free presentation layers implemented locally |

## Prediction scenario update slice

| Field | F31 record |
| --- | --- |
| Capability | Matched baseline/scenario point prediction after selected predictor or in-formula offset changes |
| Public API | `sp_prediction_update()` and its print method |
| Row policy | Logical, unique integer or exact unique row-name selection; named replacement values must exactly match the selected identifiers |
| Design policy | Rebuild both model frames, design matrices, transformations, factors, contrasts and formula offsets from the fitted terms |
| Preservation policy | No refit or fitted-object mutation; exact row order, coordinates, `sf` geometry and every unmentioned data column are verified unchanged |
| Interpretation | Conditional descriptive prediction contrast under a common fitted state; not labelled as a causal effect |
| Derivation | `docs/derivations/12_PREDICTION_SCENARIO_UPDATES.md`, equations P-01--P-05 and explicit preservation algorithm |
| Synthetic checks | Independent manual data update and dual prediction calls; reordered named values, log-exposure offset, standard errors, factor levels, retained columns, coordinates, geometry and failure contracts |
| Not included | Coordinate relocation, separately supplied offset reconstruction, refitting or causal intervention claims |
| Status | Implemented locally; package-level verification is recorded below |

## Gaussian grouped random-intercept slice

| Field | F26 record |
| --- | --- |
| Capability | One or more additive Gaussian latent random-intercept grouping variables within Gaussian, binomial-logit or Poisson-log spatial models |
| Public API | `sp_fit(..., random = ~ group1 + group2)` with `sp_control(random_variance = ...)`; ordinary prediction, simulation, variogram, validation and reporting methods |
| Covariance | Independent zero-mean Gaussian effects by observed level and grouping variable; variance contributions add to the Matérn-plus-nugget covariance |
| Estimation | Fixed positive variances; joint log-scale exact Gaussian ML or non-Gaussian Laplace estimation with spatial covariance parameters; delta-method variance-component standard errors when curvature is usable |
| Prediction policy | Conditional prediction at observed levels; unseen levels fail by default or use an explicit zero-mean marginal policy that retains fitted variance and shared-new-level dependence; the policy is explicit in direct prediction, log scores and cross-validation; fitted random variances participate in local Hessian parameter draws |
| Derivation | `docs/derivations/13_GAUSSIAN_GROUPED_RANDOM_INTERCEPTS.md`, equations Q-01--Q-16 and Q-10a--Q-10c |
| Synthetic checks | Direct grouped covariance, GLS/Laplace likelihood, known-level point/joint prediction, seeded simulation, group-aware Gaussian variogram, refit and MCML workflows, parameter-state rebuilding, reporting and failure contracts |
| Not included | Random slopes, correlated effects or non-spatial grouped fitting in this independent-intercept slice; see the fixed-covariance grouped-slope slice below |
| Status | Implemented locally; final package verification is recorded below |

## Grouped random-slope and covariance-estimation slice

| Field | F26 extension record |
| --- | --- |
| Capability | One or more additive grouped numeric slopes, optionally paired with a correlated random intercept, in Gaussian, binomial-logit and Poisson-log spatial models |
| Public API | `sp_random_slope(group, slope, covariance, intercept = TRUE, estimate = FALSE)` and `sp_fit(..., random_slope = specification)`; lists add multiple components |
| Covariance | For rows (i,j), component (r) adds `I(group_i == group_j) w_i' G_r w_j`; every supplied matrix is finite, symmetric and positive definite |
| Estimation | Each covariance may remain fixed or opt into log-Cholesky constrained estimation. Estimated natural variances/covariances receive delta-method Hessian uncertainty and fitted degrees of freedom; fixed entries remain assumption rows without invented standard errors |
| Prediction policy | Conditional covariance for observed levels; default error for unseen levels; explicit marginal policy gives zero training cross-covariance while preserving covariance between prediction rows sharing a new level |
| Workflow coverage | Gaussian ML, Laplace fitting, MCML prior, point/joint prediction, response simulation, variograms, prediction parameter draws, validation refits, simulation recovery and stable reporting |
| Derivation | `docs/derivations/17_GROUPED_RANDOM_SLOPES.md`, equations RS-01--RS-11 |
| Synthetic checks | Direct covariance, GLS estimate, likelihood, known-level prediction, shared-unseen-level joint covariance, Laplace covariance, MCML state, positive-definite fitted matrices, exact Cholesky reconstruction, analytic-versus-numerical delta Jacobian, simulation recovery, report status and identifiability/data contracts |
| Not included | Random nonlinear bases, formula shorthand, non-spatial grouped fitting or boundary-corrected covariance inference |
| Status | Implemented locally; final package verification is recorded below |

## Decay-adjusted intervention-history and impact slice

| Field | F08/W01 foundation record |
| --- | --- |
| Capability | Auditable intervention histories, fixed or constrained fitted decay-adjusted prevalence parameters, and effect display in binomial spatial models |
| Public API | `sp_intervention()`, `sp_intervention_multiplier()`, and `sp_fit(..., intervention = specification)` with ordinary MCML, prediction, simulation, validation and reporting methods |
| History contract | Unique unit/time rows, optional zero-one delivered indicator, finite numeric times and exact unit matching |
| Timing policy | Only interventions strictly before a survey contribute; a same-time survey is interpreted as occurring before rollout |
| Model | Observed prevalence equals baseline logistic prevalence times the product of one minus each exponentially decayed intervention effect |
| Parameter status | Immediate reduction `alpha` and positive decay may be fixed or jointly estimated on logit/log scales by Laplace or MCML fitting; an optional normalized beta penalty is explicit |
| Display | `sp_intervention_curve()` returns unrounded reduction/multiplier curves and delta-method pointwise intervals, with print and plot views |
| Derivation | `docs/derivations/14_DECAY_ADJUSTED_INTERVENTION_HISTORY.md`, equations I-01--I-19 |
| Synthetic checks | Hand-calculated multipliers, finite-difference latent score/curvature, ordinary-binomial limit, integrated point means, exact seeded joint-response transform and simulation, penalty reconstruction, transformed covariance, curve uncertainty, display, MCML fixed-parameter path and fold refits |
| Not included | Consumer-specific intervention decision summaries |
| Status | Implemented locally; final package verification is recorded below |

## Separable space-time fitting slice

| Field | F19 supporting record |
| --- | --- |
| Capability | Separable Matérn-exponential residual covariance for point observations with numeric survey times |
| Public API | `spacetime_cor()`, `sp_temporal()` and `sp_fit(..., time = specification)` with ordinary prediction, simulation, variogram, validation and reporting methods |
| Covariance | Matérn spatial correlation multiplied by exponential temporal correlation; ordinary nugget and optional Gaussian grouped effects remain additive |
| Estimation | Fixed or joint log-scale temporal-scale estimation for every supported family: exact Gaussian ML, non-Gaussian Laplace approximation and subsequent MCML correction |
| Prediction | Training-to-new and new-to-new covariance both use declared prediction times; point and joint uncertainty retain temporal dependence |
| Diagnostics and validation | Fitted variograms average pair-specific time-aware semivariances within distance bins; refit workflows preserve the temporal specification |
| Derivation | `docs/derivations/15_SEPARABLE_SPACETIME_CORRELATION.md`, equations T-01--T-17 |
| Synthetic checks | Correlation limits and PSD; direct space-time GLS and likelihood; new-time point/joint covariance; exact seeded simulation; pairwise variogram; fixed- and estimated-scale binomial fits; direct non-Gaussian fitted covariance; MCML propagation; cross-validation; reporting and contract failures |
| Not included | Non-separable, anisotropic, seasonal or non-stationary space-time covariance; scalable sparse approximations or boundary-corrected scale inference |
| Status | Implemented locally; final package verification is recorded below |

## Transformed parameter-interval slice

| Field | Supporting record |
| --- | --- |
| Capability | Local support-respecting Wald intervals for fitted parameters |
| Public API | `confint()` method for `spmodels_fit` |
| Transform policy | Identity for mean coefficients, log for positive parameters, logit for intervention impact |
| Missing-curvature policy | Missing standard errors remain missing limits; fixed parameters are omitted |
| Derivation | `docs/derivations/16_PARAMETER_CONFIDENCE_INTERVALS.md`, equations H-01--H-06 |
| Synthetic checks | Direct identity, log and logit reconstruction; support checks; exact selection and failure contracts |
| Not included | Likelihood profiling, boundary corrections, simultaneous intervals or analytic prediction propagation |
| Status | Implemented locally; final package verification is recorded below |

## Reference-validation preparation

The neutral case matrix in `docs/validation/01_INITIAL_CASE_MATRIX.md`
predeclares deterministic tolerances, stochastic equivalence rules, required
outputs and independent checks. No comparison has been run and no case is
recorded as a pass. The external pinned environment remains to be frozen.

## Evidence recorded on 12 October 2026

- `Rscript` execution of `tests/test-matern.R`: pass under R 4.6.1.
- Source help conversion with `R CMD Rdconv`: pass.
- `R CMD build --no-build-vignettes --no-manual`: pass under R 4.6.1.
- Clean temporary-library installation and installed smoke test: pass under
  R 4.6.1.
- `R CMD check --no-manual --no-build-vignettes`: `Status: OK` under R 4.6.1
  on Ubuntu 24.04.4 LTS. Network index access was unavailable during the
  dependency check, but the check continued and reported no NOTE, WARNING or
  ERROR.
- Reference comparison: not run; remains planned in the separate pinned
  validation environment.

## Evidence recorded on 13 October 2026

- `Rscript tests/test-matern.R`: pass under R 4.6.1.
- `Rscript tests/test-gaussian-fit.R`: pass under R 4.6.1.
- The Gaussian test independently reconstructs the exponential special case of
  the covariance, generalized least-squares coefficients, Gaussian
  log-likelihood, conditional prediction mean and both prediction variances.
- Consumer-source audit: the MBGapp, lsda and SDALGCP2 working sources were
  located elsewhere on Linux. They were inspected read-only to establish their
  identities; no call-site or migration claim was made.
- Source-package build: pass under R 4.6.1.
- Clean-library `R CMD check --no-manual --no-build-vignettes`: `Status: OK`
  under R 4.6.1 on Ubuntu 24.04.4 LTS. Network index access was unavailable
  during dependency checking, but no package NOTE, WARNING or ERROR resulted.
- Installation to the R 4.6 user library and a clean-session fitted-model and
  prediction smoke test: pass; covariance optimiser convergence code zero.
- Read-only consumer call-site audit identified immediate course requirements
  for variograms and Gaussian/binomial prediction and MBGapp requirements for
  continuous, prevalence and count fits with predictive grid draws.
- Projected point `sf` fitting and prediction, exact Gaussian simulation,
  empirical variograms, fitted variograms and residual diagnostics: implemented
  with direct tests.
- Covariance-parameter standard errors from the observed profile Hessian and
  delta method, with missing values when curvature is unusable: implemented.
- Independent Gaussian limiting route: direct agreement with `lm()` coefficient
  algebra and the ML residual-variance definition.
- Checked teaching workflow: lsda 0.2.0 `galicia`, explicitly transformed from
  WGS84 longitude/latitude to ETRS89 / UTM zone 29N; independent and spatial
  fits, variogram, diagnostics and 20 seeded simulations all run successfully.
  The data remain in lsda and an upstream bibliographic-source gap is recorded.
- Clean-library package check after this expansion: `Status: OK` under R 4.6.1
  on Ubuntu 24.04.4 LTS.
- Matérn scale gradient/curvature and binomial/Poisson latent mode: implemented
  independently with analytic, finite-difference, stationarity, curvature and
  row-permutation tests.
- Conditional MALA latent sampling and basic single-chain diagnostics:
  implemented with local seeded reproducibility, acceptance accounting, trace,
  autocorrelation and initial-positive-pair effective sample size. The
  one-dimensional Bernoulli fixture agrees with numerical-integration moments
  within predeclared stochastic margins.
- Binomial-logit and Poisson-log geostatistical fitting by a labelled Laplace
  approximation: implemented for Bernoulli, binomial count/proportion and
  Poisson count responses. The objective is independently reconstructed from
  its likelihood, prior quadratic and two determinant terms; equivalent
  binomial encodings, row permutation and Poisson offsets are tested.
- Non-Gaussian point prediction: implemented on link and response scales using
  the local Gaussian latent posterior. Poisson moments are checked against
  their lognormal identities and binomial means against direct numerical
  integration.
- Clean-library package check after non-Gaussian fitting and prediction:
  `Status: OK` under R 4.6.1 on Ubuntu 24.04.4 LTS. The installed-package
  binomial fit, response prediction, conditional sampler and sampler-diagnostic
  smoke test passed in a fresh R session.
- Joint prediction, grid and target expansion: direct covariance, seeded draws,
  nonlinear response support, weighted aggregation, exceedance, CRS and grid
  ordering tests all pass. Clean-library `R CMD check` reports `Status: OK`
  under R 4.6.1. A fresh installed-package workflow fitted a binomial spatial
  model, built a 16-point projected grid, generated 100 joint probability draws
  and calculated point and regional exceedance targets successfully.
- Polygon allocation and parameter-uncertainty expansion: projected
  point-to-region assignment, supplied point weights and direct regional
  reconstruction pass; local parameter draws are seed-reproducible and stored
  with the predictive result for Gaussian and non-Gaussian fits. The clean
  package check after this expansion reports `Status: OK` under R 4.6.1.
- Fixed-fold cross-validation and predictive-draw assessment: all source and
  installed-package tests pass, including a separately reconstructed held-out
  refit, preserved fold identifiers, local balanced random folds, Gaussian and
  binomial scores, brute-force empirical CRPS and direct exceedance checks.
  The rebuilt source tarball passes clean `R CMD check` with `Status: OK` under
  R 4.6.1 on Ubuntu 24.04.4 LTS. A fresh installed-package workflow completed
  three Gaussian refits, 100 joint draws and draw assessment successfully.
- Monte Carlo maximum-likelihood correction: the normalized importance
  estimate agrees with independent one-dimensional numerical integration; all
  binomial, Poisson, seed-isolation, direct weight, prediction and fold-wise
  correction tests pass. The clean source-package check reports `Status: OK`
  under R 4.6.1. A fresh installed-package workflow completed MCML correction,
  prediction and two-fold MCML cross-validation; its 500-draw diagnostic ESS
  was 445.2 and its estimated log-likelihood Monte Carlo standard error was
  0.0157 for that smoke-test dataset.
- Fitted non-Gaussian simulation, conditional diagnostics, geometry-driven
  folds and repeated recovery: exact seeded transforms, response support,
  direct Pearson/deviance definitions, repeated-location protection, manual
  refits and direct recovery summaries pass. Clean `R CMD check` reports
  `Status: OK` under R 4.6.1. A fresh installed-package binomial workflow
  completed simulation, diagnostics, spatial cross-validation and two
  successful simulation/refit recovery replicates.
- Reference comparison: not run.

## Evidence recorded on 14 October 2026

- Extended conditional MALA sampling from one chain to an explicit retained
  draw-by-component-by-chain array with chain-specific starts, acceptance rates
  and reproducible local RNG handling. The pooled matrix remains available for
  backwards-compatible summaries.
- Added chain selection, rank-normalized split/folded R-hat, pooled
  autocorrelation-aware ESS and mean MCSE. A fixed-chain test independently
  reconstructs R-hat, while four seeded posterior chains test convergence,
  identity, reproducibility, ESS bounds and the MCSE formula.
- Completed a read-only, revision-pinned audit of the exact fit, summary,
  prediction and grid-target schemas used by MBGapp and the week 8 course;
  recorded suitable lsda data schemas and SDALGCP2 interface evidence. No
  external project was changed.
- Added stable model, prediction and validation tables. Direct tests reproduce
  coefficient rows, scale units, draw means, exceedance frequencies, grouped
  targets, and overall/fold validation schemas from the underlying objects.
- All source tests, including `tests/test-reporting.R`, pass under R 4.6.1.
- The rebuilt source package passes clean
  `R CMD check --no-manual --no-build-vignettes` with `Status: OK` under
  R 4.6.1 on Ubuntu 24.04.4 LTS. A fresh installed-package workflow completed
  model, prediction and cross-validation table generation successfully.
- The reader guide and technical derivations PDFs were regenerated with the
  reporting chapter included and their rendered pages inspected.
- Predeclared the initial neutral reference-validation matrix and distinct
  deterministic/stochastic acceptance rules. No reference comparison has been
  run.
- Added pair-distance summaries with explicit duplicate policies and a UTM CRS
  selector with tested zone, hemisphere, boundary and high-latitude rules.
- Added permutation and fitted-model residual-variogram envelopes. Every curve,
  standardized global statistic and finite-simulation p-value is retained and
  independently reconstructed in the seeded tests.
- Added exact polygon-cell intersection-area weights, optional density weights,
  explicit unmatched/empty-region policies, and identifier-safe joins of
  unrounded target summaries to `sf` geometry.
- Added the full teaching derivation for these operations in chapter 11 and
  regenerated all three project PDFs. The four new technical-derivation pages
  were inspected after rendering; equations D-01--D-14 and the code/test map
  are legible.
- All source tests, including the new geometry, target and variogram-envelope
  tests, pass. The final rebuilt source package passes
  `R CMD check --no-manual --no-build-vignettes` with `Status: OK` under
  R 4.6.1 on Ubuntu 24.04.4 LTS.
- The final tarball was installed to the R 4.6 user library. A fresh-session
  installed-package workflow completed distance auditing, UTM selection,
  fixed-parameter spatial fitting and a 19-curve variogram envelope.
- Added multi-level predictive interval calibration with independently checked
  type-8 coverage values and Wilson bounds. Added displays for calibration,
  cross-validation, predictive assessment and one-parameter recovery; every
  view renders successfully from fixed unrounded test values.

## Evidence recorded on 15 October 2026

- Added matched prediction-scenario updates with exact named-row alignment,
  formula/design reconstruction and explicit preservation of fitted state,
  row order, coordinates, geometry and unmentioned columns. Direct tests agree
  with independently updated data and two ordinary prediction calls, including
  in-formula log-exposure offsets and standard errors.
- Added presentation-only fixed-decimal tables and an assembled report object.
  Tests confirm every raw report section is identical to its public stable
  table, small p-values are not printed as zero, and missing markers never
  replace machine-readable values.
- Added the complete prediction-update derivation and reporting algorithm,
  regenerated all three PDFs, and visually inspected the newly rendered pages
  for equations, tables and line breaks.
- All source tests, including `tests/test-prediction-update.R` and
  `tests/test-report-presentation.R`, pass under R 4.6.1.
- The rebuilt source package passes
  `R CMD check --no-manual --no-build-vignettes` with `Status: OK` under
  R 4.6.1 on Ubuntu 24.04.4 LTS.
- Added one or more additive Gaussian random intercepts with fixed or fitted
  positive variance components. Direct tests reconstruct the complete
  covariance, GLS estimates, likelihood, known-level point/joint prediction,
  seeded simulation and pair-specific variogram contribution. Refit workflows,
  table rows and explicit unseen-level/unsupported-structure failures pass.
- Added robust MCML orchestration. Independent replicate summaries reproduce
  direct column means and standard deviations while preserving the caller's
  RNG state. Precision-directed fitting stops at passing MCSE/ESS targets or
  evaluates the declared cap exactly, records every attempt and never relabels
  a failed target as a pass.
- Added the fixed decay-adjusted intervention-history foundation for binomial
  spatial models. Hand-calculated timing multipliers, finite-difference latent
  score/curvature, integrated response prediction, exact seeded prediction and
  simulation transforms, MCML and fold-wise refits all pass. Immediate impact
  and decay remain fixed inputs; their estimation is not claimed.
- All source tests, including `tests/test-random-effects.R`,
  `tests/test-mcml.R` and `tests/test-intervention.R`, pass. The rebuilt source
  package again passes clean `R CMD check --no-manual --no-build-vignettes`
  with `Status: OK` under R 4.6.1 on Ubuntu 24.04.4 LTS.
- The implemented guide and technical derivations were regenerated. The new
  MCML, grouped-effects and five intervention pages were visually inspected;
  equations C-13--C-16, Q-01--Q-12 and I-01--I-14 are legible.
- An extra eight-row Bernoulli intervention smoke run with only 100 importance
  draws ended at an unusable MCML latent state. This is retained as a numerical
  limitation rather than converted into a pass. The predeclared count-response
  intervention fixture passes source, installed-package and MCML checks with
  the same draw count; larger/replicated runs remain necessary for sparse
  binary data.
- Extended the intervention path with constrained joint Laplace estimation of
  immediate reduction and decay, an optional normalized beta penalty,
  observed-Hessian delta-method covariance and an intervention-effect curve.
  Direct tests reconstruct the penalized criterion, parameter count, curve
  gradient variance, interval bounds and rendered plot. MCML re-estimation of
  these impact parameters remains explicitly rejected rather than silently
  fixing them.
- Added separable Matérn-exponential space-time covariance throughout the
  fitted-model lifecycle. Direct tests reconstruct Gaussian GLS and likelihood,
  prediction at new times, the full joint predictive covariance, seeded
  simulation and pair-specific fitted variograms. Fixed-scale binomial,
  MCML, cross-validation and reporting paths also pass; Gaussian temporal-scale
  estimation is positive and counted explicitly in model degrees of freedom.
- Extended the same positive log-scale temporal estimation to binomial-logit
  and Poisson-log Laplace fits and MCML correction. Direct tests reconstruct
  the fitted latent covariance, verify parameter counting and reporting, and
  retain the fitted value in the temporal specification used downstream.
- Every source test, including `tests/test-spacetime-fit.R`, passes. The rebuilt
  package passes `R CMD check --no-manual --no-build-vignettes` with
  `Status: OK` under R 4.6.1 on Ubuntu 24.04.4 LTS.
- Added named multi-scenario simulation recovery with unequal replicate counts,
  deterministic per-scenario seeds, complete retained results, combined
  unpooled summaries and validation-table output. Direct tests show each
  scenario is identical to its separately seeded recovery call and preserves
  the caller's random-number state.
- Added `confint()` for fitted models with identity, log and logit local Wald
  intervals chosen by parameter support. Direct tests cover mean coefficients,
  nugget variance, intervention impact/decay and temporal scale; missing
  curvature is retained rather than replaced.
- Added seeded randomized quantile residuals for fitted binomial and Poisson
  models. The test reconstructs each discrete CDF interval, uniform draw and
  normal-quantile transform exactly and verifies caller RNG isolation.
- Replaced plug-in non-Gaussian validation scores with one-dimensional
  integration of binomial or Poisson response mass over the approximate
  Gaussian predictive link distribution. Independent binomial and Poisson
  quadrature checks and the direct Gaussian density check pass.
- Added buffered spatial folds. A direct reconstruction verifies every row
  within the fitted-unit buffer is removed from one fold's training data and
  reproduces its held-out means and standard errors; an empty buffered
  training set fails explicitly.
- Extended MCML to re-estimate constrained intervention impact and decay,
  including the optional normalized alpha penalty and delta-method covariance.
  The test reconstructs every final likelihood, Gaussian-prior, proposal and
  penalty contribution to the importance weights and verifies row multipliers.
- Extended joint predictive parameter draws to fitted intervention impact and
  decay. Every draw rebuilds its training latent state and prediction-row
  multiplier; tests reconstruct the complete returned multiplier matrix from
  the stored transformed parameter draws.
- Added binomial/Poisson latent random intercepts under the existing grouping
  contract. Tests directly reconstruct the complete latent covariance,
  known-level point and joint prediction, exact seeded response simulation,
  refitting, reporting, positive variance estimation and unsupported MCML,
  unseen-level and parameter-draw boundaries.
- Extended normalized MCML and predictive parameter draws to grouped
  non-Gaussian fits. Fixed and estimated random variances enter the Gaussian
  latent prior consistently; seeded tests cover finite MCML weights, retained
  variance curvature and covariance rebuilding during prediction draws.
- Added an explicit marginal prediction policy for grouping levels absent at
  fit time. Direct Gaussian and binomial calculations verify zero training
  cross-covariance, retained new-level variance and covariance shared only by
  prediction rows carrying the same unseen level. The safe default remains an
  error.
- Audited lsda 0.2.0 teaching-data support without changing that project. The
  `loaloa` binomial input is pinned by SHA-256 and a projected 197-village fit
  converged under R 4.6.1 with log likelihood -665.5656. Its seeded randomized
  residual, three-fold spatial validation, integrated held-out score and local
  parameter-draw prediction paths all completed; the validation RMSE was
  0.13238 and mean negative log score was 4.15339 in this checked run.
- Recorded why the available areal counts, incidence rate, undocumented grid
  alignment and absent intervention history cannot truthfully be forced into
  the remaining point-supported teaching examples.
- Corrected packaging-only regression tests that had depended on private
  namespace helpers; their direct stable log-mean-exp, Cholesky-solve and
  Euclidean-distance calculations now also run against an installed package.
- The rebuilt source package passes clean
  `R CMD check --no-manual --no-build-vignettes` with `Status: OK` under
  R 4.6.1 on Ubuntu 24.04.4 LTS after these extensions.
- Completed the recovery parameter vector: estimated grouped variances and
  temporal scale now join mean, spatial/nugget and intervention parameters in
  retained truths, replicate estimates, standard errors and all recovery
  summaries. One-replicate refits verify the qualified names and exact truths.
- Regenerated and visually inspected all project PDFs. The 68-page technical
  derivation renders the complete recovery vector and marginal unseen-level
  equations clearly; the 27-page reader guide reflects the same boundaries.
- Installed the source package into the R 4.6 user library. A fresh-session
  installed-package smoke test fitted a model with estimated grouped and
  temporal variance components, performed marginal unseen-level point and
  joint prediction, and retained both components in simulation recovery.
- Added regular spatial-block validation folds. Direct tests reconstruct every
  block index and canonical fold label, verify whole-block retention and local
  RNG isolation, exercise random allocation, and complete a block-fold
  cross-validation run.
- Added environmental-space folds with formula/name selection, default
  standardization, repeated-coordinate consistency and canonical centres.
  Tests verify equivalent interfaces, local RNG isolation and a complete
  environmental-fold cross-validation run.
- Added group-constrained folds with deterministic size-first or seeded random
  ordering followed by least-loaded allocation. Direct tests verify the greedy
  allocation, whole-group retention, RNG isolation, audit metadata and a
  complete constrained cross-validation run.
- Added optional forked parallel dispatch across named recovery scenarios on
  non-Windows systems. Scenario seeds are resolved before dispatch; direct
  tests establish exact sequential/parallel estimate equality and preserve the
  caller's RNG state.
- Added optional forked parallel dispatch across replicates within one recovery
  scenario. Responses and MCML seeds are fixed before dispatch; direct tests
  establish exact sequential/parallel simulations, estimates, standard errors
  and convergence codes while preserving the caller RNG state.
- Added precision-directed recovery orchestration. It monitors the direct
  empirical-SD-over-root-successful-count MCSE for mean bias, requires every
  parameter and a declared minimum success count, escalates geometrically and
  evaluates the exact cap without weakening criteria. Tests cover both the
  first-attempt pass and failed capped paths.
- Precision-directed attempts now pass the declared replicate-worker count
  through to complete recovery runs. A direct test verifies identical estimates
  and stopping history between one-worker and two-worker execution.
- Recorded the first installed-package dense-workflow timing baseline on the
  stated Linux hardware. Estimated Gaussian fits took 1.454 seconds at n=100
  and 40.792 seconds at n=225; the latter is now the explicit first profiling
  target. The record avoids extrapolating a complexity rate from one timing
  pair.
- Profiled the n=225 Gaussian fit: the generic Bessel evaluation dominated the
  default smoothness-0.5 covariance path. Added exact half-integer Matérn fast
  paths for smoothness 0.5, 1.5 and 2.5. A controlled 225-by-225 kernel A/B
  check measured a 6.36-fold median speedup with maximum absolute difference
  1.11e-16; the existing analytic Matérn and fitted-model tests still pass.
- Ran a separate local CRAN-style audit. All computational, help, example and
  test stages passed with no error or warning; two release notes remain for the
  development version/licence and non-standard top-level planning documents.
  The clean ordinary-check record remains preserved separately as `Status: OK`.
- Added empirical residual variograms for binomial and Poisson fits using
  Pearson, signed deviance or seeded randomized-quantile residuals. Direct
  tests reconstruct every binned pair semivariance, preserve RNG state and
  render the empirical-only plot; no fitted discrete response-scale curve is
  invented.
- Extended residual-variogram permutation envelopes and the standardized global
  curve test to binomial and Poisson residual choices. Seeded tests reconstruct
  every permuted curve and render the empirical-only envelope.
- Added fitted-model parametric envelopes for binomial and Poisson diagnostics.
  Each simulated response is refitted by the declared Laplace or MCML route
  before residual calculation; convergence codes and captured warnings are
  retained, and a seeded test independently replays the first complete curve.
- Added Cartesian parameter-grid construction for recovery experiments. The
  constructor accepts only parameters recoverable under the unchanged refit
  specification, retains exact grid order, rebuilds Gaussian or non-Gaussian
  generative covariance state, and composes with named sequential/parallel
  scenario recovery. Tests independently reconstruct the Gaussian covariance
  and the non-Gaussian Cholesky identity.
- Rebuilt and installed the final source tarball after the grid and full-refit
  envelope additions. The clean source check is `Status: OK`; a fresh-session
  installed-package smoke test exercised both APIs under R 4.6.1.
- Profiled an estimated Gaussian fit at n=400 after the half-integer fast-path
  work. Dense Cholesky factorization was the largest sampled self-time. The
  internal covariance path now avoids redundant whole-matrix validation after
  inputs have already been checked. A five-repetition, same-process 400-by-400
  kernel comparison measured medians of 1.083 and 0.477 seconds (2.27-fold)
  with maximum absolute difference zero. Specialized and generic paths have
  direct equality and attribute-preservation tests; no whole-fit speedup or
  large-n scalability claim is inferred from separate-process fit timings.
- Repeated exact estimated Gaussian fits three times per size in one Linux R
  process after optimization. Median elapsed times were 0.233 seconds at
  n=100, 1.419 seconds at n=225 and 11.888 seconds at n=400; every fit
  converged and repeated its size-specific log likelihood. This records the
  measured dense boundary without presenting one-machine timings as a
  cross-platform guarantee or selecting an approximation prematurely.
- Added opt-in grouped random-slope covariance estimation through a
  log-diagonal Cholesky parameterization. Gaussian ML, non-Gaussian Laplace,
  MCML, prediction-state rebuilding, simulation recovery and parameter-grid
  scenarios retain the fitted positive-definite matrix. Direct tests cover
  exact parameter reconstruction, the analytic natural-scale Jacobian against
  finite differences, Hessian uncertainty, fitted degrees of freedom,
  combined slope/temporal indexing and explicit minimum-design failures.
- Added joint positive log-scale temporal estimation to binomial-logit and
  Poisson-log Laplace fits and MCML correction. Tests reconstruct the fitted
  space-time prior, verify reporting/counting, independent MCML natural-scale
  summaries and the combined temporal/random-slope path. All-equal survey
  times are rejected as unidentifiable.
- Replaced near-boundary invalid-trial objective values with a finite
  `1e20` sentinel in Gaussian, Laplace and MCML optimizers. A difficult
  Poisson space-time fixture now returns an auditable fitted state and
  convergence warning instead of failing inside numerical differentiation.
- Regenerated the 28-page implemented guide, 73-page technical derivations
  and 3-page derivation standard. The new grouped-covariance and temporal
  pages were inspected at full-resolution and render clearly. The 53-page
  package help manual also builds successfully.
- Rebuilt the source archive and completed a clean check under R 4.6.1 on
  Ubuntu 24.04.4: `Status: OK`, including all examples and 22 test files. A
  separate CRAN-style audit has no error or warning and retains exactly two
  declared release notes for the development version and top-level planning
  material. The checked archive was installed into the R 4.6 user library;
  a fresh session exercised estimated grouped covariance, chunked prediction
  and estimated binomial temporal scale successfully.
- Audited normalized MCML proposals on the authorised 197-row Loa loa teaching
  data. Declared wider scale mixtures sometimes produced much lower ESS and a
  dominant weight, including a same-seed 900-draw case with ESS fraction
  0.0010 versus 0.0202 for scale one. Added opt-in discarded-pilot adaptation
  with positive shrinkage toward the Laplace moments, separate pilot/main
  diagnostics, an ESS-gated fallback to the original moments and complete
  retained audit state. Synthetic tests independently
  reconstruct the adapted mixture density and verify conditional sampling,
  positive definiteness, deterministic seeds and the full-shrinkage boundary.
- Added bounded-memory nearest-predecessor predictive draws for spatial fits.
  The method constructs only local fitted predictive submatrices, records its
  order, neighbor graph and conditional variances, and never represents those
  local conditionals as an exact full covariance. Same-seed tests reproduce
  ordinary Gaussian and non-Gaussian link draws when all predecessors are
  retained; additional checks cover reordering, small neighborhoods, caller
  RNG isolation and shared unseen grouped levels. Dense training fitting and
  exact large-grid covariance remain explicitly outside this extension.
- Extended the nearest-predecessor route with deterministic maximin ordering.
  It starts at row one, updates one running minimum-distance vector and uses
  input order for ties, avoiding a full grid distance matrix. Tests reconstruct
  its exact order and verify that an explicit user permutation takes priority.
- Recorded a 25-location synthetic approximation diagnostic with four
  predecessors and 50,000 draws. Relative covariance error was 0.04994 in
  input order and 0.03927 under maximin order; equally weighted regional-mean
  variance ratios were 1.08308 and 1.05152. This is retained as one-design
  sensitivity evidence, not a universal ordering or accuracy claim.
- Added optional fitted-Hessian parameter uncertainty to bounded-memory joint
  draws. Every parameter draw rebuilds the fitted state and local factorization;
  intervention multipliers are draw-specific and retained. Same-seed tests at
  the all-predecessor boundary reproduce ordinary Gaussian parameter draws and
  fitted-intervention response draws, while the approximate route still
  returns no invented full covariance.
- Corrected prediction-side grouped random-slope covariance under fitted-
  parameter draws. Each rebuilt Cholesky state now supplies the new-site block
  as well as the training cross-covariance; a deterministic perturbation test
  verifies exact propagation instead of relying on a potentially unusable
  unconstrained Hessian draw.
- Added a checked positive-semidefinite predecessor solve for coincident latent
  prediction rows. It uses retained eigen-directions, verifies column-space
  compatibility and adds no jitter; duplicate-row tests produce identical
  draws and a zero conditional variance at the degenerate step.
- Added direct nearest-predecessor target streaming. Numeric, polygon-centre
  and repeated fractional-area allocations are normalized sparsely; each
  point-draw vector is released after its final graph use. Same-seed tests
  reproduce post-hoc Gaussian and binomial-latent targets and the returned
  object records peak active locations rather than promising a fixed width.
- In the 25-location/50,000-draw diagnostic, one streamed target occupied
  0.382 MiB versus 9.539 MiB for all point draws. Peak active point vectors
  were 10 in input order and 16 under maximin order, illustrating that ordering
  changes both approximation quality and graph storage.
- Extended MCML replicate and precision histories with adaptive-pilot ESS
  fraction, maximum weight, relative MCSE and used/fallback status. Tests cover
  non-adaptive missing fields and an intentionally strict pilot gate across
  independent replicates and a precision attempt.
- Added an analytic nearest-predecessor covariance diagnostic for manageable
  validation subsets. It recursively composes local conditional coefficients
  and variances, optionally compares with the exact fitted covariance, and
  reports relative Frobenius and maximum absolute errors. The all-predecessor
  Gaussian/non-Gaussian and duplicate-location boundaries pass directly. On
  the recorded four-neighbor design, analytic relative errors were 0.04425 in
  input order and 0.03307 under maximin order.
- Added the explicit Gaussian nearest-predecessor training likelihood. It
  profiles the linear mean from bounded local Gaussian conditionals, retains
  order/graph/conditional audit state, supplies ordinary extraction, interval,
  table and local point-prediction methods, and rejects unsupported grouped,
  temporal and non-Gaussian extensions. Input- and maximin-order tests recover
  dense likelihood, coefficients and coefficient covariance at the all-
  predecessor boundary. A complete MSc derivation records equations
  NL-01--NL-15 and the method's primary research sources.
- Recorded moderate and process-level diagnostics without selecting a neighbor
  budget. At 2,000 rows, the fixed neighbor-20 route used 130,496 KiB peak RSS
  and retained 1.59 MB, versus 300,132 KiB and 65.22 MB dense; the neighbor
  route completed 5,000 rows at 134,112 KiB. The pure-R local route was slower
  than dense BLAS at n=225, so documentation makes a memory result only and
  predeclares recovery, stability, peak-memory and cross-platform criteria.
- Executed the statistical-accuracy/stability portion of the predeclared
  Gaussian NNL study on 225 simulated datasets and
  retained all 2,250 neighbor-fit results in
  `docs/benchmarks/data/neighbor_likelihood_study.csv`, with an independently
  reproducible summary beside it. The archive used for the run had SHA-256
  `ead290f443dd1deb433ac7684e4bdcfed5dddb362135c44ccafb4022fc14f294`.
  All datasets completed. Per-fit elapsed time and external peak RSS were not
  captured in that initial run. Scenario-
  specific small budgets passed for the
  short-range cases, but medium/long-range cases at n=225 and n=400 required
  the dense all-predecessor boundary. The exact boundary's worst observed
  likelihood error was 5.12e-13; no generally adequate budget is claimed.
- Repeated all 2,250 NNL fits under the same seeds with per-fit elapsed-time
  recording. Every deterministic outcome matched the archived accuracy run.
  The complete process took 23:08.31 and peaked at 233,124 KiB RSS under GNU
  time. Per-configuration external RSS and cross-platform profiling remain
  outstanding; machine-readable timing results and their summary are archived
  alongside the original study.
- Executed the registered GRS study under R 4.6.1 on 1,400 usable fits across
  four cases. GRS-01/02 reached the 400-fit cap without meeting the scale
  bias-MCSE threshold; GRS-03/04 met Monte Carlo precision targets but their
  mean spatial-scale estimates remained about 1 against a generating value of
  3; the scale distributions were strongly right-skewed. All 1,400 Cholesky
  reconstructions and analytic-versus-finite-difference Jacobian checks passed.
  Hessian covariance was available for 937 fits, and local-Wald coverage was
  poor. Full per-fit estimates, uncertainty flags and checksummed summary are
  in `docs/validation/03_GROUPED_SLOPE_RECOVERY_RESULTS.md` and its data files.
  This is not a recovery, nominal-coverage or reference-agreement claim.
- Diagnosed the existing GRS Hessian subset without changing acceptance
  criteria. Estimated spatial-scale medians remained 0.313, 0.336, 0.893 and
  1.201 against truth 3; median local scale standard errors were 0.752--1.197,
  and local scale-interval coverage was 35.1%--59.6%. This is a descriptive
  range-identifiability warning, not a formal profile-likelihood result.
- Executed the registered TSC non-Gaussian temporal-scale study under R 4.6.1
  using archive SHA-256
  `22d059ff921242770fc1110af93b7ce9670020a9ddee19b97193f0841431b374`.
  All four cases reached 400 attempts; only 76.8%--83.0% of fits were usable,
  and natural temporal-scale bias-MCSEs were 267.982--452.270 against targets
  0.05--0.10. Usable estimates were severely right-skewed and some hit both
  optimizer bounds. Among 25 MCML response sets per case, only 1--6 met the
  all-three-repeat gate. All cases are inconclusive; the readable report,
  machine-readable raw/summary files and seed/gate validator are in
  `docs/validation/04_TEMPORAL_SCALE_RECOVERY_RESULTS.md` and adjacent files.
  The initial run's process exited nonzero after its final data write due to
  editing the runner while it was still being read; a checksum-verified resume
  preserved the original rows and completed all remaining stages. Combined
  wall time was 1:42:09.43; peak RSS across the two processes was 197,832 KiB.
- Executed the registered synthetic MCP proposal matrix on two response
  fixtures, four proposal policies, three main sample sizes and five
  independent attempts (120 fits), using source archive SHA-256
  `dd0c7f5d1532f70c9eeadf53f62558945dfa77eec6b4674c746f10f26608fce4`.
  Both base fits and all 120 MCML optimizations converged, but no policy met
  the 4-of-5 ESS/relative-MCSE gate on either fixture at 2,000 draws. The
  authorized 197-site Loa loa teaching-data diagnostic used the documented
  lsda 0.2.0 data hash and projection; all 20 attempts converged but showed
  severe weight degeneracy, so they provide no policy recommendation. Raw
  fits, objectives at a shared parameter, fixtures, summaries, operational
  gates, validation script and case-specific interpretation are archived in
  `docs/validation/05_MCML_PROPOSAL_RESULTS.md` and adjacent files.
  Synthetic and Loa runs took 4:50.59 and 27:34.37, peaking at 150,964 and
  249,808 KiB RSS, respectively.
- Executed the registered NNP study using source archive SHA-256
  `4bd9544a4ab011b08f5d7786c1d0db526f903ecc8cc93f08b23844f71247a7e7`.
  Five model cases, three grid sizes, two orders, four neighbor budgets and
  five seeds produced 3,000 target rows with no fit/prediction errors or
  warnings. All input-order exact-boundary draws matched their ordinary route;
  the maximum exact covariance discrepancy across both orders was below
  `1.55e-15`. For Gaussian, binomial, Poisson and fixed grouped-slope cases,
  candidate budgets 5/10/20 passed at 25/100/225 rows. The space-time case
  required 10 and 40 at 25 and 100 rows, and no registered budget passed at
  225 rows. Thus there is no universal budget recommendation. Raw per-case
  results, combined summaries, resource observations and the validator are
  linked from `docs/validation/06_NEIGHBOR_PREDICTION_RESULTS.md`.
- Final verification for this development block used all 24 source test files
  and a rebuilt source archive under R 4.6.1 on Ubuntu 24.04.4. `R CMD check
  --no-manual --no-build-vignettes` completed with `Status: OK`, including
  installed-package tests, examples, namespace/Rd consistency and loading
  checks. The checked archive SHA-256 was
  `8edff44855c7581a87b75edc86080f38ad4ed27d6bf4b74e3e1ac47cf2cc1cdc`.
  The refreshed implemented guide, technical derivations and derivation
  standard contain 30, 86 and 3 pages; the new MCML, nearest-predecessor
  prediction and nearest-predecessor likelihood equation pages were visually
  inspected at full resolution. This is a Linux development check, not evidence
  of cross-platform support or numerical agreement with the separately managed
  reference implementation. Per-configuration NNL external RSS and reference
  comparisons remain open.

## Deliberately unresolved release metadata

The package name remains provisional. The owner supplied the maintainer email
on 12 October 2026. A distribution licence has not yet been selected;
`DESCRIPTION` and `LICENSE` state this honestly, so the development snapshot is
not release-ready.

