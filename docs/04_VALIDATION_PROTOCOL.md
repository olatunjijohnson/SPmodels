# SPmodels validation protocol

## Scope of this protocol

The reference comparisons below apply to equivalent models in the initial collection. Additional models need model-specific validation using independent derivations, analytic or limiting cases, simulation-based calibration or recovery as appropriate, and suitable trusted implementations where available. A comparison is meaningful only when the statistical models and targets agree. Every released model must also have a small reproducible teaching example whose outputs and interpretation are checked, using suitable lsda data where appropriate. Full MSc-accessible derivations are required for every implemented statistical and numerical method. Check dimensions, signs, derivatives, parameter conventions, limiting cases and approximation steps independently; link equations to implementation routines and validation cases. Review pedagogical completeness separately from numerical correctness.

## Meaning of the same result

The goal is to reproduce the same statistical model, estimand and observable behaviour as the selected reference implementation. Exact values are appropriate for discrete structural results; floating-point calculations require explicit numerical tolerances; random simulation requires evidence of equivalent distributions and sufficiently precise estimates. A shared random seed does not guarantee identical draws from independently implemented algorithms.

**Status: active protocol. Independent analytic/synthetic checks and local
dense-workflow benchmarks have been run and recorded; no numerical reference
parity comparison has yet been run.**

## Freeze the reference before implementation

In a separate validation environment, record the selected reference version and immutable build identifier, checksums, R and dependency versions, operating system, architecture, BLAS/LAPACK, spatial libraries, locale, thread count and RNG settings. Preserve the exact supported configuration and input checksums there. Keep comparator identity, implementation source and branding out of SPmodels project materials.

Confirm the first reference build at implementation kickoff using the owner-selected baseline and consumer requirements. A website displaying the same version number can document additional behaviour. Do not combine those surfaces silently. Preserve the original baseline when later adding comparisons against newer releases.

## Comparison record

Each case needs: case ID; feature ID; reference identifier; downstream consumer if any; synthetic fixture recipe or authorised data identifier; exact arguments; parameter mapping; seeds and RNG; expected schema; target quantity; absolute/relative tolerance or equivalence margin; execution environment; outputs; discrepancy; convergence/MCSE; elapsed time and memory; status and explanation.

Use statuses `planned`, `blocked`, `pass`, `fail` and `inconclusive`. Missing or inconclusive evidence cannot count as a pass. Store full precision results and a readable report. Keep changed reference outputs separate from the original snapshots.

## Test layers

| Layer | Examples | Required comparison |
| --- | --- | --- |
| Contracts | Arguments, row mapping, dimensions, names, factor levels, missingness, warnings | Exact discrete values or explicitly documented semantic mapping |
| Deterministic mathematics | Matérn covariance, derivatives, Gaussian conditional moments, fixed-draw summaries | Elementwise absolute plus relative tolerance and independent analytic checks |
| Optimisation | Gaussian fits, fixed-Monte-Carlo objectives, boundary solutions | Objective values, convergence, parameter scale and identifiable estimates |
| Random simulation | Latent draws, outcome simulation, MCMC, non-Gaussian estimates | Predefined equivalence margins with Monte Carlo uncertainty accounted for |
| Spatial behaviour | CRS, units, unique locations, grids, polygons, interpolation/aggregation | Aligned coordinates and row identities, numerical values, geometry/topology |
| Presentation | Tables, plots, summaries | Underlying quantities, labels, units and intervals; reviewed visual output |
| End to end | Dependent package fit to prediction to report | Same scientific conclusions and required outputs within declared criteria |

Initial deterministic tolerance candidates are `abs(a - b) <= atol + rtol * abs(b)`, with `atol = 1e-10`, `rtol = 1e-8` for well-conditioned kernels. These are planning values, not a universal promise. Define separate, justified tolerances for derivatives, ill-conditioned linear algebra, optimisation and CRS transformations before evaluating SPmodels. Check NA, NaN and infinite values separately. Matching two invalid values is not evidence of correct mathematics.

Optimisation checks must align likelihood constants, parameter transformations, starting values, constraints and convergence criteria. Compare the objective at matched parameters as well as at each fitted optimum. Weakly identified or boundary parameters need profile/objective comparisons and scientific outputs; arbitrary closeness of point estimates is insufficient.

## Stochastic equivalence

1. Define an application-relevant equivalence margin for each target before running comparisons. For prevalence and threshold probabilities use absolute probability margins; for regression and spatial parameters use meaningful units or a justified scale.
2. Run multiple independent chains and seeds for both systems, assessing stationarity, effective sample size and Monte Carlo standard error. Set a computation budget and a precision goal substantially smaller than the margin.
3. For independent estimates of a mean, use combined Monte Carlo uncertainty, such as `sqrt(MCSE_reference^2 + MCSE_new^2)`. Use quantile-specific and autocorrelation-aware uncertainty estimators for quantiles and tail probabilities.
4. Require the chosen confidence interval for the difference to lie wholly inside the equivalence interval. An interval that crosses the margin is inconclusive or a failure according to the predeclared rule; a large standard error must not make a discrepancy pass.
5. Compare marginal distributions and spatial cross-covariances, joint target distributions, tail probabilities and derived classifications. A non-significant difference test, overlapping ordinary intervals or similar means alone is insufficient.
6. Predeclare how multiple targets and repeated runs are handled. Increase precision according to a fixed schedule; do not keep rerunning until one seed passes.

When both interfaces can accept identical externally generated draws, compare deterministic downstream calculations using those draws as a separate test. This isolates arithmetic from sampler differences without requiring identical RNG streams in production.

## Required configurations

- Gaussian, binomial/Bernoulli and Poisson models; intercept-only and covariate-rich designs; continuous variables, factors, interactions and offsets.
- Binomial denominators and proportion encodings; Poisson exposure; supported custom links; supported observation weights and missing-data policies.
- Different Matérn smoothness/range/variance settings; fixed and estimated components; zero and near-zero components; duplicates and near-duplicates; singular or nearly singular cases.
- Supported grouping effects, repeated locations and GP time configurations; DAST histories, survey time versus GP time, coverage, fixed/estimated impact and penalties where available.
- Small and moderately sized datasets; reordered rows; reordered factor levels; prediction at observed and new locations; absent covariates and new group levels.
- CRS transformations and distance units; empty geometries; polygon boundaries, holes and multiple parts; grids that hit boundaries; regions containing no grid cells; zero or missing aggregation weights.
- Marginal and joint predictions; link, response, latent and observation targets; all component inclusion flags; raw draws, means, SDs, quantiles and exceedance probabilities.
- Spatial cross-validation using identical explicit folds; refitting and fixed-parameter modes; calibration and score definitions; simulation assessment; formatted summaries and plots.

Use small independent synthetic fixtures first. Include analytic limiting cases, closed-form Gaussian conditionals, covariance symmetry/positive semidefiniteness, finite-difference derivative checks, row-permutation invariance and predictable unit conversions. These guard against reproducing a reference defect.

## Downstream and release gates

The inventory rows link to feature groups here, then to concrete case IDs at implementation. No row is complete until all relevant argument branches and downstream scenarios are covered. A feature count alone does not establish behavioural coverage.

Run ordinary package tests and checks without the comparator installed. Run reference comparisons in the separate pinned validation environment. Scan direct and transitive dependencies, namespace imports, dynamic loads, namespace-qualified calls, internal calls and test/example dependencies.

When results disagree, check parameterisation, units, likelihood constants, prediction target and numerical uncertainty first. Record the discrepancy and an independent mathematical check. Do not copy a suspected reference bug, silently change the scope, or widen tolerance to conceal it; present the conflict for an explicit compatibility decision.

Complete coverage of the agreed initial scope is claimable only after every required capability and the owner's confirmed downstream workflows satisfy the gates. Cross-platform checks must report any remaining platform-specific limitations.
