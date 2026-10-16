# Initial reference-validation case matrix

## Rules

Deterministic comparisons use absolute and relative discrepancies

\[
d(a,b)=|a-b|, \qquad r(a,b)=\frac{|a-b|}{\max(1,|a|,|b|)}.
\]

Unless a row states otherwise, a deterministic scalar or array passes when
`d <= 1e-8 + 1e-6 * max(abs(reference))`. This is a pre-run engineering
criterion for ordinary well-conditioned examples, not a claim that agreement
has already been observed. Analytic fixtures in the package retain their
tighter test-specific tolerances.

Stochastic cases use independently seeded replicates. Agreement means the
difference in replicate means is no larger than three combined Monte Carlo
standard errors, and distributional coverage/quantile comparisons must also
meet the row-specific criterion. A single matching seed is not evidence of
stochastic equivalence.

## Cases

| ID | Configuration | Outputs | Criterion | Independent check |
| --- | --- | --- | --- | --- |
| RV-G-01 | Gaussian, fixed covariance, no nugget | coefficients, log likelihood, latent prediction mean/covariance | deterministic default | direct GLS and Gaussian conditioning |
| RV-G-02 | Gaussian, estimated covariance and nugget | estimates, SEs, log likelihood | deterministic default when both optimisers reach the same interior optimum; otherwise discrepancy review | profile likelihood and finite-difference curvature |
| RV-G-03 | Gaussian covariate grid | joint means, SDs, quantiles, exceedance | moments use deterministic default; empirical summaries within 3 combined MCSE | direct multivariate-normal transform |
| RV-B-01 | Binomial successes/failures, fixed covariance | Laplace estimates, latent mode, log likelihood | deterministic default after matching likelihood constants | direct score, Hessian and Laplace expression |
| RV-B-02 | Binomial MCML correction | estimates and likelihood | replicate mean within 3 combined MCSE; nominal 95% intervals cover the opposing replicate mean | numerical integration on a reduced fixture |
| RV-P-01 | Poisson with log exposure offset | estimates, fitted intensities, prediction moments | deterministic default for Laplace route | direct offset design and lognormal moments |
| RV-P-02 | Poisson MCML correction | estimates and likelihood | same stochastic rule as RV-B-02 | numerical integration on a reduced fixture |
| RV-T-01 | Weighted grouped prediction target | mean, SD, quantiles and exceedance | deterministic default when based on common stored draws | direct matrix multiplication within each draw |
| RV-V-01 | Fixed-fold Gaussian validation | row predictions and aggregate scores | deterministic default | independent fold refits and score formulas |
| RV-S-01 | fitted simulation and recovery | bias, RMSE, empirical SD and coverage | replicate mean within 3 combined MCSE; absolute coverage difference <= 0.03 with enough replicates for that precision | known generating parameters and binomial uncertainty calculation |

## Required fixtures

Each family uses a small synthetic analytic fixture and one suitable lsda
teaching dataset. Dataset version, support, row order, CRS transformation and
provenance must be recorded with the run. No data are copied into SPmodels by
this specification. Grid cases retain stable location identifiers and check
that predictor rows and geometries remain aligned.

## Not yet executed

Every case above is `planned`. Moving a case to `pass` requires a separately
stored raw result and discrepancy report. The project must not infer numerical
agreement from the package's internal tests alone.
