# Grouped random-slope covariance recovery results

## Scope and decision rule

This records the predeclared GRS study in
`02_RECOVERY_AND_APPROXIMATION_DESIGN.md`. It is an internal synthetic
recovery study, not an external reference comparison. The complete fitted
collection and scenario summaries are in `data/grouped_slope_recovery.csv` and
`data/grouped_slope_recovery_summary.csv`.

The study used the registered four covariance scenarios and data-generation
recipe, with one observation per group at each of `z = -1, 0, 1`. It began
with 100 usable fits per scenario and escalated to 200 and then 400 when the
registered Monte Carlo standard error of estimated bias (bias-MCSE) targets
were not met. A fit counted as usable for point recovery only when estimates
were finite, the fitted grouped covariance was positive definite, optimizer
convergence was zero and no warning occurred. Bias-MCSE is a precision measure
for the estimated bias; passing it does not mean bias is negligible or that
the estimator is accurate. The observed biases and local-Wald coverage are
therefore reported alongside it.

## Results

| Scenario | Groups | Attempts / usable | Median fit time (s) | Hessian covariance finite | Scale mean / median (truth 3) | Scale bias-MCSE / target | Result |
|---|---:|---:|---:|---:|---:|---:|---|
| GRS-01 | 8 | 400 / 400 | 0.312 | 231 / 400 | 9.413 / 0.157 | 6.257 / 0.150 | Inconclusive: target not met |
| GRS-02 | 8 | 400 / 400 | 0.319 | 257 / 400 | 4.886 / 0.233 | 2.289 / 0.150 | Inconclusive: target not met |
| GRS-03 | 20 | 400 / 400 | 0.512 | 303 / 400 | 1.005 / 0.442 | 0.075 / 0.150 | Bias-MCSE target met; accuracy not established |
| GRS-04 | 20 | 200 / 200 | 0.516 | 146 / 200 | 0.991 / 0.417 | 0.085 / 0.150 | Bias-MCSE target met; accuracy not established |

All attempts had finite estimates, positive-definite fitted grouped
covariances, zero optimizer convergence codes and no warnings. GRS-01 and
GRS-02 reached the 400-fit cap but missed the registered scale bias-MCSE
target; the estimated spatial scale was highly dispersed and strongly biased
upward. GRS-03 and GRS-04 met all registered bias-MCSE thresholds at their
final sample sizes, but their mean scale estimates were about 1 rather than
the generating value 3. Their precision-target status is therefore not a
recovery or numerical-agreement claim. This is a clear warning that the
registered design does not recover the spatial range reliably in these
settings, even with 20 groups. All four scale distributions were strongly
right-skewed: their medians were between 0.157 and 0.442, and the largest
observed estimates were 2,262 for GRS-01 and 737 for GRS-02. The mean bias and
its MCSE therefore do not tell the whole story; quantiles and full replicate
values remain in the raw table.

The two registered mathematical checks passed on every fitted covariance:
Cholesky reconstruction of (G=LL^{\mathsf T}) and a central-finite-difference
check of the analytic natural-scale Jacobian. Both were checked for all 1,400
rows at tolerances (10^{-12}) and (10^{-6}), respectively.

Hessian-based covariance uncertainty was finite for 937 of 1,400 fits. Local
Wald coverage is descriptive and uses only fits with finite standard errors;
the per-parameter denominators are recorded in the summary CSV. For example,
the random-intercept variance coverage was 67.1%, 75.9%, 75.9% and 76.7% for
GRS-01 through GRS-04, respectively. These results do not validate nominal
95% coverage. No interval criterion was used to remove fitted point estimates
or to relabel the bias-MCSE results.

### Spatial-range identifiability diagnostic

As a post-fit diagnostic of the already registered runs (not a new acceptance
criterion), I examined spatial-scale estimates and Hessian standard errors in
the usable fits with finite full covariance uncertainty. Counts were 231,
257, 303 and 146 for GRS-01 through GRS-04. In that subset, the fitted-scale
10th percentile / median / 90th percentile was 0.049 / 0.313 / 1.686,
0.033 / 0.336 / 1.955, 0.039 / 0.893 / 2.699, and 0.042 / 1.201 / 2.586,
respectively, against truth 3 in every case. Median local standard errors for
the scale were 0.752, 0.834, 1.004 and 1.197, with 90th percentiles 3.645,
5.076, 3.799 and 3.951. The scale interval covered truth in only 35.1%,
42.4%, 49.8% and 59.6% of these fits. Across this same subset the replicate-
wise correlation of log fitted spatial variance with log fitted scale was
0.330, 0.316, 0.292 and 0.485. This modest positive association alone is not
a formal profile-likelihood or identifiability measure; together with the
large standard errors, low medians and poor coverage, it indicates that the
registered fitting design supplies weak and unstable range information.

This is a diagnostic from the existing Hessian outputs, not evidence that a
different estimator or design would pass. A proper follow-up would require a
separately specified profile or simulation study, with the current recovery
criteria left unchanged.

The median fit time ranged from 0.312 to 0.516 seconds on this single Linux
machine. These are diagnostic timings only. The two captured process sessions
took 8:07.58 and 12:24.84 (20:32.42 combined); their largest externally
measured peak RSS was 151,272 KiB. The first session completed and saved
GRS-01/02 before a summary-file column-name error stopped the script. The
runner was corrected, the saved results were resumed for GRS-03/04, and the
summary was regenerated independently. Missing Wald standard errors in the
raw file are marked as unavailable, not as coverage failures. The fit records
were not selectively discarded or regenerated.

## Reproducibility and next work

The package source archive used for the fits had SHA-256
`ba3a6432224cb7a851550fd9826e336655eef4c436f786ac5daf671c21d1163f`.
R 4.6.1 ran on x86_64 Ubuntu with BLAS
`/usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.0`. The raw-results and
summary SHA-256 hashes are, respectively,
`8d77fc78e8494addcb12b2d217c3f1a13c0d81e529f69eece726397d2f7269f7` and
`eb46a2a6ff4a450133020fe965f244c1fe65d6e818025078f1cd6dafd85ea02b`.
Recreate the mathematical checks and result table with:

```sh
Rscript tools/validate_grouped_slope_recovery.R \
  docs/validation/data/grouped_slope_recovery.csv
Rscript tools/summarize_grouped_slope_recovery.R \
  docs/validation/data/grouped_slope_recovery.csv \
  docs/validation/data/grouped_slope_recovery_summary.csv
```

Next, diagnose spatial-range identifiability and Hessian uncertainty under the
predeclared scenarios without changing their pass criteria. Any new data
design or estimator study must be separately specified and seeded before its
results are examined. The remaining TSC, MCP and NNP studies are still open.
