# Nearest-predecessor prediction study results

## Scope and status

The preregistered NNP study was executed on 16 October 2026. It evaluates the
bounded-memory nearest-predecessor joint-draw approximation against the
ordinary exact fitted joint predictive distribution. It is a synthetic
accuracy and resource study, not a numerical comparison with the separately
managed reference implementation and not a general performance guarantee.

The study contains five cases (Gaussian, binomial-logit, Poisson-log,
separable space-time Gaussian, and Gaussian with a fixed grouped random-slope
covariance), prediction grids of 25, 100 and 225 rows, two orderings, four
registered neighbor budgets and five independent response seeds. Each
configuration uses 20,000 draws. All 3,000 retained target rows (600
configurations × five targets) completed; there were no fit/prediction errors
or warnings. All criteria and synthetic fixtures were fixed before results
were examined in
[`02_RECOVERY_AND_APPROXIMATION_DESIGN.md`](02_RECOVERY_AND_APPROXIMATION_DESIGN.md).

The checksummed SPmodels source archive was
`4bd9544a4ab011b08f5d7786c1d0db526f903ecc8cc93f08b23844f71247a7e7`. Runs
used R 4.6.1 on Ubuntu 24.04.4, `x86_64-pc-linux-gnu`, UTF-8, and BLAS
`/usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.0`. Timing and memory
observations are specific to this Linux host.

## Case-specific results

The table gives the smallest budget for which all five predeclared targets
passed the relative-SD and threshold-probability criteria in both input and
maximin order. The word “candidate” is deliberate: these are results for the
registered synthetic designs, fitted states, thresholds and targets only.

| Case | 25 rows | 100 rows | 225 rows |
| --- | ---: | ---: | ---: |
| Gaussian | 5 | 10 | 20 |
| Binomial-logit | 5 | 10 | 20 |
| Poisson-log | 5 | 10 | 20 |
| Fixed grouped random slope | 5 | 10 | 20 |
| Separable space-time Gaussian | 10 | 40 | No registered budget passed |

No universal neighbor budget is recommended. At 225 rows in the space-time
case, even 40 predecessors failed the overall target's relative-SD criterion
for input ordering (8.2% error); the threshold-probability interval also
remained inconclusive for some targets. Spatial-only neighbor selection does
not guarantee retention of the dependence structure needed by larger joint
targets. The exact joint-draw route remains the choice when these registered
accuracy criteria are not met.

For the four non-temporal cases, the registered candidate budgets passed all
five target checks in both orderings. This does not establish performance for
other covariance parameters, grids, orderings, target definitions or optional
parameter draws. The covariance components were fixed; the exact predictive
covariance still includes the fixed-effect uncertainty from each model's
ordinary joint prediction route.

The full input-order all-predecessor boundary reproduced the ordinary draws
with the same seed for every case. Its maximum relative Frobenius covariance
error was `1.55e-15`, and its maximum absolute marginal-variance error was
`2.23e-16`. The maximin all-predecessor factorization reproduced the exact
covariance to the same numerical precision; same-seed draw equality is
recorded only for input order because the permutation changes the random-number
assignment.

Order sensitivity was material in some configurations: across the entire
registered matrix, the maximum absolute input-versus-maximin difference was
0.0209 for target means, 0.1278 for target SDs and 0.0502 for threshold
probabilities. These are maxima over cases, grid sizes, budgets, seeds and
targets, not typical differences. They reinforce the need to check ordering
when an approximate budget is used.

## Resource observations

The streamed target output contains five target draws. At 20,000 draws, its
matrix occupied about 0.76 MiB for each grid size. The corresponding complete
point-draw matrices occupied about 3.82, 15.26 and 34.35 MiB at 25, 100 and
225 rows. Streamed target matrices therefore scale with the five requested
targets rather than the complete prediction grid. The largest observed
simultaneous active-point width was 155 of 225 locations; maximin ordering
often kept more vectors live than input ordering.

Per-configuration streamed prediction elapsed time ranged from 0.10 to 21.41
seconds, with median 0.80 seconds across this study. The longest observation
is retained rather than discarded. These timings are not cross-platform
comparisons, and the run did not measure external peak RSS separately for
each NNP configuration.

## Reproducibility and interpretation

The raw per-case files retain every seed/ordering/budget/target row, warning,
error, covariance discrepancy, target error, time and active-width measure:

- [`neighbor_prediction_gaussian.csv`](data/neighbor_prediction_gaussian.csv)
- [`neighbor_prediction_binomial.csv`](data/neighbor_prediction_binomial.csv)
- [`neighbor_prediction_poisson.csv`](data/neighbor_prediction_poisson.csv)
- [`neighbor_prediction_spacetime.csv`](data/neighbor_prediction_spacetime.csv)
- [`neighbor_prediction_grouped_slope.csv`](data/neighbor_prediction_grouped_slope.csv)

The combined target summary, case-specific budget decisions and order-paired
sensitivity records are in the corresponding `neighbor_prediction_study_*`
CSV files in `data/`. The registered validator checked all expected rows,
seeds, unique configurations and exact-boundary outcomes. Reproduce or
revalidate with:

```sh
Rscript tools/validate_neighbor_prediction_study.R \
  docs/validation/data/neighbor_prediction_gaussian.csv \
  docs/validation/data/neighbor_prediction_binomial.csv \
  docs/validation/data/neighbor_prediction_poisson.csv \
  docs/validation/data/neighbor_prediction_spacetime.csv \
  docs/validation/data/neighbor_prediction_grouped_slope.csv

Rscript tools/summarize_neighbor_prediction_study.R \
  --output-prefix=docs/validation/data/neighbor_prediction_study \
  docs/validation/data/neighbor_prediction_gaussian.csv \
  docs/validation/data/neighbor_prediction_binomial.csv \
  docs/validation/data/neighbor_prediction_poisson.csv \
  docs/validation/data/neighbor_prediction_spacetime.csv \
  docs/validation/data/neighbor_prediction_grouped_slope.csv
```

The results establish only internal accuracy against the exact predictive
distribution produced by the same fitted SPmodels object. They do not
establish reference agreement, independent scientific validity of the fitted
models, macOS behavior or a generally safe production budget.
