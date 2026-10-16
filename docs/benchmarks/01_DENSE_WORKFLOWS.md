# Dense-workflow baseline benchmark

## Purpose and boundary

This is a reproducible first timing baseline for the current dense numerical
implementation. It identifies scaling pressure before optimization; it is not
a cross-platform performance guarantee, a competitor comparison or a claim
about production-sized data.

## Environment

- date: 15 October 2026;
- SPmodels: 0.0.0.9000 installed source build;
- R: 4.6.1 on Ubuntu 24.04.4 LTS, x86-64;
- CPU: Intel Core i7-9850H at 2.60 GHz, 6 physical/12 logical cores;
- memory: 62 GiB total, approximately 14 GiB available during the run;
- execution: one cold elapsed-time measurement per case, ordinary R process.

Single measurements are deliberately reported to only three decimals and
should be rerun with repeated warm timings before fine performance decisions.

## Synthetic designs

The Gaussian cases used 10 by 10 and 15 by 15 unit-spaced grids, one linear
predictor, a Gaussian response generated with seed `1000 + n`, and estimated
spatial variance, scale and nugget variance from starts 0.5, 2 and 0.25. The
binomial case used a 10 by 10 grid, 100 Bernoulli observations generated with
seed 1200, and fixed covariance values 0.5, 2 and 0.3. All distances were
labelled kilometres.

The joint-draw case produced 100 fitted-location Gaussian draws over 100
locations. The validation case performed three complete estimated Gaussian
training-fold refits at the 100-row design.

## Results

| Case | Elapsed seconds | Convergence |
|---|---:|---:|
| Estimated Gaussian fit, n = 100 | 1.454 | 0 |
| Estimated Gaussian fit, n = 225 | 40.792 | 0 |
| Fixed-covariance binomial Laplace fit, n = 100 | 0.631 | 0 |
| 100 joint Gaussian draws at 100 locations | 0.020 | not applicable |
| Three-fold estimated Gaussian cross-validation, n = 100 | 2.063 | all folds completed |

## Interpretation

The 100-to-225 Gaussian increase is much larger than the row-count increase.
That is consistent with repeated dense covariance factorizations inside
optimization, but one timing pair cannot estimate a reliable complexity
exponent. Dense fitting at a few hundred locations is therefore the first
profiling target. Joint drawing at 100 locations is not the immediate
bottleneck in this baseline.

Before changing algorithms, rerun each case repeatedly under controlled CPU
load and profile the estimated n = 225 fit by function. Any sparse,
low-rank or compiled alternative must retain the existing small-case algebra
and validation tolerances; faster output is not acceptable if it changes the
declared model silently.

## Profile-guided half-integer Matérn optimization

A second n = 225 run under `Rprof(interval = 0.01)` sampled 8.18 seconds. Of
that sampled time, 84.7% was under `matern_cor()` and 65.5% was self-time in
`besselK()`. The model used the default smoothness 0.5, whose exact correlation
is \(\exp(-d/\phi)\). SPmodels therefore now dispatches smoothness 0.5, 1.5
and 2.5 to their already derived closed forms and retains the scaled-Bessel
path for every other positive smoothness.

In a same-process A/B kernel check on the 225 by 225 distance matrix, 50
evaluations had median elapsed time 0.231 seconds through the new 0.5 fast path
and 1.470 seconds through the previous generic scaled-Bessel expression: a
6.36-fold kernel speedup. Maximum absolute numerical difference was
`1.11e-16`. A complete estimated n = 225 fit then took 2.914 seconds and
converged with code 0. Because the earlier 40.792-second cold run and this warm
run experienced different system conditions, their ratio is not presented as
a controlled whole-fit speedup.

A subsequent n = 400 profile separated the remaining costs more clearly. The
pre-change fit took 6.363 seconds and converged with code 0; sampled self-time
was 43.7% in dense Cholesky factorization and 27.0% in Matérn construction.
The internal covariance builder was therefore changed to call an
already-validated Matérn kernel: it retains the same closed forms and generic
Bessel fallback, but does not rescan the complete distance matrix at every
optimizer evaluation. A separate post-change fit took 5.175 seconds and also
converged with code 0. These are separate-process observations, not a
controlled whole-fit speedup claim.

The controlled same-process check evaluated each kernel 100 times on the same
400 by 400 distance matrix in five repetitions. Median elapsed time was 1.083
seconds for the validating public function and 0.477 seconds for the internal
already-validated route, a 2.27-fold kernel speedup, with maximum absolute
difference zero. Direct tests cover all three specialized smoothness values
and generic smoothness values, including matrix dimensions and labels. Dense
Cholesky factorization is now the principal measured bottleneck at this size.
This evidence supports retaining exact dense fitting for current cases; it
does not justify silently substituting an approximation.

## Repeated post-optimization Linux timings

A later controlled local run used square grids of 100, 225 and 400 rows, the
same formula structure and explicit covariance starts, and three consecutive
fits of each unchanged dataset in one R process. The response for each size
used seed `1000 + n`; every repetition converged with code zero and reproduced
its size-specific log likelihood exactly at the printed precision.

| Rows | Elapsed seconds (three repetitions) | Median seconds |
|---:|---:|---:|
| 100 | 0.896, 0.233, 0.176 | 0.233 |
| 225 | 1.419, 1.415, 1.432 | 1.419 |
| 400 | 11.888, 10.860, 12.183 | 11.888 |

The first 100-row observation includes visible cold-start cost; the larger
cases are much more stable within size. These measurements reinforce the
profile result: exact dense fitting becomes expensive by 400 observations and
does not scale like the row count alone. They are a Linux development-machine
baseline, not a cross-platform performance guarantee. Broader platform and
real-consumer timings remain necessary before selecting a compiled, sparse or
low-rank replacement.

## Bounded-memory point-grid path

The prediction audit identified a separate allocation that scales as the
number of fitted rows times the number of grid rows. `sp_predict_chunked()` now
evaluates exact pointwise means and marginal standard errors in chunks of at
most $c$ rows. The largest cross-covariance allocation is thereby reduced
from $nm$ to at most $nc$ doubles. For example, $n=1{,}000$,
$m=100{,}000$ requires about 800 MB for one unchunked double matrix before
temporary calculations, whereas $c=1{,}000$ reduces that matrix to about
8 MB.

This is a storage bound, not a timing claim. All $nm$ distances are still
evaluated and dense fitting remains unchanged. Direct tests with chunk sizes
that do not divide the grid length reproduce one-shot Gaussian and binomial
means and standard errors, including correlated grouped slopes and marginal
unseen levels. Joint covariance across grid rows is intentionally unavailable
from this route.

