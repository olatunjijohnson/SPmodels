# Gaussian nearest-predecessor likelihood diagnostic

## Scope

This is a reproducible Linux development diagnostic for the first
nearest-predecessor training likelihood. It is not a cross-platform speed
claim and does not select a default neighbor count. The public dense fitter
remains unchanged; `sp_fit_neighbor()` is explicit opt-in.

The run used R 4.6.1 on Ubuntu 24.04.4, seed 460001, a 15 by 15 regular grid
(`n = 225`), one standardized linear covariate, exponential Matérn correlation
(`smoothness = 0.5`), spatial variance 0.8, scale 2.5 km and nugget variance
0.2. Covariance values were fixed so the diagnostic isolates one likelihood
evaluation plus graph construction. After one warm-up per configuration, five
complete fits were timed in one R process. Object sizes are retained fitted-
object sizes, not peak-memory measurements.

## Recorded results

All approximate fits used deterministic maximin ordering.

| Predecessors | Median elapsed (s) | Retained object (bytes) | Log-likelihood difference from dense | Maximum coefficient difference |
| ---: | ---: | ---: | ---: | ---: |
| Dense | 0.009 | 1,006,752 | 0 | 0 |
| 5 | 0.113 | 188,192 | 0.3605825 | 0.05650961 |
| 10 | 0.114 | 208,960 | 0.7434373 | 0.03136671 |
| 20 | 0.125 | 232,272 | 0.04392653 | 0.01861897 |
| 224 (all) | 0.579 | 474,048 | -2.84e-14 | 3.33e-16 |

The five dense timings were 0.062, 0.009, 0.009, 0.008 and 0.008 seconds.
Corresponding neighbor-fit timings were 0.109--0.132 seconds for 5,
0.109--0.120 for 10, 0.121--0.145 for 20, and 0.557--0.624 for all
predecessors.

An additional estimated-covariance diagnostic used seed 460101 and 100 random
locations. Every optimizer returned convergence code zero. Dense ML took
0.663 seconds. Neighbor fits took 3.978, 4.603 and 9.967 seconds for budgets
10, 20 and 99. The all-predecessor covariance estimates, log likelihood and
coefficients agreed with dense ML to at most `6.89e-13`. Budgets 10 and 20 had
maximum coefficient differences 0.04138 and 0.01339 from the dense fit.

## Separate-process peak-memory diagnostic

`/usr/bin/time -v` measured separate fresh R processes using fixed covariance,
input ordering, budget 20 and deterministic standard-normal responses. This is
a computational stress diagnostic rather than a model-recovery experiment;
the same source loading overhead is present in each process.

| Training rows | Route | Elapsed (s) | Maximum resident set (KiB) | Retained object (bytes) |
| ---: | --- | ---: | ---: | ---: |
| 1,000 | dense | 1.14 | 131,960 | 16,643,152 |
| 1,000 | neighbor 20 | 1.13 | 91,784 | 824,912 |
| 2,000 | dense | 3.08 | 300,132 | 65,219,152 |
| 2,000 | neighbor 20 | 2.37 | 130,496 | 1,588,912 |
| 5,000 | neighbor 20 | 6.25 | 134,112 | 3,880,912 |

The 5,000-row dense route was not run merely to manufacture a dramatic memory
contrast. The smaller runs already show the expected dense growth, while the
neighbor route completed the larger case without an `n` by `n` covariance.
These are single runs on one machine and are not substitutes for the
predeclared repeated cross-platform study.

## Predeclared recovery and stability study

The statistical-accuracy and budget-stability portion of the NNL design in
`docs/validation/02_RECOVERY_AND_APPROXIMATION_DESIGN.md` was run for 25
response replicates at each of three training sizes
(100, 225, 400), three Matérn scales (0.05, 0.20, 0.50 times the unit-square
diameter), two orders and five budgets (5, 10, 20, 40, all predecessors).
This produced 225 independent datasets and 2,250 neighbor fits. Covariance
parameters were fixed at their known generating values to isolate the
approximation; each dataset included 50 held-out coordinates generated jointly
with the training response. The checksummed source archive was
`ead290f443dd1deb433ac7684e4bdcfed5dddb362135c44ccafb4022fc14f294`, under R
4.6.1 on x86_64 Ubuntu. All 225 datasets and fits completed; no failed fit was
excluded. This first run did not capture timing or memory, so a second complete
run captured per-fit elapsed times and externally measured peak resident memory
for the complete study process. That performance run used archive
`e67e3d3f31cb78e43ba20b274c5cee400f701060611b5975597aa4955a36a431`; all
accuracy outcomes matched the first run across 2,250 rows. GNU time measured
23:08.31 elapsed and 233,124 KiB maximum RSS for the full run. Per-fit timings
are recorded at `docs/benchmarks/data/neighbor_likelihood_performance_study.csv`
and summarized in
`docs/benchmarks/data/neighbor_likelihood_performance_summary.csv`. Their
SHA-256 hashes are respectively
`a27d92958ad2a37183ca6877fc228766d6ea8a83a56afe06ac7fb5b0914714f4` and
`f8079ef1b841f8771508eaf0fc3fcb804fb656c7f6d94fd7ca28c39b022cf6a9`. The single
process-level RSS is not a per-configuration peak-memory comparison; the
separate-process diagnostics above remain the only route-specific external
memory measurements. Median `sp_fit_neighbor()` elapsed time, pooled across
range, ordering and replicate settings, was:

| Training rows | 5 neighbors (s) | 20 neighbors (s) | All predecessors (s) |
| ---: | ---: | ---: | ---: |
| 100 | 0.047 | 0.051 | 0.077 |
| 225 | 0.119 | 0.133 | 0.632 |
| 400 | 0.219 | 0.248 | 3.801 |

These timings are from one Ubuntu process and describe fit calls only; they
exclude dense-baseline fitting and data generation. Raw accuracy results are in
`docs/benchmarks/data/neighbor_likelihood_study.csv`, and the reproducible
summary is `docs/benchmarks/data/neighbor_likelihood_summary.csv`. Recompute the
summary with `Rscript tools/summarize_neighbor_likelihood_study.R`.

A replicate passed when both regression coefficients differed from the dense
fit by less than 0.05 dense-fit standard errors and held-out prediction RMSE
differed by no more than 2% on the log ratio scale. A configuration passed a
scenario when at least 90% of its replicates passed. The smallest tested
budget passing under both orders was:

| Training rows | Scale fraction | Smallest passing budget | Result |
| ---: | ---: | ---: | --- |
| 100 | 0.05 | 10 | Scenario-specific budget passed |
| 225 | 0.05 | 20 | Scenario-specific budget passed |
| 400 | 0.05 | 20 | Scenario-specific budget passed |
| 100 | 0.20 | 40 | Scenario-specific budget passed |
| 225 | 0.20 | 224 | Dense boundary required |
| 400 | 0.20 | 399 | Dense boundary required |
| 100 | 0.50 | 40 | Scenario-specific budget passed |
| 225 | 0.50 | 224 | Dense boundary required |
| 400 | 0.50 | 399 | Dense boundary required |

At those scenario-specific settings, both orders met the pass rule in 23--25
of 25 replicates (92--100%). At every exact boundary, all 25 replicates passed
for both orders. Maximum absolute exact-boundary discrepancies across the full
study were `5.12e-13` for log likelihood, `5.58e-14` for the largest
coefficient error in dense-fit standard errors, and `3.00e-14` for the
log-RMSE ratio.

These known-range simulation results do not support one general budget for
unknown data. At the two larger sample sizes, the medium- and long-range cases
needed every predecessor to satisfy the thresholds. No default or broad
production recommendation is made. The declared simulation protocol and exact
statistical-accuracy, per-fit timing and exact-boundary checks passed, and
whole-run peak RSS was measured. Per-configuration memory, reference-
comparison, and cross-platform gates remain open. Approximate-budget selection
remains case-specific.

## Interpretation

The exact boundary is confirmed. Small budgets greatly reduce the retained
object size and the algorithm never constructs an `n` by `n` covariance
matrix, but the present R loops are slower than optimized dense linear algebra
at these moderate sizes. Approximation error is not monotone in the displayed
log likelihood. Therefore:

- no neighbor count is declared generally adequate for unknown data;
- the implementation is described as memory-scalable, not faster;
- compiled/indexed neighbor search and local factorization are future
  performance work; and
- case-specific validation is required before any teaching or production use.

