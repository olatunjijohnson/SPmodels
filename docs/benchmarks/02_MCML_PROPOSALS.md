# MCML proposal diagnostic benchmark

## Purpose and boundary

This record examines how declared Gaussian proposal scales behave on the
authorised Loa loa teaching data. It motivated an opt-in, discarded-pilot
adaptation route. It is not a comparison with another package, not a universal
choice of proposal settings, and not evidence that adaptation always helps.

## Data and environment

- date: 16 October 2026;
- SPmodels development source under R 4.6.1 on Ubuntu 24.04.4;
- data: `loaloa` from the owner's local `lsda` teaching data, 197 rows;
- preliminary model: binomial-logit Laplace spatial fit;
- the data object was loaded read-only and was not copied into SPmodels;
- every reported MCML run used a local seed and recorded ESS fraction,
  relative likelihood MCSE and maximum normalized weight.

The preliminary Laplace fit took 15.638 seconds and converged. These timings
are single-machine observations, so the diagnostic values—not small timing
differences—are the relevant evidence.

## Declared-scale results

The first runs used 1,000 retained importance draws with different seeds.

| Proposal scales | Seed | Seconds | ESS fraction | Relative MCSE | Maximum weight |
|---|---:|---:|---:|---:|---:|
| 1 | 811 | 45.738 | 0.230324 | 0.057837 | 0.020501 |
| 0.75, 1.5 | 812 | 57.131 | 0.002725 | 0.604766 | 0.571487 |
| 0.6, 1, 1.8 | 813 | 34.901 | 0.008242 | 0.344509 | 0.260051 |

Because the seeds differ, this first table is an operational stress check and
not a paired causal comparison. A second check held the seed and retained draw
count fixed at seed 900 and 900 draws.

| Proposal scales | Seconds | Convergence | ESS fraction | Relative MCSE | Maximum weight | MC log likelihood |
|---|---:|---:|---:|---:|---:|---:|
| 1 | 30.631 | 0 | 0.020199 | 0.220352 | 0.146239 | -664.979796 |
| 0.9, 1.1 | 84.105 | 0 | 0.001022 | 0.988961 | 0.988971 | -662.249825 |

The wider mixture was markedly less stable in this same-seed check: nearly all
normalized weight was carried by one draw. A normalized mixture is a valid
importance distribution, but validity does not imply low variance.

## Development consequence

`sp_mcml_control(adapt = TRUE)` now offers a separate pilot stage. It computes
importance-weighted pilot moments, shrinks them toward the positive-definite
Laplace moments, discards the pilot from likelihood estimation, and draws a
new retained sample from the complete normalized adapted proposal. Both pilot
and main diagnostics are retained. The derivation explains why the estimator
continues to target the marginal likelihood conditional on the pilot.

The benchmark does not select adaptation as the default. Proposal performance
must be judged from independent repeats, ESS, maximum weight and relative
MCSE. A low-ESS pilot can itself be unstable, and the implementation warns
and retains the original proposal moments when the pilot fails the declared
ESS threshold. The retained status does not silently call such an attempt an
adaptation success.

## Same-seed adaptive sensitivity

An additional 900-pilot/900-main run used the same preliminary fit, seed 900,
scale one and shrinkage 0.50. Its pilot ESS fraction was 0.057079, relative
MCSE was 0.135557 and maximum weight was 0.118371. The independently drawn main
sample was worse: ESS fraction 0.002439, relative MCSE 0.674471 and maximum
weight 0.645035. The optimizer converged, but that convergence does not repair
the weak importance sample. These sensitivity runs deliberately used
`ess_warning = 0.001`, below the measured pilot ESS, to exercise the adaptation
rather than the new default fallback gate.

Increasing shrinkage did not rescue this particular random stream. At
shrinkage 0.80, main ESS fraction was 0.001140 and the maximum weight was
0.987232; at shrinkage 0.95 they were 0.001190 and 0.965787. These main draws
occur after the pilot, so they are not draw-for-draw paired with the earlier
non-adaptive run. They nevertheless demonstrate why the feature stays off by
default and why both pilot and main diagnostics are exposed in `summary()`.

