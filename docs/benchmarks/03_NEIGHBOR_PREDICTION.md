# Nearest-predecessor prediction diagnostic

## Purpose and limits

This small synthetic check measures approximation error for the new
bounded-memory joint-draw route. It is a development diagnostic, not a general
accuracy guarantee or a replacement for application-specific sensitivity
analysis.

## Design

- date: 16 October 2026;
- R 4.6.1 on the Linux development machine;
- 16 training locations on a 4 by 4 unit grid;
- 25 prediction locations on a 5 by 5 grid over the same domain;
- fixed Gaussian covariance: spatial variance 0.8, scale 1.8 and independent
  variance 0.15;
- latent-link target, 50,000 draws, four predecessors and draw seed 701;
- comparison against the complete fitted 25 by 25 predictive covariance;
- input-row and deterministic maximin orders evaluated separately.

## Results

| Ordering | Seconds | Relative Frobenius covariance error | Largest marginal-variance error | Regional-mean variance ratio |
|---|---:|---:|---:|---:|
| input | 0.327 | 0.04994 | 0.00390 | 1.08308 |
| maximin | 0.211 | 0.03927 | 0.00292 | 1.05152 |

The relative Frobenius error is

\[
 \|\widehat S-S\|_F/\|S\|_F,
\]

and the final column divides the empirical variance of the equally weighted
draw mean by its exact fitted-model value. Monte Carlo noise remains in every
empirical entry. Maximin ordering was somewhat closer in this one design, but
that does not prove it will dominate input or scientifically informed orders.

The subsequently implemented analytic covariance recursion removes Monte
Carlo noise from the comparison. On the identical four-predecessor design it
gave:

| Ordering | Relative Frobenius error | Maximum absolute covariance error | Maximum marginal-variance error | Regional-mean variance ratio |
|---|---:|---:|---:|---:|
| input | 0.044247 | 0.007160 | 0.0000776 | 1.075866 |
| maximin | 0.033065 | 0.006824 | 0.0000804 | 1.045203 |

These analytic values are the preferred approximation-error record. The
50,000-draw values above remain as an independent simulation check of the
recursion and of target aggregation.

The executable tests provide the stronger algebraic boundary: when every
predecessor is retained, same-seed sequential draws equal the ordinary exact
joint draws to numerical precision. With smaller neighborhoods, users should
vary both the order and neighbor count for the targets that matter to their
analysis.

## Direct target-streaming storage

The same 50,000-draw, 25-location design was aggregated to one equally
weighted target during sequential simulation. The returned target draw matrix
occupied 0.382 MiB, compared with 9.539 MiB for the location-level draw matrix.
Input ordering retained at most 10 point-draw vectors simultaneously and took
0.382 seconds; maximin ordering retained at most 16 and took 0.209 seconds.

These are observed active widths for this graph, not fixed bounds. Maximin can
improve covariance approximation while keeping early dispersed sites active
for longer. The implementation reports the peak so users can see that
trade-off. The target matrix itself scales with draws times targets and cannot
be discarded because it defines quantiles and exceedance probabilities.

