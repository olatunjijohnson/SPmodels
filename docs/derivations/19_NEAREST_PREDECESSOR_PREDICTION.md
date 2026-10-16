# Bounded-memory nearest-predecessor predictive draws

## Purpose and approximation boundary

Chapter 6 derives the fitted joint Gaussian predictive distribution

\[
 \boldsymbol\eta_*\mid\boldsymbol y
 \ \dot\sim\ N(\boldsymbol m,\boldsymbol S),                    \tag{NP-01}
\]

where the dot is needed for non-Gaussian fits because their latent prediction
already uses a Laplace approximation. Ordinary joint drawing constructs the
complete \(M\times M\) matrix \(\boldsymbol S\). That is exact for the declared
predictive distribution but requires quadratic grid memory.

`sp_predict_draws_neighbor()` instead defines a sparse sequential
factorization. It is exact when every predecessor is retained and otherwise a
clearly labelled, order-dependent approximation. It reduces prediction-grid
working memory; it does not make the dense training fit scalable and does not
include fitted-parameter uncertainty.

## From the Gaussian chain rule to local conditioning

Choose a permutation \(\pi=(\pi_1,\ldots,\pi_M)\) of prediction rows. The
probability chain rule is the identity

\[
 p(\eta_{\pi_1},\ldots,\eta_{\pi_M})=
 \prod_{j=1}^M p(\eta_{\pi_j}\mid
   \eta_{\pi_1},\ldots,\eta_{\pi_{j-1}}).                   \tag{NP-02}
\]

For row \(i=\pi_j\), let \(P_i=\{\pi_1,\ldots,\pi_{j-1}\}\). If all of
\(P_i\) is retained, the multivariate-normal conditioning identity gives

\[
 \eta_i\mid\boldsymbol\eta_{P_i}\sim N(a_i,v_i),               \tag{NP-03}
\]

\[
 a_i=m_i+\boldsymbol S_{iP_i}\boldsymbol S_{P_iP_i}^{-1}
 (\boldsymbol\eta_{P_i}-\boldsymbol m_{P_i}),                  \tag{NP-04}
\]

\[
 v_i=S_{ii}-\boldsymbol S_{iP_i}
 \boldsymbol S_{P_iP_i}^{-1}\boldsymbol S_{P_i i}.             \tag{NP-05}
\]

Multiplying (NP-03) over all rows recovers exactly the density in (NP-01).
Sequential simulation from these conditionals is therefore another
factorization of the same Gaussian distribution, analogous to Cholesky
simulation.

## Nearest-predecessor approximation

Fix a positive neighbor budget \(r\). For every row \(i\), SPmodels selects

\[
 N_i=\text{the at most }r\text{ members of }P_i
 \text{ having smallest Euclidean coordinate distance to }i.  \tag{NP-06}
\]

Ties retain predecessor order. It then replaces (NP-03) by

\[
 \widetilde p(\boldsymbol\eta_*)=
 \prod_{j=1}^M p(\eta_{\pi_j}\mid\boldsymbol\eta_{N_{\pi_j}}),\tag{NP-07}
\]

using (NP-04)--(NP-05) with \(P_i\) replaced by \(N_i\). For independent
standard-normal values \(z_i^{(b)}\), draw \(b\) is generated recursively as

\[
 \eta_i^{(b)}=\widetilde a_i^{(b)}+\sqrt{\widetilde v_i}\,
 z_i^{(b)}.                                                     \tag{NP-08}
\]

The predictive submatrix for \(N_i\cup\{i\}\) is calculated directly from
the same equations as the ordinary joint method. This retains spatial,
temporal and supported grouped covariance contributions inside each local
conditional. Cholesky solves are used in (NP-04)--(NP-05), not explicit matrix
inversion. Coincident latent prediction rows can make the predecessor block
positive semidefinite. If Cholesky fails, a checked symmetric eigendecomposition
uses the Moore--Penrose solution on eigenvalues above a scale-dependent
tolerance and verifies that the cross-covariance lies in the predecessor
column space. It does not add artificial jitter. A materially negative
eigenvalue, incompatible cross-covariance or conditional variance is an error;
only round-off-sized negative conditional variance is set to zero.

If \(r\geq M-1\), then \(N_i=P_i\) for every row and (NP-07) equals (NP-02).
This exact boundary is both a mathematical check and an executable test. For
smaller \(r\), omitted long-range conditional links change the joint
distribution. Marginal means are retained, but neither exact covariances nor
exact regional uncertainty should be claimed.

## Ordering, grouped effects and interpretation

The default permutation is input-row order. A caller may supply an explicit
permutation or request a deterministic maximin order. Maximin starts with input
row one. If \(A_{j-1}\) is the set already selected, the next row is

\[
 \pi_j=\mathop{\arg\max}_{i\notin A_{j-1}}
 \left\{\min_{k\in A_{j-1}}\|\boldsymbol c_i-
 \boldsymbol c_k\|_2\right\},                                  \tag{NP-06a}
\]

with input order breaking ties. The running minimum distance requires only
\(O(M)\) storage, although selecting the complete order uses \(O(M^2)\)
distance work. Returned draw columns always use original row order. Because
the neighbor sets in (NP-06) depend on the permutation, approximate results
can change with ordering. Sensitivity to defensible orders and neighbor counts
should be assessed for scientific analyses.

Neighbor selection uses spatial coordinate distance. A temporal or grouped
effect still enters the selected local covariance, but it does not change which
predecessors are selected. In particular, a far-away row sharing a grouped
level may be omitted. Increasing \(r\), choosing an informed order or using the
ordinary exact joint method is necessary when such long-range structure drives
the target.

Gaussian response/link draws use (NP-08) directly. Binomial-logit and
Poisson-log link draws use their Chapter 6 local Gaussian predictive
approximation, followed by the same inverse-link and optional response-sampling
steps as ordinary joint draws.

## Optional fitted-parameter uncertainty

Chapter 6 equation (J-15) defines the fitted-Hessian parameter draw
\(\boldsymbol\psi^{(b)}\). When `parameter_uncertainty = TRUE`, SPmodels
rebuilds the fitted covariance or Laplace state for each such draw. Denote its
predictive moments by \(\boldsymbol m^{(b)}\) and
\(\boldsymbol S^{(b)}\). Equations (NP-04)--(NP-08) are then evaluated using
those draw-specific moments:

\[
 \eta_i^{(b)}=m_i^{(b)}+
 S_{iN_i}^{(b)}\{S_{N_iN_i}^{(b)}\}^{-1}
 (\boldsymbol\eta_{N_i}^{(b)}-\boldsymbol m_{N_i}^{(b)})+
 \sqrt{\widetilde v_i^{(b)}}z_i^{(b)}.                         \tag{NP-09}
\]

Estimated intervention impact and decay rebuild the row multipliers for every
response-scale draw. The parameter draws and multiplier draws are retained.
The reported conditional variance is the average of the draw-specific local
conditional variances; no complete predictive covariance is manufactured.
When every predecessor is retained, this uses the same parameter and prediction
random-number order as ordinary `sp_predict_draws()`, providing a draw-for-draw
exact boundary. With a small neighbor budget it combines the local joint
approximation with the separate local-Hessian parameter approximation.

## Streaming grouped targets

When only grouped weighted targets are required, storing all \(BM\) point
draws is unnecessary. `sp_predict_target_neighbor()` normalizes the raw target
weights exactly as in Chapter 6 and updates

\[
 T_g^{(b)}\leftarrow T_g^{(b)}+w_{ig}D_i^{(b)}.               \tag{NP-10}
\]

as soon as point \(i\) is simulated. A point may contribute to several targets,
as occurs with fractional polygon-cell allocations.

Let \(\ell(i)\) be the last processing step at which point \(i\) appears in a
later neighbor set, or its own step if it is never reused. Its \(B\)-draw link
vector is deleted immediately after step \(\ell(i)\). If \(W\) is the largest
number of simultaneously active point vectors and \(G\) the number of targets,
draw storage is of order

\[
 O\{B(W+G)\},                                                  \tag{NP-11}
\]

rather than automatically \(O\{B(M+G)\}\). The neighbor graph itself uses
\(O(Mr)\) indices. The result records \(W\), the graph, order, conditional
variances and sparse normalized allocation rows. It does not claim that
\(W\) is always small: an early site used by many later conditionals remains
active, and retaining all predecessors can give \(W=M\).

Streaming currently conditions on fitted parameters. For observation-scale
binomial or Poisson targets, response sampling occurs when each point is
visited, so a fixed seed need not give the same individual response draws as a
post-hoc workflow that consumes all Gaussian random numbers first. The
declared distribution is unchanged; Gaussian and latent-response cases have
direct same-seed aggregation checks.

## Analytic implied covariance diagnostic

For a manageable validation subset, the approximation can be examined without
Monte Carlo noise. Define the local coefficient row

\[
 \boldsymbol A_i=\boldsymbol S_{iN_i}
 \boldsymbol S_{N_iN_i}^{-1}.                                 \tag{NP-12}
\]

After approximate covariance rows for all predecessors \(P_i\) have been
constructed, the recursion implied by (NP-07) is

\[
 \widetilde{\boldsymbol S}_{iP_i}=
 \boldsymbol A_i\widetilde{\boldsymbol S}_{N_iP_i},            \tag{NP-13}
\]

\[
 \widetilde S_{ii}=\widetilde v_i+
 \boldsymbol A_i\widetilde{\boldsymbol S}_{N_iN_i}
 \boldsymbol A_i^{\mathsf T}.                                 \tag{NP-14}
\]

`sp_neighbor_covariance()` applies (NP-12)--(NP-14) and can compare the result
with the ordinary exact fitted predictive covariance. It reports relative
Frobenius error, maximum absolute covariance error and maximum marginal-
variance error. This diagnostic intentionally allocates \(M^2\) entries; it is
for selecting/checking a budget on manageable subsets, not the production
bounded-memory route. With all predecessors, the recursion recovers the exact
fitted covariance.

## Storage and computation

At a single step, the largest newly constructed cross-covariance has order
\(n(r+1)\), and the local prediction covariance has order \((r+1)^2\), where
\(n\) is the number of training rows. The complete \(M^2\) prediction
covariance is never returned or allocated. The unavoidable returned draw
matrix contains \(BM\) values for \(B\) draws.

This is a memory result, not a speed guarantee. With parameter uncertainty,
the local conditionals are rebuilt for every one of the \(B\) parameter draws.
Repeated dense training solves can be expensive, and the original fitted model still required dense training
covariance factorization. Sparse or low-rank training likelihoods remain a
separate development and validation problem.

## Algorithm and validation map

1. Validate the spatial fit, draw count, neighbor budget and row permutation.
2. Build the prediction design once and generate all call-local standard
   normals.
3. Resolve input, maximin (NP-06a) or user-supplied order; visit rows in that
   order and find (NP-06).
4. Construct only the predictive distribution for \(N_i\cup\{i\}\).
5. Calculate (NP-04)--(NP-05) by Cholesky solves and draw with (NP-08), or
   rebuild the fitted state and use (NP-09) for every parameter draw.
6. Restore original column order, apply any response transformation, and
   return neighbor indices, conditional variances and order instead of a full
   covariance matrix.
7. For direct targets, apply (NP-10), release point vectors after their last
   graph use and report (NP-11)'s observed active width.
8. For an explicit covariance diagnostic, compose (NP-12)--(NP-14) and,
   optionally, compare it with the exact fitted covariance.

| Equations | Implementation | Independent checks |
|---|---|---|
| NP-01--NP-05 | existing joint predictive equations plus sequential loop | Same-seed equality to ordinary Gaussian and non-Gaussian link draws when all predecessors are retained |
| NP-06--NP-08 | `sp_predict_draws_neighbor()` | Neighbor-count bound, deterministic maximin/user ordering, reordered moments, exact full-predecessor boundary and reproducibility |
| grouped covariance | local design subsetting | Same-seed exact boundary with a shared unseen marginal group level |
| NP-09 | parameter-draw loop | Same-seed Gaussian and fitted-intervention equality to ordinary all-predecessor joint draws; retained parameter/multiplier draws |
| NP-10--NP-11 | `sp_predict_target_neighbor()` | Same-seed Gaussian and binomial-latent post-hoc target equality, repeated polygon allocation, sparse normalized weights, peak active width and caller-seed preservation |
| NP-12--NP-14 | `sp_neighbor_covariance()` | Exact all-predecessor Gaussian/non-Gaussian equality, PSD duplicate boundary, finite small-neighbor errors and symmetric positive-semidefinite result |
| failure and audit contracts | public validation and returned metadata | Invalid budget/order, local seed preservation, no full covariance and printed approximation label |
| singular predecessor boundary | `.sp_psd_solve()` | Coincident latent rows give identical draws and zero second conditional variance without jitter |

The executable checks are in `tests/test-neighbor-prediction.R`.

