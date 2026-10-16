\newpage

# Gaussian nearest-predecessor likelihood

## Purpose and evidence boundary

Dense Gaussian fitting factorizes an \(n\)-by-\(n\) covariance matrix. Its
storage grows as \(O(n^2)\), and a general Cholesky factorization takes
\(O(n^3)\) arithmetic. This chapter derives the explicit approximation used by
`sp_fit_neighbor()`. It replaces the dense joint density by ordered local
conditional densities. It does not change `sp_fit()`, and it is not claimed to
be exact unless every predecessor is retained.

The initial implementation supports a Gaussian identity-link point model with
one stationary isotropic Matérn field and independent nugget noise. Grouped,
random-slope, temporal and non-Gaussian likelihoods require separate
derivations and validation and are rejected rather than approximated silently.

The ordered conditional approximation follows the statistical construction of
Vecchia (1988). The implementation and equations below were derived directly
for SPmodels; no external software implementation was used. Guinness (2018)
explains why permutation is part of the approximation and studies ordering and
grouping strategies. SPmodels presently implements the stated input and
deterministic maximin orders, not the paper's broader grouping methods.

## Model and assumptions

For observations \(i=1,\ldots,n\), let

\[
 y_i=\boldsymbol x_i^{\mathsf T}\boldsymbol\beta+o_i+w_i+e_i,
 \qquad e_i\overset{ind}{\sim}N(0,\tau^2),                 \tag{NL-01}
\]

where \(o_i\) is a known offset and

\[
 \operatorname{Cov}(w_i,w_j)=
 \sigma^2\rho_\nu(\|\boldsymbol s_i-\boldsymbol s_j\|;\phi).
                                                                  \tag{NL-02}
\]

The coordinates are projected, both axes use the recorded distance unit,
\(\sigma^2>0\), \(\phi>0\), \(\tau^2>0\), and the smoothness \(\nu>0\) is
fixed. The mean design has full column rank. Write
\(\boldsymbol y_o=\boldsymbol y-\boldsymbol o\) and let
\(V_{ij}=\operatorname{Cov}(y_i,y_j)\). Thus

\[
 V_{ij}=\sigma^2\rho_\nu(d_{ij};\phi)+
 \tau^2\,\mathbb 1(i=j).                                   \tag{NL-03}
\]

## Ordered conditional factorization

Choose a permutation \(\pi=(\pi_1,\ldots,\pi_n)\). The exact density always
factorizes as

\[
 p(\boldsymbol y)=p(y_{\pi_1})
 \prod_{r=2}^{n}p(y_{\pi_r}\mid
 y_{\pi_1},\ldots,y_{\pi_{r-1}}).                         \tag{NL-04}
\]

For row \(i=\pi_r\), let \(P_i=\{\pi_1,\ldots,\pi_{r-1}\}\). Retain the at
most \(m\) members of \(P_i\) nearest to \(\boldsymbol s_i\); call that set
\(N_i\). The approximation is

\[
 \widetilde p(\boldsymbol y)=
 \prod_{r=1}^{n}p(y_{\pi_r}\mid\boldsymbol y_{N_{\pi_r}}).
                                                                  \tag{NL-05}
\]

For input ordering, \(\pi_r=r\). Deterministic maximin ordering starts with
row one and repeatedly selects the unselected location farthest from the
selected set. A user-supplied permutation is also permitted. Order and graph
are retained because they are part of the approximation.

## One local Gaussian conditional

For a current observation \(i\), define

\[
 \boldsymbol a_i=V_{iN_i}V_{N_iN_i}^{-1},\qquad
 d_i=V_{ii}-\boldsymbol a_iV_{N_ii}.                         \tag{NL-06}
\]

Standard multivariate-normal conditioning gives

\[
 y_{o,i}\mid\boldsymbol y_{o,N_i}\sim N\left(
 \boldsymbol x_i^{\mathsf T}\boldsymbol\beta+
 \boldsymbol a_i(\boldsymbol y_{o,N_i}-
                  X_{N_i}\boldsymbol\beta),\ d_i\right).   \tag{NL-07}
\]

No inverse is formed in code. The system in (NL-06) is solved by Cholesky
factorization, with the checked positive-semidefinite solve used for an exact
singular boundary. A non-positive conditional variance fails the candidate
parameter state.

Define transformed data

\[
 \widetilde y_i=y_{o,i}-\boldsymbol a_i\boldsymbol y_{o,N_i},
 \qquad
 \widetilde{\boldsymbol x}_i^{\mathsf T}=
 \boldsymbol x_i^{\mathsf T}-\boldsymbol a_iX_{N_i}.         \tag{NL-08}
\]

Equation (NL-07) is now an independent weighted regression statement under the
approximate density:

\[
 \widetilde y_i\sim N(
 \widetilde{\boldsymbol x}_i^{\mathsf T}\boldsymbol\beta,d_i).
                                                                  \tag{NL-09}
\]

## Profiled mean and approximate log likelihood

Let \(D=\operatorname{diag}(d_1,\ldots,d_n)\), stack the transformed rows in
\(\widetilde X\), and stack the responses in \(\widetilde{\boldsymbol y}\).
For fixed covariance parameters \(\boldsymbol\theta=(\sigma^2,\phi,\tau^2)\),
the generalized least-squares estimate is

\[
 \widehat{\boldsymbol\beta}(\boldsymbol\theta)=
 (\widetilde X^{\mathsf T}D^{-1}\widetilde X)^{-1}
 \widetilde X^{\mathsf T}D^{-1}\widetilde{\boldsymbol y}. \tag{NL-10}
\]

Substitution gives the profiled approximate log likelihood

\[
 \widetilde\ell(\boldsymbol\theta)=
 -\frac12\sum_{i=1}^{n}\left\{
 \log(2\pi)+\log d_i+
 \frac{(\widetilde y_i-
 \widetilde{\boldsymbol x}_i^{\mathsf T}
 \widehat{\boldsymbol\beta})^2}{d_i}\right\}.              \tag{NL-11}
\]

The positive covariance parameters are optimized on their log scales with
bounded L-BFGS-B. The conditional coefficient covariance reported at the final
state is

\[
 \widehat{\operatorname{Var}}(\widehat{\boldsymbol\beta}\mid
 \widehat{\boldsymbol\theta})=
 (\widetilde X^{\mathsf T}D^{-1}\widetilde X)^{-1}.          \tag{NL-12}
\]

When the optimizer Hessian in log covariance parameters is positive definite,
its inverse is transformed to the natural scale by the delta method. This is a
local curvature approximation; it does not include finite-neighbor model
error.

## Exact all-predecessor boundary

If \(m\ge n-1\), then \(N_i=P_i\) for every row. Equation (NL-05) is exactly
(NL-04), whatever valid order is used. Consequently (NL-10)--(NL-12) equal
the dense Gaussian ML log likelihood, profiled regression estimate and
conditional coefficient covariance, apart from floating-point roundoff. This
is an executable boundary check, not an empirical expectation.

For \(m<n-1\), the output is an approximation. Increasing \(m\) need not make
every finite-sample parameter estimate move monotonically, so selection must
use a predeclared stability study rather than choosing the most convenient
result.

## Local prediction

For a new location \(0\), select its at most \(m\) nearest training
observations \(N_0\). With

\[
 \boldsymbol a_0=V_{0N_0}V_{N_0N_0}^{-1},                  \tag{NL-13}
\]

the plug-in conditional mean is

\[
 \widehat y_0=\boldsymbol x_0^{\mathsf T}\widehat{\boldsymbol\beta}+o_0+
 \boldsymbol a_0(\boldsymbol y_{o,N_0}-
 X_{N_0}\widehat{\boldsymbol\beta}).                       \tag{NL-14}
\]

Let \(\boldsymbol g_0^{\mathsf T}=\boldsymbol x_0^{\mathsf T}-
\boldsymbol a_0X_{N_0}\). The approximate latent prediction variance is

\[
 \sigma^2-V_{0N_0}V_{N_0N_0}^{-1}V_{N_00}+
 \boldsymbol g_0^{\mathsf T}
 \widehat{\operatorname{Var}}(\widehat{\boldsymbol\beta})
 \boldsymbol g_0.                                           \tag{NL-15}
\]

Add \(\tau^2\) for a future observation. Each point is conditioned locally;
the method does not return or imply a joint prediction covariance. Use the
separate nearest-predecessor draw and covariance-audit functions when joint
approximate prediction is required.

At a training row, the same calculation uses that row's recorded predecessor
set. Its mean-estimation gradient is the corresponding transformed design row
from (NL-08); this contribution is included in the reported sequential fitted-
value standard error. Consequently the observation and latent variances differ
by exactly \(\tau^2\) at a common row.

## Storage, computation and limitations

The retained graph uses \(O(nm)\) indices. Each likelihood contribution
factorizes at most an \(m\)-by-\(m\) covariance block, so covariance workspace
is \(O(m^2)\) rather than \(O(n^2)\). The direct implementation takes up to
\(O(nm^3)\) factorization work per objective evaluation; neighbor discovery is
currently an exact distance scan with \(O(n^2)\) time but bounded working
memory. This is a memory-scalable statistical foundation, not yet a compiled
or indexed speed implementation.

The approximation does not remove the need for recovery studies, neighbor-
budget sensitivity checks, cross-platform timings and comparison against a
separately frozen reference environment. Non-Gaussian, grouped and temporal
extensions remain unsupported until their own derivations and tests exist.

## Algorithm and validation map

1. Build the model frame, projected coordinates and full-rank mean design.
2. Resolve and retain input, maximin or supplied order and at most \(m\)
   nearest predecessors per row.
3. For each covariance candidate, compute (NL-06), transform with (NL-08), and
   profile the mean by (NL-10).
4. Optimize (NL-11) on log covariance scales or evaluate supplied fixed values.
5. Retain the graph, local coefficients, conditional variances, optimizer and
   unrounded estimates.
6. Predict pointwise with (NL-13)--(NL-15), explicitly labelling the result as
   nearest-predecessor approximate.

| Equations | Implementation | Independent checks |
| --- | --- | --- |
| NL-01--NL-03 | existing Matérn and nugget covariance | fixed-parameter simulated fixture |
| NL-04--NL-09 | local graph and Gaussian-state builders | graph budget/order, positive local variances, caller RNG unchanged |
| NL-10--NL-12 | profiled local GLS and optimizer | all-predecessor dense likelihood, coefficient and covariance equality; finite estimated state |
| NL-13--NL-15 | neighbor-fit prediction method | finite latent/observation predictions and observation variance ordering |
| failure contract | `sp_fit_neighbor()` validation | invalid budget/order, coincident-only coordinates and unsupported grouped covariance fail explicitly |

The executable checks are in `tests/test-neighbor-fit.R`.

## References

Vecchia, A. V. (1988). Estimation and model identification for continuous
spatial processes. *Journal of the Royal Statistical Society: Series B*,
50(2), 297--312. <https://doi.org/10.1111/j.2517-6161.1988.tb01729.x>

Guinness, J. (2018). Permutation and grouping methods for sharpening Gaussian
process approximations. *Technometrics*, 60(4), 415--429.
<https://doi.org/10.1080/00401706.2018.1437476>
