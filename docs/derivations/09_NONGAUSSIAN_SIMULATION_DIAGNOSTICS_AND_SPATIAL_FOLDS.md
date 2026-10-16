# Non-Gaussian simulation, residuals and spatial folds

## Purpose and prerequisites

This chapter derives fitted response simulation and residual diagnostics for
binomial-logit and Poisson-log spatial models, then specifies the coordinate
clustering used to construct spatial cross-validation folds. The reader should
know conditional probability, binomial and Poisson means and variances, and
squared Euclidean distance.

## Fitted-model response simulation

For fitted parameters \(\widehat{\boldsymbol\beta}\) and
\(\widehat{\boldsymbol\psi}\), let

\[
\boldsymbol s^{(b)}\sim
N(\boldsymbol0,\boldsymbol\Sigma_{\widehat{\boldsymbol\psi}}),
\qquad b=1,\ldots,B.                                             \tag{S-01}
\]

If \(\boldsymbol R^{\mathsf T}\boldsymbol R=
\boldsymbol\Sigma_{\widehat{\boldsymbol\psi}}\), SPmodels generates

\[
\boldsymbol z^{(b)}\sim N(\boldsymbol0,\boldsymbol I),
\qquad
\boldsymbol s^{(b)}=\boldsymbol R^{\mathsf T}\boldsymbol z^{(b)}.\tag{S-02}
\]

The simulated linear predictor is

\[
\eta_i^{(b)}=\boldsymbol x_i^{\mathsf T}
\widehat{\boldsymbol\beta}+o_i+s_i^{(b)}.                        \tag{S-03}
\]

For a binomial model with fitted trial total \(m_i\), simulation continues as

\[
p_i^{(b)}=\operatorname{logit}^{-1}(\eta_i^{(b)}),
\qquad
Y_i^{(b)}\mid s_i^{(b)}\sim
\operatorname{Binomial}(m_i,p_i^{(b)}).                          \tag{S-04}
\]

The returned value is a success count, including when the original response
was entered as a proportion. For Poisson-log models,

\[
\lambda_i^{(b)}=\exp(\eta_i^{(b)}),
\qquad
Y_i^{(b)}\mid s_i^{(b)}\sim
\operatorname{Poisson}(\lambda_i^{(b)}).                        \tag{S-05}
\]

Equations (S-01)--(S-05) are unconditional on the observed responses and
conditional on fitted parameters. They are appropriate for fitted-model checks
and simulation studies. They do not draw parameter uncertainty. A supplied
seed is scoped to the function and the caller's random-number state is restored.

## Conditional residuals at the latent mode

Let \(\widehat\mu_i\) be the fitted conditional response mean at the final
latent mode. The raw response residual is

\[
r_i=y_i-\widehat\mu_i.                                           \tag{S-06}
\]

The conditional response variance is

\[
V_i=
\begin{cases}
m_i\widehat p_i(1-\widehat p_i), & \text{binomial},\\
\widehat\mu_i, & \text{Poisson},
\end{cases}
\qquad \widehat p_i=\widehat\mu_i/m_i.                           \tag{S-07}
\]

When \(V_i>0\), the Pearson residual is

\[
r_i^{P}=\frac{y_i-\widehat\mu_i}{\sqrt{V_i}}.                    \tag{S-08}
\]

Zero conditional variance has no finite Pearson standardization and is
reported as missing rather than replaced by an arbitrary number.

For Poisson data, the deviance contribution is

\[
d_i^2=2\left[y_i\log\left(\frac{y_i}{\widehat\mu_i}\right)
-(y_i-\widehat\mu_i)\right],                                    \tag{S-09}
\]

with \(0\log0=0\). For binomial counts,

\[
d_i^2=2\left[
y_i\log\left(\frac{y_i}{\widehat\mu_i}\right)
+(m_i-y_i)\log\left(
\frac{m_i-y_i}{m_i-\widehat\mu_i}\right)
\right],                                                        \tag{S-10}
\]

with the analogous zero-count conventions. The signed deviance residual is

\[
r_i^D=\operatorname{sign}(r_i)\sqrt{d_i^2}.                      \tag{S-11}
\]

For a discrete fitted conditional CDF \(F_i\), define

\[
a_i=F_i(y_i-1),\qquad b_i=F_i(y_i).                            \tag{S-11a}
\]

Draw a caller-seeded value and transform it:

\[
U_i\sim\operatorname{Uniform}(a_i,b_i),\qquad
r_i^Q=\Phi^{-1}(U_i).                                          \tag{S-11b}
\]

Under the correct discrete conditional distribution, the probability
intervals partition \([0,1]\). Mixing over them makes \(U_i\) marginally
uniform, hence

\[
r_i^Q\sim N(0,1).                                              \tag{S-11c}
\]

Probabilities are clipped only at machine epsilon before \(\Phi^{-1}\),
preventing infinite computer values. These residuals remain conditional on the
fitted latent mode and parameters; spatial dependence and parameter fitting are
not removed. They are distinct from Gaussian Cholesky-decorrelated residuals.

For any chosen non-Gaussian residual \(r_i\) above, the empirical residual
variogram in distance bin \(P_k\) is

\[
\widehat\gamma_k^{(r)}=
\frac{1}{|P_k|}\sum_{(i,j)\in P_k}\frac12(r_i-r_j)^2.          \tag{S-11d}
\]

SPmodels supports Pearson, signed deviance and caller-seeded randomized-
quantile residuals. It returns pair counts and the empirical curve but no
fitted response-scale curve: the latent link covariance is not silently
treated as a semivariogram for discrete residuals.

For a permutation envelope, replicate \(b\) applies a uniformly sampled
permutation \(\pi_b\) to the selected residual vector and recomputes (S-11d):

\[
\widehat\gamma_{kb}^{\mathrm{perm}}=
\frac{1}{|P_k|}\sum_{(i,j)\in P_k}
\frac12(r_{\pi_b(i)}-r_{\pi_b(j)})^2.                         \tag{S-11e}
\]

Pointwise type-8 quantiles, the standardized global curve statistic and the
finite Monte Carlo p-value then follow equations (D-07)--(D-10). The null is
exchangeability of the chosen residual values across fixed locations. This is
a conditional diagnostic, not a proof of model adequacy; heterogeneous or
poorly calibrated residual distributions weaken that interpretation.

For a fitted-model parametric envelope, replicate $b$ instead follows the
complete discrete-data workflow

\[
Y_b^*\sim p(\,\cdot\mid\widehat\theta),\qquad
\widehat\theta_b^*=\operatorname{fit}(Y_b^*),\qquad
\widehat\gamma_{kb}^{\mathrm{par}}=
\frac{1}{|P_k|}\sum_{(i,j)\in P_k}
\frac12(r_{ib}^*-r_{jb}^*)^2.                                \tag{S-11f}
\]

Here `fit` repeats the model's declared Laplace or Monte Carlo ML route, and
$r_{ib}^*$ is the chosen Pearson, deviance or randomized-quantile residual
from that refit. Thus the envelope includes response simulation and refitting
variation rather than reusing the original latent mode. Every optimizer
convergence code and captured warning is retained. A non-zero code is reported
to the caller and must be examined; SPmodels does not relabel it as success.
The pointwise and global summaries again use (D-07)--(D-10).

## Geometry-driven spatial folds

Random folds can place nearby or repeated locations on both sides of a
validation split. To prevent exact-location leakage, define coordinate
equivalence classes

\[
i\sim j \quad\Longleftrightarrow\quad
\boldsymbol c_i=\boldsymbol c_j.                                 \tag{S-12}
\]

SPmodels clusters the \(U\) unique coordinate vectors, not all \(n\) rows.
For requested fold count \(K\le U\), k-means selects clusters
\(C_1,\ldots,C_K\) and centroids \(\boldsymbol a_1,\ldots,
\boldsymbol a_K\) to reduce

\[
\sum_{k=1}^K\sum_{u\in C_k}
\|\boldsymbol c_u-\boldsymbol a_k\|^2.                           \tag{S-13}
\]

All original rows in one equivalence class inherit the same cluster. Cluster
numbers from k-means are arbitrary, so SPmodels canonicalizes labels by sorting
centroids from largest to smallest second coordinate, then smallest to largest
first coordinate. This makes labels interpretable as top-to-bottom and
left-to-right; the partition itself can still depend on initialization.

### Regular spatial blocks

For declared positive axis widths \(w_1,w_2\) and grid origin
\(\boldsymbol o=(o_1,o_2)\), row \(i\) belongs to occupied block

\[
B_i=\left(\left\lfloor\frac{c_{i1}-o_1}{w_1}\right\rfloor,
          \left\lfloor\frac{c_{i2}-o_2}{w_2}\right\rfloor\right). \tag{S-13c}
\]

Every row with the same \(B_i\), including repeated coordinates, receives one
fold label. Occupied blocks are ordered top-to-bottom then left-to-right and
assigned round-robin, or the occupied-block order is locally seeded and
permuted before the same balanced allocation:

\[
f_i=a(B_i),\qquad a(B)\in\{1,\ldots,K\}.                      \tag{S-13d}
\]

The operation balances numbers of occupied blocks, not necessarily numbers of
observations or outcomes. It retains block indices, row counts, widths,
origin, units and allocation rule so the split is auditable. Blocks control
which supports are withheld together; they do not impose a distance gap.

### Environmental-space folds

Let \(e_{ij}\) be environmental variable \(j\) at unique location \(i\).
With the default standardization, clustering uses

\[
z_{ij}=\frac{e_{ij}-\bar e_j}{s_j},                            \tag{S-13e}
\]

where a zero-variance variable is rejected. K-means then chooses environmental
clusters \(G_1,\ldots,G_K\) and centres \(\boldsymbol m_k\) to reduce

\[
\sum_{k=1}^K\sum_{i\in G_k}
\|\boldsymbol z_i-\boldsymbol m_k\|^2.                        \tag{S-13f}
\]

Rows at an identical coordinate must have identical declared environmental
values and inherit one cluster, preventing exact-location leakage. Centres are
canonically ordered lexicographically, and the seed is local. These folds test
transfer across environmental regimes; they need not be geographically
compact and do not impose a distance buffer.

### Group-constrained folds

Let $g(i)\in\{1,\ldots,G\}$ identify the declared clinic, survey, school or
other grouping unit for row $i$, and let

\[
n_g=|\{i:g(i)=g\}|.                                           \tag{S-13g}
\]

Every group is indivisible. In deterministic mode, groups are considered from
largest to smallest $n_g$, with first appearance breaking ties. In random
mode their order is locally seeded and permuted. If $L_k$ is the number of
rows already assigned to fold $k$, the next group goes to

\[
k^*=\operatorname*{arg\,min}_{1\le k\le K}L_k,\qquad
L_{k^*}\leftarrow L_{k^*}+n_g.                               \tag{S-13h}
\]

The first fold breaks load ties. This greedy rule targets row balance without
ever splitting a group; exact equality is not promised when group sizes differ.
Returned metadata records every group label, row count and allocation.

### Buffered training sets

For test fold \(\mathcal H_k\) and positive buffer distance \(b\), define the
fold-specific training set

\[
\mathcal T_k(b)=\left\{i:f_i\ne k,\thinspace
\min_{j\in\mathcal H_k}d_{ij}>b\right\}.                       \tag{S-13a}
\]

Rows with \(f_i\ne k\) that fall inside the buffer are neither used for
training nor scored as test observations in fold \(k\). The held-out set stays
\(\mathcal H_k\). The exclusion set is therefore

\[
\mathcal E_k(b)=\{i:f_i\ne k\}\setminus\mathcal T_k(b).       \tag{S-13b}
\]

The buffer uses the fitted coordinate units and is applied from every test
location. This reduces short-range leakage but can greatly reduce the training
sample. A fold with no remaining training rows fails explicitly.

Multiple random starts reduce sensitivity to poor local k-means solutions. A
supplied seed is local. These folds are spatially compact but are not promised
to have equal observation counts, equal response prevalence or a minimum
between-fold separation. Scientific applications should inspect the returned
centres and fold sizes.

## Repeated simulation and parameter recovery

The same fitted-model generator supports a repeated recovery experiment. Let
\(\widehat{\boldsymbol\theta}\) be the supplied fitted parameter vector. For
replicate \(b\), generate

\[
\boldsymbol Y^{(b)}\sim
p(\boldsymbol y\mid\widehat{\boldsymbol\theta})                 \tag{S-14}
\]

using the family-specific hierarchy above, then refit the same declared model
to obtain \(\widehat{\boldsymbol\theta}^{(b)}\). If the supplied fit used MCML,
both the preliminary Laplace stage and MCML correction are repeated.

The recovery vector contains every parameter estimated by the supplied route:
mean coefficients; spatial and nugget covariance parameters; grouped
random-intercept variances; temporal scale; and intervention impact/decay when
their corresponding estimation flags are true. Fixed quantities are retained
by each refit but are not mislabelled as recovered estimates. Qualified
`random::` names prevent a group variance from colliding with a mean or spatial
parameter name. For parameter \(j\), define replicate error

\[
e_{bj}=\widehat\theta_j^{(b)}-\widehat\theta_j.                   \tag{S-15}
\]

Across \(B_j\) successful refits, estimated bias and RMSE are

\[
\widehat{\operatorname{Bias}}_j=
\frac1{B_j}\sum_b e_{bj},
\qquad
\widehat{\operatorname{RMSE}}_j=
\left(\frac1{B_j}\sum_b e_{bj}^2\right)^{1/2}.                  \tag{S-16}
\]

The empirical standard deviation is

\[
s_j=\left\{\frac1{B_j-1}\sum_b
(\widehat\theta_j^{(b)}-\overline{\widehat\theta}_j)^2
\right\}^{1/2}.                                                  \tag{S-17}
\]

If replicate \(b\) reports usable standard error
\(\widehat{\operatorname{se}}_{bj}\), SPmodels also calculates their mean and
the Wald coverage indicator

\[
I_{bj}=\mathbb 1\left[
\widehat\theta_j^{(b)}-z_{(1+\gamma)/2}
\widehat{\operatorname{se}}_{bj}
\le \widehat\theta_j \le
\widehat\theta_j^{(b)}+z_{(1+\gamma)/2}
\widehat{\operatorname{se}}_{bj}\right],                        \tag{S-18}
\]

and coverage

\[
\widehat C_j=\frac1{B_j^{\mathrm{se}}}
\sum_b I_{bj}.                                                   \tag{S-19}
\]

Failed refits are never converted to successful values. They are recorded with
their messages, and every parameter summary reports its successful count. The
caller may instead request immediate failure.

The response matrix is generated before any optional forked dispatch. For
replicate $b$, MCML receives the deterministic seed

\[
s_b=(s+b)\bmod M,                                             \tag{S-19a}
\]

where $s$ is the declared base seed and $M$ is R's maximum integer. Therefore
worker scheduling cannot change simulated responses or stochastic refit
streams, and sequential and parallel estimates agree exactly. Parallel MCML
requires an explicit base seed; Windows currently uses one worker. The caller's
random-number state is restored before workers are launched.

This experiment treats the fitted values as truth. It measures recovery near
one fitted scenario, conditional on its design, coordinates and offsets. It is
not a substitute for a predeclared simulation grid spanning scientifically
important parameter values.

## Multiple fitted scenarios

Let \(q=1,\ldots,Q\) index named fitted models. Each model defines its own
truth \(\widehat\theta_q\), design, covariance, family and replicate count
\(B_q\). SPmodels runs the preceding experiment separately:

\[
\mathcal R_q=\operatorname{Recovery}
(\widehat\theta_q,B_q,s_q),\qquad s_q=s_1+q-1,                 \tag{S-20}
\]

where integer seeds wrap at R's maximum integer. There is no pooling across
scenarios: rows are stacked only after attaching the scenario identifier.
Every `spmodels_simulation_recovery` object, simulated response, failure and
estimate remains stored. This orchestration makes a predeclared simulation
design easier to audit, but does not make an arbitrary list of scenarios
scientifically comprehensive. Scenario seeds are resolved before dispatch.
Consequently, the sequential route and optional non-Windows forked scenario
workers produce identical per-scenario experiments; parallel scheduling does
not choose the random streams.

## Precision-directed replicate escalation

For attempt \(a\), let \(B_a\) be its declared replicate count and let
\(B_{aj}\) estimates of parameter \(j\) be finite. The Monte Carlo standard
error of estimated mean bias is estimated by

\[
\widehat{\operatorname{MCSE}}_{aj}=
\frac{s_{aj}}{\sqrt{B_{aj}}},                                  \tag{S-21}
\]

where \(s_{aj}\) is the empirical standard deviation in (S-17). An attempt
passes only if every recovered parameter has at least the declared minimum
successful count and

\[
\widehat{\operatorname{MCSE}}_{aj}\le\epsilon
\quad\text{for every }j.                                      \tag{S-22}
\]

Otherwise the requested count grows geometrically and is capped exactly at
\(B_{\max}\). Every attempt is a complete independently seeded recovery
experiment; its count, seed, minimum successes, worst finite MCSE and status
are retained. A failed capped attempt remains a failure. The threshold is
never relaxed to obtain a pass. Each complete attempt may use the replicate
workers in (S-19a); tests require its estimates and stopping history to equal
the sequential route exactly.

## Cartesian parameter scenarios

Let parameter $j$ have a finite declared set of truths

\[
\Theta_j=\{\theta_{j1},\ldots,\theta_{jm_j}\}.
                                                               \tag{S-23}
\]

`sp_simulation_grid()` constructs the Cartesian product

\[
\mathcal G=\Theta_1\times\cdots\times\Theta_q,
\qquad |\mathcal G|=\prod_{j=1}^{q}m_j.                       \tag{S-24}
\]

R's first grid column varies fastest, and SPmodels records the resulting table
as an attribute in exactly the returned scenario order. Only parameters that
the unchanged refit specification can recover are accepted. This prevents a
misleading experiment in which data are generated under one fixed covariance
but every replicate is refitted under a different fixed covariance.

For each row of \(\mathcal G\), SPmodels replaces the declared coefficient,
covariance, random-variance, temporal or intervention values and rebuilds the
generative covariance factor. For Gaussian models this is the Cholesky factor
of \(V\) in (S-01). For binomial and Poisson models it is the Cholesky factor
of \(\Sigma\), followed by recomputation of the latent mode so the scenario is
internally coherent. The original formula, locations, response family, trial
totals and fitting controls are retained. These objects define simulation
truths; they are not presented as new fits to the observed response.

## Code and validation map

| Equations | Implementation | Independent checks |
|---|---|---|
| S-01--S-05 | `simulate.spmodels_fit()` | Exact reconstruction from seeded normal and response-family RNG draws; response support |
| S-06--S-11c | `sp_diagnostics()` | Direct mean, variance, Pearson and deviance calculations; exact seeded CDF-interval uniform and normal transform |
| S-11d | `sp_variogram()` | Direct binning of seeded randomized-quantile pair semivariances, pair counts, RNG isolation and rendered empirical-only display |
| S-11e--S-11f and D-07--D-10 | `sp_variogram_envelope()` | Exact reconstruction of every seeded non-Gaussian permutation curve and the first full simulate-refit parametric curve, deterministic repeats, retained refit convergence/warnings, finite p-value bounds and empirical-only display |
| S-12--S-13 | `sp_spatial_folds()` | Duplicate coordinates remain together, local seeded reproducibility, canonical labels and cross-validation use |
| S-13c--S-13d | `sp_block_folds()`, `sp_cross_validate()` | Direct block-index reconstruction, whole-block retention, deterministic/random allocation, RNG isolation and cross-validation use |
| S-13e--S-13f | `sp_environmental_folds()`, `sp_cross_validate()` | Formula/name equivalence, standardized centres, local RNG isolation and environmental-fold cross-validation use |
| S-13g--S-13h | `sp_group_folds()`, `sp_cross_validate()` | Direct greedy allocation, whole-group retention, local seeded reproducibility, retained metadata and complete cross-validation use |
| S-13a--S-13b | `sp_buffered_folds()`, `sp_cross_validate()` | Every near-test row excluded; one fold reproduced by a direct manually buffered refit; oversized buffer failure |
| S-14--S-19a | `sp_simulation_recovery()` | First replicate matched to an independent refit; direct bias/RMSE reconstruction; binomial count/trial support; exact sequential/parallel result equality and RNG isolation |
| S-20 | `sp_simulation_scenarios()` | Each sequential or parallel scenario exactly matches its direct seeded recovery call; input order, unequal replicate counts, RNG isolation and combined labels |
| S-21--S-22 | `sp_simulation_precision()` | Direct MCSE reconstruction, first-attempt pass, geometric escalation, exact cap, RNG isolation and exact sequential/parallel attempt equality |
| S-23--S-24 | `sp_simulation_grid()` | Exact Cartesian order, parameter substitution, direct Gaussian covariance reconstruction, non-Gaussian Cholesky reconstruction, response support and rejection of fixed/non-recoverable parameters |

## Limitations

Simulation is dense and uses the fitted observation locations and trial totals.
It supports declared new parameter values over Cartesian grids but does not
yet simulate new spatial supports.
Non-Gaussian diagnostics include empirical residual variograms and permutation
or full simulate-refit parametric envelopes/global curve tests. Spatial fold test groups can use coordinate
k-means, regular axis-aligned blocks or environmental-space clusters;
group-constrained allocation is also available, while arbitrary user-supplied
fold identifiers remain accepted. Further domain-specific constructors should
be added only for a declared use case. Buffered training exclusion is implemented
around geometry-driven k-means test groups.
Repeated recovery supports named response variables, two-column binomial
responses, Cartesian parameter grids and multiple named fitted scenarios with
unequal replicate counts. Precision-directed replicate
escalation is available for one fitted scenario. Named scenarios and individual
replicates can run in parallel on non-Windows systems with preassigned seeds.
