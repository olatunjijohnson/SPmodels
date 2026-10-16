# Separable space-time geostatistical models

## Purpose, data and assumptions

Repeated surveys can be close in space but separated in time. SPmodels
represents this residual dependence by multiplying Matérn spatial correlation
by exponential temporal correlation. Each point observation needs one
two-dimensional projected location and one finite numeric time.

The model assumes stationarity and isotropy in space, stationarity in time,
separability, common covariance parameters across the study, and independent
observation-specific nugget variation. It does not claim non-separable
interaction, directional anisotropy, seasonality or non-stationarity.

## From coordinates and times to correlation

For locations \(s_i,s_j\in\mathbb R^2\) and times \(t_i,t_j\), define

\[
d_{ij}=\lVert s_i-s_j\rVert,\qquad u_{ij}=|t_i-t_j|.           \tag{T-01}
\]

Let \(R_M(d;\phi,\nu)\) be the Matérn correlation derived in Chapter 1.
With temporal scale \(\lambda>0\),

\[
R_T(u;\lambda)=\exp(-u/\lambda),                              \tag{T-02}
\]

and the separable correlation is

\[
R_{ST}\{(s_i,t_i),(s_j,t_j)\}
=R_M(d_{ij};\phi,\nu)R_T(u_{ij};\lambda).                     \tag{T-03}
\]

For finite correlation matrices \(R_S\) and \(R_T\),

\[
R_{ST}=R_S\circ R_T.                                          \tag{T-04}
\]

The Schur product theorem makes this elementwise product positive
semidefinite. Its diagonal is one, so it is a valid correlation matrix. At
zero time lag it reduces to Matérn correlation; at zero spatial distance it
reduces to exponential temporal correlation.

## Gaussian model and likelihood

For a Gaussian response,

\[
Y=X\beta+o+Zb+S+\epsilon,                                    \tag{T-05}
\]

where \(Zb\) contains any grouped random intercepts. With spatial variance
\(\sigma_s^2\), nugget variance \(\sigma_e^2\), and grouped-effect
covariance \(ZGZ^{\mathsf T}\),

\[
V=\sigma_s^2R_{ST}+ZGZ^{\mathsf T}+\sigma_e^2I.               \tag{T-06}
\]

For fixed covariance parameters, generalized least squares gives

\[
\widehat\beta=(X^{\mathsf T}V^{-1}X)^{-1}
X^{\mathsf T}V^{-1}(Y-o).                                    \tag{T-07}
\]

Writing \(e=Y-o-X\widehat\beta\), the profiled ordinary log-likelihood is

\[
\ell=-\tfrac12\{n\log(2\pi)+\log|V|+e^{\mathsf T}V^{-1}e\}. \tag{T-08}
\]

SPmodels uses Cholesky solves. If temporal scale is estimated, it optimizes
\(\log\lambda\) jointly with log spatial scale and log variances. The delta
method transforms usable observed-Hessian uncertainty back to its ordinary
scale:

\[
\operatorname{Var}(\widehat\lambda)\approx\widehat\lambda^2
\operatorname{Var}\{\log(\widehat\lambda)\}.                 \tag{T-09}
\]

Indefinite curvature is reported as missing uncertainty, not repaired by
weakening a numerical criterion.

## Prediction at new locations and times

For new design \(X_*\), times \(t_*\), and locations \(s_*\), the
training-to-new covariance has elements

\[
(C_*)_{ij}=\sigma_s^2R_M(\lVert s_i-s_{*j}\rVert;\phi,\nu)
\exp(-|t_i-t_{*j}|/\lambda),                                  \tag{T-10}
\]

plus any known-level grouped-effect covariance. The conditional mean is

\[
\widehat Y_*=X_*\widehat\beta+o_*+C_*^{\mathsf T}V^{-1}e.     \tag{T-11}
\]

Let \(K_{**}\) contain all new-to-new space-time covariances and define

\[
Q=X_*-C_*^{\mathsf T}V^{-1}X.                                \tag{T-12}
\]

The joint latent predictive covariance is

\[
\Sigma_*=K_{**}-C_*^{\mathsf T}V^{-1}C_*
+Q(X^{\mathsf T}V^{-1}X)^{-1}Q^{\mathsf T}.                  \tag{T-13}
\]

Future-observation prediction adds \(\sigma_e^2I\). Point standard errors
use the diagonal of (T-13); joint draws use its full matrix. Optional Gaussian
parameter draws include \(\log\lambda\) when its Hessian covariance is usable.

## Binomial and Poisson models

For binomial-logit and Poisson-log fits, the latent Gaussian prior is

\[
\Sigma=\sigma_s^2R_{ST}+\sigma_e^2I.                          \tag{T-14}
\]

The Chapter 5 Laplace derivation is unchanged after substituting (T-14). For
each trial value of \(\log\lambda\), SPmodels rebuilds (T-14), finds the latent
mode, and includes both the response-plus-prior value and determinant
correction in the Laplace objective. Thus temporal scale is estimated jointly
with the regression and spatial covariance parameters rather than from a
separate empirical correlation step. Equation (T-09) gives its natural-scale
local curvature uncertainty when the full observed Hessian is usable.

MCML keeps the same transformed parameter in its deterministic fixed-draw
objective. At every trial value it evaluates the exact Gaussian prior density
with the corresponding space-time covariance before dividing by the fixed
normalized proposal density. The final latent state, prediction, simulation,
validation refits and parameter-recovery summaries all use the resulting
\(\widehat\lambda\). This extends the likelihood calculation; it does not
remove possible weak identification caused by sparse or confounded survey
times.

## Simulation, refitting and variograms

Gaussian simulation at fitted sites uses

\[
Y^{(b)}=X\widehat\beta+o+L^{\mathsf T}z^{(b)},\quad
z^{(b)}\sim N(0,I),\quad L^{\mathsf T}L=V.                    \tag{T-15}
\]

Non-Gaussian simulation transforms a draw from (T-14) through its response
model. Validation and simulation-refitting pass the temporal specification to
every refit and read held-out times from prediction data.

For Gaussian pair \((i,j)\), the fitted process-plus-nugget semivariance is

\[
\gamma_{ij}=\sigma_e^2+\sigma_s^2
\{1-R_M(d_{ij};\phi,\nu)\exp(-u_{ij}/\lambda)\}.              \tag{T-16}
\]

A distance bin can contain different time lags, so SPmodels averages the
pair-specific fitted values rather than using only a distance midpoint:

\[
\widehat\gamma_k=|B_k|^{-1}\sum_{(i,j)\in B_k}\gamma_{ij}.    \tag{T-17}
\]

Grouped random-intercept pair contributions are added as in Chapter 13.
Parametric variogram envelopes simulate from the fitted space-time covariance.

## Units, interpretation and identifiable limits

Spatial scale uses coordinate units; temporal scale uses the declared
`time_units`. Multiplying all times and \(\lambda\) by the same positive
constant leaves correlation unchanged. A large temporal scale means slow
decay, not a hard duration. If all times coincide, temporal scale is not
identifiable, so an estimation request is rejected explicitly. Sparse time
coverage or distance-time confounding can likewise
produce weak curvature, so estimates require study-design judgment.

## Implementation and independent checks

| Equations | Implementation | Independent checks |
|---|---|---|
| T-01--T-04 | `spacetime_cor()`, `sp_temporal()` | Exponential and limiting cases, unit rescaling, symmetry and eigenvalue PSD |
| T-05--T-09 | `sp_fit()`, covariance/state helpers | Direct dense GLS, determinant likelihood, positive scale and Hessian schema |
| T-10--T-13 | `predict()`, `sp_predict_draws()` | Direct cross/new covariance, point variance and joint covariance at new times |
| T-14 | non-Gaussian fitting, prediction and MCML | Direct fitted covariance, positive log-scale estimate, fitted degree of freedom, MCML propagation and finite response/joint prediction |
| T-15 | `simulate()` | Exact seeded Cholesky reconstruction |
| T-16--T-17 | `sp_variogram()`, envelope code | Direct pairwise time-aware bin averages |
| refitting/reporting | validation and table methods | Specification retained; scale status and units exposed |

Primitive checks are in `tests/test-spacetime.R`; fitted-workflow checks are in
`tests/test-spacetime-fit.R`. These establish the stated equations for small,
well-conditioned fixtures. They do not establish external reference agreement,
broad scientific validity or large-data performance.
