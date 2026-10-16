# Grouped random intercepts

## Purpose, assumptions and current boundary

Measurements can share both spatial dependence and membership of a clinic,
survey, laboratory or other grouping unit. A Gaussian random intercept assigns
one common zero-mean deviation to every observation in the same level. SPmodels
supports one or more additive grouping variables in Gaussian, binomial-logit
and Poisson-log spatial models through `random = ~ group1 + group2`.

The interface derived in this chapter assumes independent Gaussian intercepts
across levels and grouping variables. Fixed grouped slopes and correlated
intercept-slope effects use the separate explicit interface derived in Chapter
17; neither interface claims arbitrary mixed-model formula syntax or
non-spatial grouped fits. Non-Gaussian grouped models support MCML correction. Prediction at a
new level is available only through an explicit marginal policy; it is never
silently treated as a known level.

## Model and indicator matrices

For observations \(i=1,\ldots,n\), let \(x_i\) contain fixed-effect
covariates and let \(S_i\) be the Matérn spatial effect. Suppose grouping
variable \(r\) has \(L_r\) observed levels. Its indicator matrix \(Z_r\) has
entries

\[
(Z_r)_{i\ell}=\mathbb 1\{\text{observation }i
                    \text{ belongs to level }\ell\}.           \tag{Q-01}
\]

Every row has exactly one unit entry for each grouping variable. The model is

\[
Y=X\beta+o+S+\sum_{r=1}^R Z_r b_r+\varepsilon,                 \tag{Q-02}
\]

with mutually independent components

\[
S\sim N(0,\sigma_s^2R_\theta),\qquad
b_r\sim N(0,\omega_r^2I_{L_r}),\qquad
\varepsilon\sim N(0,\tau^2I_n).                               \tag{Q-03}
\]

Here \(\omega_r^2\) is the variance of the random intercept for grouping
variable \(r\). It is not the variance of a group mean and is not a residual
variance. Integrating out the random intercepts gives

\[
V=\sigma_s^2R_\theta+
  \sum_{r=1}^R\omega_r^2Z_rZ_r^{\mathsf T}+\tau^2I_n.          \tag{Q-04}
\]

Thus two observations in the same level of grouping variable \(r\) receive
additional covariance \(\omega_r^2\); observations in different levels do
not. Contributions from multiple grouping variables add.

## Maximum likelihood

For fixed variance parameters, the generalized least-squares estimate is

\[
\widehat\beta=(X^{\mathsf T}V^{-1}X)^{-1}
                X^{\mathsf T}V^{-1}(y-o).                      \tag{Q-05}
\]

Substitution into the Gaussian log likelihood gives

\[
\ell=-\frac12\left[
n\log(2\pi)+\log|V|+
(y-o-X\widehat\beta)^{\mathsf T}V^{-1}
(y-o-X\widehat\beta)
\right].                                                       \tag{Q-06}
\]

SPmodels evaluates (Q-05)--(Q-06) using Cholesky solves, never an explicit
matrix inverse. When estimation is requested, every strictly positive
covariance parameter, including each \(\omega_r^2\), is optimized on its log
scale. With fixed covariance controls, one positive random variance must be
provided for each grouping variable. A scalar is repeated; a named vector must
match the grouping-variable names exactly.

The observed optimizer Hessian is inverted on the log scale when positive
definite. The delta method maps the random-variance block back to the ordinary
variance scale. An unusable Hessian produces missing standard errors, not
invented certainty.

## Prediction at known grouping levels

For new prediction rows with fixed-effect matrix \(X_*\), coordinates and
known-level indicator matrices \(Z_{r*}\), define the training-to-new
covariance

\[
C=\sigma_s^2R_{\theta,\,*}+
  \sum_{r=1}^R\omega_r^2Z_rZ_{r*}^{\mathsf T},                 \tag{Q-07}
\]

and the new-to-new latent covariance

\[
V_*^{(0)}=\sigma_s^2R_{\theta,**}+
  \sum_{r=1}^R\omega_r^2Z_{r*}Z_{r*}^{\mathsf T}.              \tag{Q-08}
\]

The universal-kriging mean remains

\[
\widehat y_*=X_*\widehat\beta+o_*+
C^{\mathsf T}V^{-1}(y-o-X\widehat\beta),                      \tag{Q-09}
\]

and the joint latent prediction covariance is

\[
\begin{aligned}
\Sigma_*={}&V_*^{(0)}-C^{\mathsf T}V^{-1}C\\
&+(X_*-C^{\mathsf T}V^{-1}X)
 (X^{\mathsf T}V^{-1}X)^{-1}
 (X_*-C^{\mathsf T}V^{-1}X)^{\mathsf T}.
\end{aligned}                                                   \tag{Q-10}
\]

Observation prediction adds \(\tau^2I\). Joint draws, pointwise standard
errors and conditional means all use the same covariance. By default, new
group levels are rejected because their conditional relationship to fitted
levels differs from known-level prediction.

## Explicit marginal prediction at new grouping levels

Under `unseen = "marginal"`, a previously unseen level receives a fresh
zero-mean random intercept from the fitted distribution. For grouping variable
\(r\), split the new-row indicators into \(Z_{r*}^{K}\), whose columns refer
to fitted levels, and \(U_{r*}\), whose columns refer only to distinct unseen
levels in the prediction data. Then

\[
C=\sigma_s^2R_{\theta,*}+
  \sum_{r=1}^R\omega_r^2Z_r(Z_{r*}^{K})^{\mathsf T},          \tag{Q-10a}
\]

because a new random intercept is independent of every training intercept,
whereas the new-to-new covariance becomes

\[
V_*^{(0)}=\sigma_s^2R_{\theta,**}+
\sum_{r=1}^R\omega_r^2\left\lbrace
Z_{r*}^{K}(Z_{r*}^{K})^{\mathsf T}+U_{r*}U_{r*}^{\mathsf T}
\right\rbrace.                                               \tag{Q-10b}
\]

Consequently, a new level does not shift the conditional mean: its marginal
mean is zero. Its variance \(\omega_r^2\) is retained, however, and two
prediction rows carrying the same new level have covariance \(\omega_r^2\).
Rows carrying distinct new levels have zero covariance from that grouping
variable. Substituting (Q-10a)--(Q-10b) into (Q-10) or (Q-15) gives the Gaussian
or non-Gaussian predictive covariance respectively:

\[
\operatorname{Cov}(b_{r*}^{\mathrm{new}})=
\omega_r^2U_{r*}U_{r*}^{\mathsf T}.                          \tag{Q-10c}
\]

This is a population-distribution prediction for new levels, not an estimate
of any particular unseen level's realized intercept. Parameter uncertainty,
when requested, rebuilds (Q-10a)--(Q-10c) with the variance from each parameter
draw.

## Simulation and variograms

Unconditional fitted simulation draws

\[
Y^{(b)}=X\widehat\beta+o+L^{\mathsf T}z^{(b)},\qquad
L^{\mathsf T}L=V,\quad z^{(b)}\sim N(0,I).                     \tag{Q-11}
\]

For two distinct observations \(i,j\), grouping variable \(r\) contributes
to the model semivariogram

\[
\gamma_{r,ij}=\omega_r^2
  \mathbb 1\{g_{ir}\ne g_{jr}\}.                              \tag{Q-12}
\]

Because this term depends on group membership rather than distance alone,
SPmodels averages (Q-12) over the actual pairs in each distance bin and adds it
to the Matérn-plus-nugget fitted curve. This is more faithful than pretending a
group effect is an isotropic function of distance.

## Binomial and Poisson latent models

For a supported non-Gaussian response, grouped effects enter the link-scale
latent predictor:

\[
\eta=X\beta+o+S+\sum_{r=1}^R Z_rb_r.                           \tag{Q-13}
\]

After integrating the independent Gaussian intercepts, the latent covariance
is

\[
\Sigma=\sigma_s^2R_\theta+
\sum_{r=1}^R\omega_r^2Z_rZ_r^{\mathsf T}+\tau^2I.             \tag{Q-14}
\]

The Laplace mode, curvature and approximate likelihood in Chapter 5 follow
after replacing its covariance by (Q-14). Fixed random variances are supplied
in `sp_control()`; when covariance estimation is enabled their logarithms are
optimized jointly, and usable Hessian curvature is delta-transformed.

For known new levels, cross and new-new covariance use the grouped terms in
(Q-07)--(Q-08). If \(C\) is the resulting cross-covariance and \(H^{-1}\) the
Laplace mode covariance, the approximate link covariance is

\[
\Sigma_* = V_*^{(0)}-C^{\mathsf T}\Sigma^{-1}C
+(\Sigma^{-1}C)^{\mathsf T}H^{-1}(\Sigma^{-1}C).               \tag{Q-15}
\]

Unconditional simulation draws the complete latent vector then applies its
response distribution:

\[
S_Z^{(b)}\sim N(0,\Sigma),\qquad
Y_i^{(b)}\mid S_{Z,i}^{(b)}\sim
F_i\{x_i^{\mathsf T}\widehat\beta+o_i+S_{Z,i}^{(b)}\}.         \tag{Q-16}
\]

MCML uses the same covariance (Q-14), with estimated random variances included
on the log scale in its common-random-number objective. Predictive parameter
draws rebuild the complete latent covariance for every draw. Laplace fitting,
MCML, prediction, simulation, diagnostics, validation and reporting therefore
share one grouping contract.

## Code and validation map

| Equations | Implementation | Independent checks |
|---|---|---|
| Q-01--Q-04 | `.sp_random_training()`, `.sp_random_covariance()`, `.sp_covariance()` | Two grouping variables; direct indicator matrices and full covariance equality |
| Q-05--Q-06 | `sp_fit()` | Direct GLS coefficients and Gaussian log likelihood |
| Q-07--Q-10 | `predict()` and `sp_predict_draws()` | Direct known-level cross-covariance, mean, point variance and joint covariance |
| Q-11 | `simulate()` | Exact seeded Cholesky reconstruction |
| Q-12 | `sp_variogram()` | Direct pairwise group-difference contribution in every bin |
| Q-13--Q-14 | non-Gaussian `sp_fit()` | Direct latent covariance with fixed variance; positive fitted variance and Hessian schema |
| Q-15 | non-Gaussian prediction methods | Direct point mean/variance and full joint link covariance at known levels |
| Q-16 | `simulate()` | Exact seeded latent-covariance and binomial transform |
| Q-10a--Q-10c | `.sp_random_new()`, `.sp_random_new_covariance()`, `predict()`, `sp_predict_draws()` | Direct Gaussian and binomial marginal means, point variances and joint covariance for shared unseen levels |
| stochastic extensions | `sp_mcml()`, `sp_predict_draws()` | Fixed and estimated grouped MCML paths; finite seeded parameter draws with rebuilt covariance |
| Workflow retention | validation and reporting methods | Fold-wise refitting, simulation recovery, random-variance table rows, safe default failures and explicit marginal unseen-level prediction |

The primary checks are in `tests/test-random-effects.R` and
`tests/test-nongaussian-random-effects.R`.
