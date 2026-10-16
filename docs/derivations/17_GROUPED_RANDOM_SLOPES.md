# Grouped random slopes and constrained covariance estimation

## Scientific question and current scope

Repeated measurements from the same clinic, school, village or instrument may
share more than a common level shift. Their relationship with a numeric
covariate can also vary by group. For example, clinics may differ in both
baseline prevalence and the change associated with age. SPmodels represents
that variation with a Gaussian random intercept and slope whose covariance is
declared explicitly with `sp_random_slope()`. The matrix may remain fixed or
serve as the starting value for constrained maximum-likelihood estimation.

This implementation supports one or more additive grouped components in all
three spatial response families. Each component may contain a slope alone or
a correlated intercept and slope. Fixed matrices remain explicit sensitivity
assumptions. Estimated matrices use a Cholesky parameterization that is
positive definite at every optimizer evaluation; their natural-scale entries,
delta-method standard errors and fitted degrees of freedom are reported.

Prerequisites are matrix multiplication, covariance, multivariate normal
distributions, generalized least squares, and the prediction derivations in
Chapters 2, 5, 6 and 13.

## Notation and model

Let there be \(n\) observations and \(R\) grouped slope components. Component
\(r\) has grouping value \(g_{ir}\), numeric slope value \(q_{ir}\), and
coefficient vector

\[
b_{r\ell}=\begin{cases}
(b_{0,r\ell},b_{1,r\ell})^{\mathsf T},&\text{intercept included},\\
b_{1,r\ell},&\text{slope only},
\end{cases}
\qquad b_{r\ell}\sim N(0,G_r).                              \tag{RS-01}
\]

Different levels \(\ell\), different components, the spatial process and the
nugget are mutually independent. The row vector multiplying the coefficient
is

\[
w_{ir}^{\mathsf T}=\begin{cases}(1,q_{ir}),&\text{intercept included},\\
q_{ir},&\text{slope only}.\end{cases}                         \tag{RS-02}
\]

For an intercept and slope,

\[
G_r=\begin{pmatrix}
\omega_{0r}^2&\omega_{01,r}\\
\omega_{01,r}&\omega_{1r}^2
\end{pmatrix}.                                                \tag{RS-03}
\]

SPmodels requires \(G_r\) to be finite, symmetric and positive definite.
Positive definiteness implies positive marginal variances and

\[
|\omega_{01,r}|<\omega_{0r}\omega_{1r},                       \tag{RS-04}
\]

so the implied correlation lies strictly between minus one and one.

The Gaussian response model is

\[
Y_i=x_i^{\mathsf T}\beta+o_i+S_i+
\sum_{r=1}^R w_{ir}^{\mathsf T}b_{r,g_{ir}}+\varepsilon_i.    \tag{RS-05}
\]

## Deriving the grouped covariance

The mean of a grouped contribution is zero. For observations \(i,j\), if
\(g_{ir}\ne g_{jr}\), their coefficient vectors are independent and the
covariance is zero. If they share a level, they multiply the same vector.
Using \(\operatorname{Cov}(Ab,Cb)=A\operatorname{Cov}(b)C^{\mathsf T}\),

\[
\begin{aligned}
K_{r,ij}
&=\operatorname{Cov}(w_{ir}^{\mathsf T}b_{r,g_{ir}},
                     w_{jr}^{\mathsf T}b_{r,g_{jr}})\\
&=\mathbb 1(g_{ir}=g_{jr})w_{ir}^{\mathsf T}G_rw_{jr}.
\end{aligned}                                                \tag{RS-06}
\]

For the two-coefficient case, expanding the multiplication gives

\[
K_{r,ij}=\mathbb 1(g_{ir}=g_{jr})
\{\omega_{0r}^2+(q_{ir}+q_{jr})\omega_{01,r}
+q_{ir}q_{jr}\omega_{1r}^2\}.                               \tag{RS-07}
\]

This intermediate expression is useful for hand checks. It shows why the
group contribution is not constant when slopes are present.

Let \(K_r\) be the matrix with entries (RS-06). Independence lets covariance
components add, so the integrated Gaussian covariance is

\[
V=\sigma_s^2R_\theta+\sum_r K_r+
\sum_a\omega_a^2Z_aZ_a^{\mathsf T}+\tau^2I.                  \tag{RS-08}
\]

The final sum permits the independent random intercepts from Chapter 13 to be
used simultaneously. Because every \(G_r\) is positive definite, each \(K_r\)
is positive semidefinite: it can be written as a block design matrix times
\(I\otimes G_r\) times its transpose. Adding it therefore preserves a valid
covariance model.

## Likelihood and computation

After integration, the exact Gaussian log likelihood and GLS coefficient
estimate are precisely equations Q-05 and Q-06 with \(V\) from (RS-08).
SPmodels forms each \(K_r\) directly from the equality matrix
\(\mathbb 1(g_{ir}=g_{jr})\) and the product \(W_rG_rW_r^{\mathsf T}\), then
adds the two element by element. The complete covariance is factorized by
Cholesky decomposition for solves and its log determinant.

When `estimate = FALSE`, the supplied \(G_r\) entries do not enter the
optimizer. Consequently there is no delta-method block for that fixed matrix,
and confidence intervals are not fabricated for a supplied assumption.

### Positive-definite estimation

For an estimated intercept-slope matrix, write

\[
L_r=\begin{pmatrix}
\exp(a_r)&0\\ c_r&\exp(d_r)
\end{pmatrix},\qquad G_r=L_rL_r^{\mathsf T}.                 \tag{RS-08a}
\]

The diagonal exponentials are strictly positive. Because \(L_r\) is
nonsingular, \(G_r\) is positive definite for every finite
\((a_r,c_r,d_r)\); the optimizer cannot step to a negative variance or an
impossible correlation. Multiplication gives the natural parameters

\[
\omega_{0r}^2=e^{2a_r},\qquad
\omega_{01,r}=e^{a_r}c_r,\qquad
\omega_{1r}^2=c_r^2+e^{2d_r}.                                \tag{RS-08b}
\]

A slope-only component uses \(G_r=e^{2a_r}\). The supplied covariance is
converted to its unique Cholesky factor to obtain optimizer starting values.
For an exact Gaussian model, these transformed parameters join the other
covariance parameters in the profiled likelihood. For a latent non-Gaussian
model they join the regression and covariance parameters in the Laplace or
fixed-draw MCML objective.

Let \(\psi_r=(a_r,c_r,d_r)^{\mathsf T}\) and
\(g_r=(\omega_{0r}^2,\omega_{01,r},\omega_{1r}^2)^{\mathsf T}\).
The Jacobian used to return Hessian uncertainty to the readable scale is

\[
J_r=\frac{\partial g_r}{\partial\psi_r^{\mathsf T}}=
\begin{pmatrix}
2e^{2a_r}&0&0\\
e^{a_r}c_r&e^{a_r}&0\\
0&2c_r&2e^{2d_r}
\end{pmatrix}.                                                \tag{RS-08c}
\]

If \(H^{-1}_{\psi_r}\) is the relevant block of usable inverse observed
curvature, the local natural-scale covariance is

\[
\operatorname{Var}(\widehat g_r)\approx
J_rH^{-1}_{\psi_r}J_r^{\mathsf T}.                            \tag{RS-08d}
\]

Variance intervals use a log transform; the covariance entry, which may be
negative, uses a local identity-scale interval. These remain local Wald
summaries and can be unreliable near a boundary. For identifiability, an
estimated intercept-slope component requires at least three observed groups
and at least two distinct slope values within every group. It cannot coexist
with a separately estimated random intercept for the same grouping variable,
and only one estimated slope component may use that grouping variable. A
slope-only component requires at least two groups and at least two non-zero
slope observations within each group, so its shared grouped covariance is
actually represented in the data.

For binomial-logit and Poisson-log models, the grouped contribution is on the
latent link scale. Replacing the latent covariance in Chapter 5 by (RS-08)
produces the implemented Laplace mode, curvature and objective. MCML uses the
same covariance in the Gaussian prior density; no hidden secondary grouping
calculation is used.

## Prediction at known and unseen levels

For training row \(i\) and prediction row \(j\), define

\[
C_{r,ij}=\mathbb 1(g_{ir}=g_{jr}^{*})
w_{ir}^{\mathsf T}G_rw_{jr}^{*}.                             \tag{RS-09}
\]

The new-to-new component is

\[
K_{r,jk}^{*}=\mathbb 1(g_{jr}^{*}=g_{kr}^{*})
(w_{jr}^{*})^{\mathsf T}G_rw_{kr}^{*}.                       \tag{RS-10}
\]

These terms are added to the spatial cross and new covariance before applying
the Gaussian universal-kriging equations Q-09--Q-10 or the non-Gaussian
approximation Q-15.

An unseen level is rejected by default. With `unseen = "marginal"`, it is a
fresh zero-mean coefficient vector independent of all training groups.
Therefore its column in (RS-09) is zero. Equation (RS-10) still applies among
prediction rows: two rows with the same new label share one new coefficient
vector, while distinct new labels are independent. The slope values of both
rows remain in the covariance, so shared-new-level dependence is generally
not constant.

The same policy is explicit in cross-validation. If a group-constrained fold
withholds an entire fitted grouping level, the default prediction error is a
useful guard against silently treating it as known. Passing
`unseen = "marginal"` to `sp_cross_validate()` instead evaluates the declared
new-group prediction problem and carries that policy into integrated
non-Gaussian log scores.

## Simulation and variograms

Unconditional simulation uses one draw from the complete covariance
(RS-08). This is distributionally identical to drawing each group coefficient
vector separately and adding \(w_{ir}^{\mathsf T}b_{r,g_{ir}}\), because sums
of independent Gaussian vectors remain Gaussian and have additive covariance.

For any pair, the grouped semivariance contribution is

\[
\gamma_{r,ij}=\tfrac12(K_{r,ii}+K_{r,jj})-K_{r,ij}.           \tag{RS-11}
\]

It depends on group membership and both covariate values, not only distance.
SPmodels therefore evaluates (RS-11) for the actual observation pairs and
averages it within each variogram distance bin.

## Small worked calculation

Suppose two observations in clinic A have slope values 1 and 2 and

\[
G=\begin{pmatrix}4&1\\1&9\end{pmatrix}.
\]

Their row vectors are \(w_1=(1,1)\) and \(w_2=(1,2)\). Equation (RS-06) gives

\[
K_{12}=(1,1)G(1,2)^{\mathsf T}
=(1,1)(6,19)^{\mathsf T}=25.
\]

The individual variances are \(K_{11}=15\) and \(K_{22}=44\), so their grouped
semivariance is \((15+44)/2-25=4.5\). A row in clinic B has covariance zero
with both, even if its slope value is identical.

## Assumptions, checks and limitations

- Group labels and numeric slope values must be observed without missingness.
- Covariates should be centered when the intercept-slope covariance is to have
  a useful interpretation; shifting the zero of a slope changes the meaning
  of the random intercept and their covariance.
- A fixed \(G_r\) is a sensitivity assumption. Analysts should repeat fits over
  scientifically defensible matrices when it is uncertain. An estimated
  \(G_r\) needs enough independent groups and within-group covariate variation;
  passing the minimum checks does not by itself prove strong identification.
- The current method has one numeric slope per component, no random nonlinear
  basis, and no non-spatial fit.
- Agreement tests independently construct (RS-06)--(RS-10), the GLS estimate,
  likelihood, joint prediction covariance, Laplace covariance and MCML state.
- For fitted-parameter predictive draws, every sampled Cholesky parameter set
  rebuilds both the training and prediction-side copies of (G_r); using the
  fitted-centre matrix for the new-site block would mix two parameter states.

## Practice

1. Expand (RS-06) to obtain (RS-07).
2. Show that a slope-only model has covariance
   \(\mathbb 1(g_i=g_j)q_iq_j\omega_1^2\).
3. Explain why a new level changes predictive variance but not conditional
   mean under marginal prediction.
4. Recalculate the worked semivariance after centering the slope values at
   1.5, remembering that a scientifically equivalent covariance matrix would
   also need reparameterisation.

## Code and validation map

| Equations | Implementation | Independent check |
|---|---|---|
| RS-01--RS-04 | `sp_random_slope()` | dimensions, symmetry and Cholesky failures |
| RS-05--RS-08 | `.sp_random_slope_training()`, `.sp_random_slope_covariance()`, `sp_fit()` | direct full covariance, GLS and log likelihood |
| RS-08a--RS-08d | Cholesky parameter helpers, Gaussian/Laplace/MCML fitting, `confint()` | reconstruction, positive definiteness, analytic-versus-numerical Jacobian, natural-scale uncertainty and fitted degrees of freedom |
| RS-09--RS-10 | `predict()`, `sp_predict_draws()` | direct cross and joint covariance; shared unseen levels |
| parameter-state propagation | `.sp_prediction_random_slope()` | deterministic perturbed Cholesky state appears identically in the prediction-side covariance |
| RS-08 in latent models | Laplace fitting and `sp_mcml()` | exact latent-prior covariance equality |
| RS-11 | `sp_variogram()` | direct pairwise semivariance contribution |
| reporting contract | `summary()`, `sp_model_table()` | fixed versus estimated status, natural-scale standard errors and likelihood degrees of freedom |

The main executable checks are in `tests/test-random-slopes.R`.
