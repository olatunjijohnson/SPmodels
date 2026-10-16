# Decay-adjusted intervention histories and impact estimation

## Purpose, evidence and boundary

Repeated prevalence surveys may occur before and after interventions in the
same implementation unit. Treating intervention count as an ordinary linear
covariate does not express an immediate reduction followed by decay of that
reduction. This chapter derives the history adjustment and constrained impact estimation implemented
by `sp_intervention()`, `sp_intervention_multiplier()` and the binomial path in
`sp_fit()`.

The model is based on the mathematical specification in Giorgi, Fronterre and
Diggle (2026), *A decay-adjusted spatio-temporal model to account for the impact
of mass drug administration on neglected tropical disease prevalence*, JRSS A,
doi:10.1093/jrsssa/qnag100; the version checked during implementation was
arXiv:2512.03760v3. SPmodels' interface, derivation and code are independently
written.

The immediate reduction and decay time may be fixed or estimated jointly with
regression and spatial covariance parameters by either the labelled Laplace or
MCML route. A normalized beta penalty for the immediate reduction is optional.
Separable residual temporal correlation is available as derived in Chapter 15.

## Data, timing and assumptions

For survey row \(i=1,\ldots,n\), define:

- \(Y_i\): positive individuals among \(m_i\) examined;
- \(x_i\): fixed-effect covariates;
- \(t_i\): finite numeric survey time;
- \(g_i\): implementation-unit identifier;
- \(u_{gj}\): time of intervention \(j\) in unit \(g\);
- \(I_{gj}\in\{0,1\}\): whether that intervention was delivered.

Time differences and the decay parameter must use the same units. An
intervention contributes only when \(u_{gj}<t_i\). Consequently a survey and
intervention recorded at the same time are interpreted as survey first,
intervention second. Unit/time pairs in the history must be unique.

Without intervention, baseline prevalence is

\[
 p_i^*=\operatorname{logit}^{-1}(x_i^{\mathsf T}\beta+o_i+S_i),
 \qquad \boldsymbol S\sim N(\boldsymbol0,\Sigma_\psi).          \tag{I-01}
\]

Here \(o_i\) is a known formula offset and \(\Sigma_\psi\) is the existing
Matérn-plus-independent latent covariance.

## One intervention and the accumulated multiplier

Let \(\alpha\in[0,1]\) be the immediate proportional reduction and
\(\gamma>0\) the decay time. At lag \(v>0\), the remaining reduction is

\[
 f(v)=\alpha\exp(-v/\gamma).                                    \tag{I-02}
\]

The corresponding prevalence multiplier is (1-f(v)). For survey row (i),
SPmodels multiplies the contributions from all earlier delivered
interventions in the same unit:

\[
 A_i=\prod_{j:u_{g_i j}<t_i}
 \left[1-\alpha\exp\{-(t_i-u_{g_i j})/\gamma\}\right]^{I_{g_i j}}.
                                                                    \tag{I-03}
\]

An empty product is one. Since every lag in the product is strictly positive,
each factor is in ((0,1]), including when \(\alpha=1\); hence
\(0<A_i\le1\). The observed prevalence and sampling model are

\[
 p_i=A_i p_i^*,\qquad Y_i\mid S_i\sim\operatorname{Binomial}(m_i,p_i).
                                                                    \tag{I-04}
\]

This is not equivalent to adding \(\log A_i\) to a logit predictor because
\(A_i\operatorname{logit}^{-1}(\eta_i)\ne
\operatorname{logit}^{-1}(\eta_i+\log A_i)\).

## Latent score and curvature

Write \(q_i=\operatorname{logit}^{-1}(\eta_i)\), \(p_i=A_iq_i\), and
\(q_i'=q_i(1-q_i)\), where
\(\eta_i=x_i^{\mathsf T}\beta+o_i+S_i\). Apart from constants, observation
(i)'s log likelihood is

\[
 \ell_i=Y_i\{\log A_i+\log q_i\}
 +(m_i-Y_i)\log(1-A_iq_i).                                      \tag{I-05}
\]

Differentiating once gives

\[
 \frac{\partial\ell_i}{\partial\eta_i}
 =Y_i(1-q_i)
 -(m_i-Y_i)\frac{A_iq_i(1-q_i)}{1-A_iq_i}.                     \tag{I-06}
\]

For the second term in (I-06), the quotient derivative is

\[
 \frac{\partial}{\partial\eta_i}
 \frac{A_iq_i(1-q_i)}{1-A_iq_i}
 =\frac{A_iq_i(1-q_i)\{1-2q_i+A_iq_i^2\}}
 {(1-A_iq_i)^2}.                                                 \tag{I-07}
\]

Therefore the negative observation curvature is

\[
 w_i=Y_iq_i(1-q_i)+(m_i-Y_i)
 \frac{A_iq_i(1-q_i)\{1-2q_i+A_iq_i^2\}}
 {(1-A_iq_i)^2}.                                                 \tag{I-08}
\]

When \(A_i=1\), (I-06) becomes \(Y_i-m_iq_i\) and (I-08) becomes
\(m_iq_i(1-q_i)\), recovering ordinary binomial-logit expressions.

Let \(P=\Sigma_\psi^{-1}\). The posterior score and negative Hessian for the
latent vector are

\[
 g_i=\frac{\partial\ell_i}{\partial\eta_i}-(P\boldsymbol S)_i,
 \qquad H=P+\operatorname{diag}(w_1,\ldots,w_n).                 \tag{I-09}
\]

SPmodels uses the existing damped Newton iteration when \(H\) is positive
definite and accepts a step only when the log posterior increases. The direct
test differentiates (I-05) numerically and checks both (I-06) and (I-08).

## Laplace fitting and MCML

For trial values of \(\alpha\) and \(\gamma\), the latent mode
\(\widehat{\boldsymbol S}\), local covariance \(H^{-1}\), and Laplace
log likelihood are calculated exactly as in the non-Gaussian derivation, but
using (I-04)--(I-09). Regression and, when requested, covariance parameters are
then optimized. Fixed intervention parameters add no fitted degrees of freedom.

For MCML proposal draw \(\boldsymbol S_b\), the observation part of the log
weight becomes

\[
 \sum_i\log\operatorname{Binomial}
 \{Y_i;m_i,A_i\operatorname{logit}^{-1}(x_i^{\mathsf T}\beta+o_i+S_{bi})\}.
                                                                    \tag{I-10}
\]

The normalized Gaussian proposal, prior density, common-random-number rule and
diagnostics remain those in the MCML chapter.

When intervention parameters are estimated, MCML uses

\[
q_\alpha=\operatorname{logit}(\alpha),\qquad q_\gamma=\log\gamma. \tag{I-10a}
\]

It recalculates every \(A_i\) for every trial parameter vector and maximizes

\[
\widehat\ell_{MC}
=\log\left\{B^{-1}\sum_{b=1}^B w_b(\beta,\psi,\alpha,\gamma)\right\}
+\log\pi(\alpha),                                               \tag{I-10b}
\]

where \(\log\pi(\alpha)\) is zero when no beta penalty was declared. Adding
the normalized penalty to every log weight is algebraically identical to the
last term in (I-10b). Observed MCML curvature is transformed back with

\[
J=\operatorname{diag}\{\alpha(1-\alpha),\gamma\},\qquad
V_{\alpha\gamma}=J V_{q_\alpha q_\gamma}J^{\mathsf T}.          \tag{I-10c}
\]

## Prediction and simulation

For a prediction row, SPmodels first reconstructs \(A_{\star}\) from its unit,
survey time and the retained history. Let the approximate predictive link be
\(Z_{\star}\sim N(\mu_{\star},v_{\star})\). The latent response mean is

\[
 E(P_{\star})=A_{\star}E\{\operatorname{logit}^{-1}(Z_{\star})\}, \tag{I-11}
\]

evaluated by the same numerical Gaussian integration used for ordinary
binomial prediction. Its latent response variance uses

\[
 \operatorname{Var}(P_{\star})=A_{\star}^2
 \left[E\{\operatorname{logit}^{-1}(Z_{\star})^2\}
 -E\{\operatorname{logit}^{-1}(Z_{\star})\}^2\right].          \tag{I-12}
\]

For a future Bernoulli observation, the marginal variance is
\(E(P_{\star})\{1-E(P_{\star})\}\). Joint latent response draws transform each
Gaussian link draw \(z_{\star}^{(b)}\) as

\[
 p_{\star}^{(b)}=A_{\star}\operatorname{logit}^{-1}
 (z_{\star}^{(b)}).                                              \tag{I-13}
\]

When parameter uncertainty is requested, each Hessian draw contains
\((\beta,q_\alpha,q_\gamma)\) and any estimated covariance parameters. For
draw \(b\), SPmodels reconstructs

\[
\alpha^{(b)}=\operatorname{logit}^{-1}(q_\alpha^{(b)}),\qquad
\gamma^{(b)}=\exp(q_\gamma^{(b)}),\qquad
A_*^{(b)}=A_*(\alpha^{(b)},\gamma^{(b)}),                       \tag{I-13a}
\]

then rebuilds the fitted latent state and multiplies that draw's response
probability by \(A_*^{(b)}\). The returned multiplier matrix makes this
transformation auditable.

At fitted survey sites, unconditional simulation draws
\(\boldsymbol S^{(b)}\sim N(0,\Sigma_\psi)\) and then

\[
 Y_i^{(b)}\sim\operatorname{Binomial}
 \left[m_i,A_i\operatorname{logit}^{-1}
 \{x_i^{\mathsf T}\widehat\beta+o_i+S_i^{(b)}\}\right].      \tag{I-14}
\]

Cross-validation and simulation recovery retain the complete history while
excluding held-out responses from fitting. Thus intervention schedules are
treated as known design information, not outcomes that leak across folds.

## Algorithm and validation map

1. Validate and normalize unique unit/time intervention rows.
2. Match survey unit and numeric time columns by their declared names.
3. Apply the strict-before rule and calculate (I-03) without rounding.
4. Evaluate the adjusted latent mode through (I-05)--(I-09).
5. Fit regression, covariance and optionally constrained intervention
   parameters by the declared Laplace or normalized MCML criterion.
6. Recalculate multipliers for new prediction rows; apply (I-11)--(I-14).
7. Retain the normalized history, parameters, timing rule and row multipliers
   in auditable objects.

| Equations | Implementation | Independent checks |
|---|---|---|
| I-02--I-04 | history constructor and multiplier | Hand-calculated no-event, same-time, one-event, repeated-event and untreated cases |
| I-05--I-09 | `sp_latent_mode()` intervention path | Centred finite-difference score and Hessian; ordinary-binomial limiting case |
| I-10--I-10c | `sp_mcml()` intervention path | Direct response, prior, proposal and beta-penalty log weights; constraints, multiplier update and delta covariance |
| I-11--I-13a | `predict()`, `sp_predict_draws()` | Direct integration, exact seeded transform and parameter-draw-specific multiplier reconstruction |
| I-14 | `simulate()` | Exact seeded latent Gaussian and binomial reconstruction |
| I-15--I-17 | constrained Laplace optimizer and reporting | Parameter constraints, normalized penalty reconstruction, degrees of freedom and delta-method covariance |
| I-18--I-19 | effect-curve calculation and plot | Direct curve values, gradient variance, bounds and rendered display |
| Workflow | validation and reporting methods | Retained intervention specification, fitted/fixed table rows and explicit invalid-input failures |

## Constrained impact estimation and uncertainty

When estimation is requested, supplied \(\alpha\in(0,1)\) and \(\gamma>0\)
are starting values. SPmodels optimizes

\[
 \rho=\log\{\alpha/(1-\alpha)\},\qquad \xi=\log\gamma,
 \quad \alpha=\operatorname{logit}^{-1}(\rho),\quad \gamma=e^\xi.
                                                                    \tag{I-15}
\]

Finite computational bounds keep \(\alpha\) away from exact zero and one and
scale broad positive decay bounds to the observed history/survey time span.
If positive shapes \(a,b\) are supplied, the optimized criterion is

\[
 \ell_{\mathrm{pen}}=\ell_{\mathrm{Laplace}}+
 \log\left\{\frac{\alpha^{a-1}(1-\alpha)^{b-1}}{B(a,b)}\right\}.
                                                                    \tag{I-16}
\]

The beta normalizing constant is retained. The result is labelled *Penalized
Laplace*, not an unpenalized marginal likelihood. Let \(V_{\rho\xi}\) be the
inverse observed-Hessian block. Ordinary-scale uncertainty is

\[
 V_{\alpha\gamma}=J V_{\rho\xi}J^{\mathsf T},\qquad
 J=\operatorname{diag}\{\alpha(1-\alpha),\gamma\}.              \tag{I-17}
\]

An unusable Hessian produces missing standard errors rather than certainty.

## Intervention-effect curve and display

For a single intervention and lag \(v\ge0\), the displayed reduction and
prevalence multiplier are

\[
 r(v)=\alpha e^{-v/\gamma},\qquad m(v)=1-r(v).                  \tag{I-18}
\]

For an estimated fit,

\[
 \nabla r(v)=\left(e^{-v/\gamma},\quad
 r(v)\frac{v}{\gamma^2}\right)^{\mathsf T},\qquad
 \operatorname{Var}\{r(v)\}\approx
 \nabla r(v)^{\mathsf T}V_{\alpha\gamma}\nabla r(v).             \tag{I-19}
\]

Pointwise normal intervals are clipped to \([0,1]\); multiplier limits reverse
because \(m=1-r\). The returned table retains unrounded values, and the plot is
only a view over that table.

The primary tests are in `tests/test-intervention.R`. They reconstruct the
penalized objective, constraints, degrees of freedom, curve uncertainty and
display. A fixed or fitted intervention history can coexist with the separable
temporal covariance in Chapter 15. The MCML checks reconstruct every final
importance log weight and the fitted row multipliers.
