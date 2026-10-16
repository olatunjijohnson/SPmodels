# Monte Carlo maximum likelihood for non-Gaussian spatial models

## Purpose and prerequisites

The Laplace fit replaces a difficult integral by a local quadratic
approximation. This chapter derives the optional Monte Carlo maximum-likelihood
(MCML) correction implemented by `sp_mcml()`. It estimates the same marginal
likelihood integral directly by normalized importance sampling. The method is
independently implemented from the equations below.

The reader should know probability densities, the multivariate normal
distribution, logarithms, sample means and elementary maximum likelihood. The
latent model and Laplace mode are derived in Chapters 4 and 5.

## The marginal likelihood problem

Let \(\boldsymbol y=(y_1,\ldots,y_n)^{\mathsf T}\), let
\(\boldsymbol s\) be the latent Gaussian vector, and write

\[
\boldsymbol s\sim N(\boldsymbol0,\boldsymbol\Sigma_{\boldsymbol\psi}),
\qquad
\eta_i=\boldsymbol x_i^{\mathsf T}\boldsymbol\beta+o_i+s_i.       \tag{C-01}
\]

Here \(\boldsymbol\psi\) contains spatial variance, scale and independent
latent variance, while Matérn smoothness remains fixed. Conditional on
\(\boldsymbol s\), observations are independent binomial-logit or Poisson-log.
For \(\boldsymbol\theta=(\boldsymbol\beta,\boldsymbol\psi)\), define the joint
density

\[
h_{\boldsymbol\theta}(\boldsymbol s)=
\left\{\prod_{i=1}^n p(y_i\mid s_i,\boldsymbol\beta)\right\}
\phi_n(\boldsymbol s;\boldsymbol0,
       \boldsymbol\Sigma_{\boldsymbol\psi}).                    \tag{C-02}
\]

The marginal likelihood is the integral

\[
L(\boldsymbol\theta;\boldsymbol y)=
\int_{\mathbb R^n}h_{\boldsymbol\theta}(\boldsymbol s)
\,d\boldsymbol s.                                                \tag{C-03}
\]

This integral usually has no closed form. Laplace fitting approximates it near
one mode. MCML instead evaluates it as an expectation.

## A normalized importance distribution

At the supplied preliminary fit, let \(\widehat{\boldsymbol s}_0\) be the
latent posterior mode and let \(\boldsymbol C_0\) be the inverse negative
posterior Hessian. For a positive proposal-scale constant \(a\), SPmodels uses
the fully normalized density

\[
q(\boldsymbol s)=
\phi_n(\boldsymbol s;\widehat{\boldsymbol s}_0,
       a^2\boldsymbol C_0).                                      \tag{C-04}
\]

This is a proposal distribution, not a claim that the posterior is exactly
Gaussian. If \(\boldsymbol R^{\mathsf T}\boldsymbol R=a^2\boldsymbol C_0\)
and \(\boldsymbol z_b\sim N(\boldsymbol0,\boldsymbol I_n)\), independent
proposal draws are

\[
\boldsymbol s_b=\widehat{\boldsymbol s}_0+
\boldsymbol R^{\mathsf T}\boldsymbol z_b,
\qquad b=1,\ldots,B.                                             \tag{C-05}
\]

One scale can be too narrow for the posterior tails or unnecessarily wide
near its centre. With declared distinct scales \(a_1,\ldots,a_K\), SPmodels
therefore also supports the normalized mixture

\[
q_{\mathrm{mix}}(\boldsymbol s)=\sum_{k=1}^K\alpha_k
\phi_n(\boldsymbol s;\widehat{\boldsymbol s}_0,
a_k^2\boldsymbol C_0),\qquad \sum_k\alpha_k=1.                 \tag{C-05a}
\]

It assigns \(B_k\) draws deterministically to component \(k\), balanced so
the counts differ by at most one, and sets \(\alpha_k=B_k/B\). Every draw is
weighted by the complete mixture density in (C-05a), not merely by the density
of the component that generated it. The deterministic-mixture estimator is

\[
\widehat L_B=\sum_{k=1}^K\alpha_k\frac1{B_k}
\sum_{b\in k}\frac{h_{\boldsymbol\theta}(\boldsymbol s_b)}
{q_{\mathrm{mix}}(\boldsymbol s_b)}.                            \tag{C-05b}
\]

Because \(\alpha_k=B_k/B\), (C-05b) reduces algebraically to the ordinary
mean in (C-08). Evaluation uses row-wise log-sum-exp over the component
densities. This preserves normalization and avoids tail-density underflow.

### Optional discarded-pilot adaptation

A local Gaussian approximation can have a useful shape but a poor centre or
dispersion. When `adapt = TRUE`, SPmodels first draws an independent pilot
sample \(\boldsymbol S^{(0)}_1,\ldots,\boldsymbol S^{(0)}_{B_0}\) from the
original normalized proposal \(q_0\). At the preliminary parameter vector
\(\widehat{\boldsymbol\theta}_0\), it calculates

\[
 u_b=\frac{h_{\widehat{\boldsymbol\theta}_0}
 (\boldsymbol S^{(0)}_b)}{q_0(\boldsymbol S^{(0)}_b)},\qquad
 \widetilde u_b=\frac{u_b}{\sum_{c=1}^{B_0}u_c}.               \tag{C-05c}
\]

The importance-weighted pilot moments are

\[
 \overline{\boldsymbol s}_w=\sum_{b=1}^{B_0}\widetilde u_b
 \boldsymbol S^{(0)}_b,                                      \tag{C-05d}
\]

\[
 \boldsymbol C_w=
 \frac{\sum_{b=1}^{B_0}\widetilde u_b
 (\boldsymbol S^{(0)}_b-\overline{\boldsymbol s}_w)
 (\boldsymbol S^{(0)}_b-\overline{\boldsymbol s}_w)^{\mathsf T}}
 {1-\sum_{b=1}^{B_0}\widetilde u_b^2}.                        \tag{C-05e}
\]

The denominator is the usual normalized-weight correction. It approaches
zero when one pilot draw carries almost all the weight, which is why the pilot
ESS is reported and the calculation fails explicitly if the covariance cannot
be identified. For a declared shrinkage weight \(0<\lambda\leq1\), the adapted
base moments are

\[
 \boldsymbol m_1=\lambda\widehat{\boldsymbol s}_0+
 (1-\lambda)\overline{\boldsymbol s}_w,\qquad
 \boldsymbol C_1=\lambda\boldsymbol C_0+(1-\lambda)
 \boldsymbol C_w.                                             \tag{C-05f}
\]

Positive \(\lambda\) regularizes the empirical covariance toward the positive-
definite Laplace covariance. The retained proposal replaces the centre and
base covariance in (C-05a) by \(\boldsymbol m_1\) and \(\boldsymbol C_1\),
while preserving the declared scale factors and complete mixture density.

The pilot draws are not reused in the likelihood estimate. Conditional on the
pilot, SPmodels draws a new independent main sample from the resulting
normalized density \(q_1\). Consequently

\[
 E\left\{\left.
 \frac{h_{\boldsymbol\theta}(\boldsymbol S)}
 {q_1(\boldsymbol S)}\right|\text{pilot}\right\}
 =\int h_{\boldsymbol\theta}(\boldsymbol s)\,d\boldsymbol s
 =L(\boldsymbol\theta;\boldsymbol y).                          \tag{C-05g}
\]

Thus adaptation changes the conditional variance, not the target integral.
It is opt-in because noisy pilot weights can make moment estimates worse. If
the pilot ESS fraction is below the declared `ess_warning` threshold, SPmodels
warns and retains the original proposal moments; lowering that threshold is an
explicit decision to trust a weaker pilot. The used/fallback decision and its
reason, original and adapted moments, pilot draws, pilot weights and pilot
diagnostics are retained so the decision can be audited.

Because \(q\) integrates to one and is positive wherever the Gaussian prior is
positive, multiply and divide (C-03) by \(q\):

\[
L(\boldsymbol\theta;\boldsymbol y)
=\int \frac{h_{\boldsymbol\theta}(\boldsymbol s)}{q(\boldsymbol s)}
q(\boldsymbol s)\,d\boldsymbol s
=E_q\{w_{\boldsymbol\theta}(\boldsymbol S)\},                    \tag{C-06}
\]

where

\[
w_{\boldsymbol\theta}(\boldsymbol s)=
\frac{h_{\boldsymbol\theta}(\boldsymbol s)}{q(\boldsymbol s)}.  \tag{C-07}
\]

Therefore the Monte Carlo likelihood estimator is

\[
\widehat L_B(\boldsymbol\theta)=
\frac1B\sum_{b=1}^B w_{\boldsymbol\theta}(\boldsymbol s_b).      \tag{C-08}
\]

No unknown normalizing constant has been discarded. Consequently `logLik()`
can report \(\log\widehat L_B\), rather than only a relative likelihood.

## The implemented log weight

For each proposal draw, SPmodels evaluates

\[
\begin{split}
\log w_{\boldsymbol\theta}(\boldsymbol s_b)
={}&\sum_{i=1}^n\log p(y_i\mid s_{bi},\boldsymbol\beta)\\
&-\frac12\left\{n\log(2\pi)+
\log|\boldsymbol\Sigma_{\boldsymbol\psi}|+
\boldsymbol s_b^{\mathsf T}
\boldsymbol\Sigma_{\boldsymbol\psi}^{-1}\boldsymbol s_b\right\}\\
&-\log q(\boldsymbol s_b).
\end{split}                                                       \tag{C-09}
\]

Cholesky solves compute the quadratic form; SPmodels does not explicitly
invert the covariance for every draw. If \(m=\max_b\log w_b\), the stable
log-likelihood calculation is

\[
\log\widehat L_B=m+
\log\left\{\frac1B\sum_{b=1}^B\exp(\log w_b-m)\right\}.          \tag{C-10}
\]

The same standard-normal draws are retained throughout one optimization. This
common-random-number rule makes the approximated objective a deterministic
function of \(\boldsymbol\theta\) during that call. Regression parameters are
optimized on their ordinary scale. Positive covariance parameters are
optimized on the log scale with the same finite numerical bounds as the
preliminary fitter.

If a trial covariance cannot be factorized or its fixed-draw log weights are
not finite, the minimized objective returns the finite sentinel \(10^{20}\).
This prevents numerical finite-difference overflow at an invalid boundary;
the final parameter vector must still yield a usable latent state, and the
optimizer convergence code and importance diagnostics remain visible.

## Importance-sampling diagnostics

Let normalized weights be

\[
\widetilde w_b=\frac{w_b}{\sum_{c=1}^B w_c}.
\]

The importance effective sample size and fraction are

\[
B_{\mathrm{eff}}=\frac{1}{\sum_b\widetilde w_b^2},
\qquad f_{\mathrm{eff}}=\frac{B_{\mathrm{eff}}}{B}.              \tag{C-11}
\]

They range from \(1\) to \(B\), and from \(1/B\) to \(1\), respectively.
SPmodels also records \(\max_b\widetilde w_b\). A small effective fraction or
one dominant weight means that too few draws determine the answer.

For independent raw weights \(w_1,\ldots,w_B\), let \(\bar w\) and \(s_w\)
be their sample mean and standard deviation. The relative Monte Carlo standard
error is

\[
r_{\mathrm{MC}}=\frac{s_w}{\sqrt B\,\bar w}.                     \tag{C-12}
\]

The delta method gives the same first-order value as the Monte Carlo standard
error of \(\log\widehat L_B\). Rescaling every weight by a common constant does
not change (C-11) or (C-12), so both are calculated safely after subtracting
the largest log weight.

For deterministic mixture sampling, the likelihood mean remains the same but
its variance estimate respects the proposal-component strata:

\[
\widehat{\operatorname{Var}}(\widehat L_B)=
\sum_{k=1}^K\alpha_k^2\frac{s_{w,k}^2}{B_k}.                  \tag{C-12a}
\]

SPmodels uses the square root of (C-12a) divided by \(\widehat L_B\) for the
mixture relative MCSE. Treating all stratified weights as identically
distributed would use the wrong variance formula.

These diagnostics detect many failures but cannot prove that the proposal has
covered every important part of the integrand. Results should be repeated with
larger \(B\), a changed proposal scale, and a refitted proposal when scientific
conclusions are sensitive.

## Independent replication

A diagnostic calculated from one importance sample cannot show how much the
reported optimum moves when that sample changes. Let
\(\widehat{\boldsymbol\theta}_r\) be the optimizer from independent replicate
\(r=1,\ldots,R\), with the proposal rule and sample size held fixed. For
parameter component \(j\), SPmodels reports

\[
 \overline\theta_j=\frac1R\sum_{r=1}^R\widehat\theta_{rj},
 \qquad
 s_{\mathrm{between},j}=
 \left\{\frac1{R-1}\sum_{r=1}^R
 (\widehat\theta_{rj}-\overline\theta_j)^2\right\}^{1/2}.       \tag{C-13}
\]

It also retains each replicate's log likelihood, effective sample size,
largest normalized weight and relative MCSE. This is a sensitivity analysis:
the replicate mean is not substituted for an MCML fit, and
\(s_{\mathrm{between},j}\) is not the sampling standard error of a parameter.
`sp_mcml_replicates()` returns all fitted objects so no convergence failure or
diagnostic is hidden by the summary.

When discarded-pilot adaptation is requested, the replicate diagnostics also
carry pilot ESS fraction, maximum weight, relative MCSE and whether the pilot
passed the gate in (C-05c)--(C-05f). Precision-escalation histories record the
same fields at every attempt. A fallback is therefore visible in the compact
history as well as in the retained fitted object.

Let \(\widehat{\boldsymbol V}_r\) be the fixed-objective Hessian covariance
from replicate \(r\). SPmodels additionally records

\[
\widehat{\boldsymbol V}_{\mathrm{MC}}=
\frac1{R-1}\sum_{r=1}^R
(\widehat{\boldsymbol\theta}_r-\overline{\boldsymbol\theta})
(\widehat{\boldsymbol\theta}_r-\overline{\boldsymbol\theta})^{\mathsf T},
                                                                  \tag{C-13a}
\]

\[
\widehat{\boldsymbol V}_{\mathrm{sampling}}=
\frac1R\sum_{r=1}^R\widehat{\boldsymbol V}_r,
\qquad
\widehat{\boldsymbol V}_{\mathrm{total}}=
\widehat{\boldsymbol V}_{\mathrm{sampling}}+
\widehat{\boldsymbol V}_{\mathrm{MC}}.                         \tag{C-13b}
\]

Diagonal square roots give sampling, Monte Carlo and combined standard errors
on the optimizer's parameter scale. SPmodels also transforms every replicate
to natural model parameters and applies the diagonal Jacobian containing 1 for
mean coefficients, the parameter value for log-scale components, and
\(\alpha(1-\alpha)\) for a logit intervention impact. Equations C-13a--C-13b
are then repeated on that interpretable scale. The Monte Carlo standard error
of the replicate mean is smaller by \(\sqrt R\), but the package does not silently
substitute that mean for a fitted model. If any replicate Hessian covariance
is unusable, sampling and total entries remain missing while empirical Monte
Carlo covariance is retained.

## Precision-directed sample-size escalation

For an initial sample size \(B_1\), maximum \(B_{\max}\), and growth factor
\(g>1\), define

\[
 B_{k+1}=\min\{B_{\max},\max(B_k+1,\lceil gB_k\rceil)\}.        \tag{C-14}
\]

Attempt \(k\) meets the declared targets when

\[
 r_{\mathrm{MC},k}\le r_\star
 \quad\hbox{and}\quad
 f_{\mathrm{eff},k}\ge f_\star,                               \tag{C-15}
\]

where \(r_\star>0\) and \(0<f_\star<1\). `sp_mcml_precision()` stops at the
first success. If no attempt succeeds, it evaluates \(B_{\max}\), returns that
fit, records that the targets were not met and warns. It never relaxes a target
to manufacture success.

Every attempt uses a new independent standard-normal sample. If proposal
refresh is requested, the preceding estimate supplies a newly calculated
latent mode \(\widehat{\boldsymbol s}_{k}\) and local covariance
\(\boldsymbol C_k\), giving

\[
 q_{k+1}(\boldsymbol s)=
 \phi_n(\boldsymbol s;\widehat{\boldsymbol s}_{k},
 a^2\boldsymbol C_k).                                          \tag{C-16}
\]

Without refresh, every attempt uses the original proposal. Refresh can improve
local overlap, but neither policy proves that all important regions were
visited. The returned history records \(B_k\), seed, both sides of (C-15),
largest weight, log likelihood and the stopping decision for every attempt.

## Curvature and downstream prediction

The observed Hessian of the fixed-draw objective supplies a local parameter-
covariance approximation when it is positive definite. It describes curvature
of that particular Monte Carlo objective. It does not automatically add the
simulation error in (C-12) to every parameter standard error. Independent
replicates are required before the separate and combined covariance estimates
in (C-13a)--(C-13b) are available.

After optimization, SPmodels recomputes the latent mode and local covariance at
the MCML estimates. Existing point and joint prediction can therefore run at
the corrected parameter estimates, but their conditional latent calculations
still use a local Gaussian approximation. MCML correction of parameter
estimation does not make downstream latent prediction exact.

Cross-validation of an MCML fit repeats both the preliminary Laplace fit and
the MCML correction inside every training fold. Held-out observations never
enter either stage.

## Algorithm

1. Require a fitted binomial-logit or Poisson-log spatial model.
2. Construct (C-04), or the declared mixture (C-05a), from its latent mode and
   local covariance.
3. If adaptation is requested, draw the independent pilot and diagnose it. If
   its ESS threshold passes, calculate (C-05c)--(C-05f) and construct the
   normalized adapted proposal; otherwise retain the original moments and
   record the fallback.
4. Discard the pilot from likelihood estimation and generate the independent,
   retained component-balanced sample. Equation (C-05g) preserves the target.
5. For every trial parameter vector, construct the Matérn covariance and
   evaluate (C-09)--(C-10).
6. Maximize the fixed-draw log likelihood using bounded log covariance
   parameters.
7. Recalculate final weights, (C-11)--(C-12), Hessian information and the local
   latent state.
8. Warn, without silently changing the criterion, when the effective-sample
   fraction is below the declared threshold.
9. For independent replication, repeat the preceding steps with independent
   draws and
   calculate (C-13)--(C-13b).
10. For precision escalation, apply (C-14)--(C-16), retaining every attempt and
   stopping only at the declared targets or maximum.

## Validation map

| Equations | Implementation | Independent checks |
|---|---|---|
| C-01--C-03 | existing non-Gaussian model state | Previously reconstructed likelihood and prior terms |
| C-04--C-10 | `sp_mcml()` and `.sp_mcml_log_weights()` | One-dimensional binomial integration; direct normalized mixture density, balanced allocation and finite weights |
| C-05c--C-05g | `.sp_mcml_adapt_proposal()` and `sp_mcml()` | Exact seeded pilot/main separation, weighted moments, positive-definite regularization, direct adapted-mixture density and full-shrinkage boundary |
| C-11--C-12a | `.sp_mcml_diagnostics()` | Bounds, stratified mixture error, direct likelihood reconstruction and exact seeded reproducibility |
| C-13--C-13b | `sp_mcml_replicates()` | Direct replicate means, MC covariance, mean MCSE and combined covariance; independent seeds and caller-seed preservation |
| C-14--C-16 | `sp_mcml_precision()` | First-attempt success, exact capped escalation, history and invalid-control failures |
| Complete algorithm | `sp_mcml()` | Binomial and Poisson fits; local-seed preservation; downstream prediction; MCML repeated inside folds |

## Scope and references

The implementation is dense and intended for small to moderate point datasets.
Each attempt may use one Gaussian proposal or a normalized same-centre mixture
of declared scales. Independent replicates, proposal refresh between attempts,
parameter-level Monte Carlo covariance and automatic sample-size escalation
are available. Optional within-attempt adaptation uses an independently drawn,
discarded pilot and regularized importance-weighted moments; it is not assumed
to improve every proposal. Bridge sampling and large-data sparse matrix methods
are not implemented.

The general Monte Carlo maximum-likelihood justification follows Geyer (1994),
*Journal of the Royal Statistical Society: Series B*, 56(1), 261--274,
doi:10.1111/j.2517-6161.1994.tb01976.x. Spatial generalized-linear-model use is
discussed by Christensen (2004), *Journal of Computational and Graphical
Statistics*, 13(3), 702--718, doi:10.1198/106186004X2525. These citations
support the statistical method; the derivation, interface and implementation
here are SPmodels' own.
