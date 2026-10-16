# Binomial and Poisson geostatistical fitting by Laplace approximation

## Scope and assumptions

This chapter derives the first fitted non-Gaussian models in SPmodels. It is
important to name the approximation: this implementation maximises a
**Laplace approximation** to the marginal likelihood. It is not Monte Carlo
maximum likelihood.

The observations are point supported. Conditional on a latent Gaussian vector,
they are independent. Coordinates lie in a two-dimensional Euclidean system
with declared units. The mean-model design matrix has full column rank.
Matérn smoothness is fixed.

## Linear predictor and covariance

For observation \(i\),

\[
\eta_i=\boldsymbol{x}_i^{\mathsf T}\boldsymbol\beta+o_i+s_i,      \tag{L-01}
\]

where \(o_i\) is a known offset. The latent vector has distribution

\[
\boldsymbol{s}\sim N(\boldsymbol 0,\boldsymbol\Sigma),           \tag{L-02}
\]

with entries

\[
\Sigma_{ij}=\sigma_s^2R(d_{ij};\rho,\nu)
 +\sigma_e^2\mathbb 1(i=j).                                     \tag{L-03}
\]

Here \(\sigma_s^2>0\) is spatial variance, \(\rho>0\) is Matérn
scale, fixed \(\nu>0\) is smoothness, and \(\sigma_e^2>0\) is an
independent latent-effect variance. Unlike Gaussian measurement noise, this
independent effect enters before the inverse link.

## Response likelihoods

For binomial-logit data,

\[
Y_i\mid s_i\sim\operatorname{Binomial}(m_i,p_i),\qquad
p_i=\{1+\exp(-\eta_i)\}^{-1}.                                    \tag{L-04}
\]

SPmodels accepts successes and failures as a two-column response, a proportion
with integer trial weights, or a zero-one Bernoulli response. Each encoding is
converted to successes \(y_i\) and trials \(m_i\), so it defines the same log
likelihood

\[
\ell_B(\boldsymbol s)=\sum_i\left[
\log {m_i\choose y_i}+y_i\log p_i+(m_i-y_i)\log(1-p_i)
\right].                                                         \tag{L-05}
\]

For Poisson-log data,

\[
Y_i\mid s_i\sim\operatorname{Poisson}(\mu_i),\qquad
\mu_i=\exp(\eta_i),                                              \tag{L-06}
\]

and

\[
\ell_P(\boldsymbol s)=\sum_i
\{y_i\log\mu_i-\mu_i-\log(y_i!)\}.                             \tag{L-07}
\]

Known exposure belongs in \(o_i=\log(e_i)\). Generic Poisson weights are not
silently reinterpreted as exposure.

## The marginal likelihood problem

Let \(\boldsymbol\theta=(\boldsymbol\beta,\sigma_s^2,\rho,
\sigma_e^2)\), and write \(\ell_y\) for (L-05) or (L-07). The exact
marginal likelihood integrates out the latent vector:

\[
L(\boldsymbol\theta)=\int
\exp\{\ell_y(\boldsymbol s)\}
\phi_n(\boldsymbol s;\boldsymbol 0,\boldsymbol\Sigma)
\,d\boldsymbol s.                                                \tag{L-08}
\]

For nonlinear links this integral generally has no closed form. Define

\[
h(\boldsymbol s)=\ell_y(\boldsymbol s)
-\frac12\boldsymbol s^{\mathsf T}\boldsymbol\Sigma^{-1}
\boldsymbol s.                                                    \tag{L-09}
\]

Chapter 4 derives the unique mode \(\widehat{\boldsymbol s}\) and negative
Hessian

\[
\boldsymbol H=\boldsymbol\Sigma^{-1}+\boldsymbol W              \tag{L-10}
\]

at that mode.

## Laplace approximation

A second-order Taylor expansion around the mode is

\[
h(\boldsymbol s)\approx h(\widehat{\boldsymbol s})
-\frac12(\boldsymbol s-\widehat{\boldsymbol s})^{\mathsf T}
\boldsymbol H
(\boldsymbol s-\widehat{\boldsymbol s}).                         \tag{L-11}
\]

Integrating the Gaussian quadratic in (L-11) and retaining every normalising
term that depends on the parameters gives

\[
\boxed{
\widetilde\ell(\boldsymbol\theta)=
\ell_y(\widehat{\boldsymbol s})
-\frac12\widehat{\boldsymbol s}^{\mathsf T}
 \boldsymbol\Sigma^{-1}\widehat{\boldsymbol s}
-\frac12\log|\boldsymbol\Sigma|
-\frac12\log|\boldsymbol H|.
}                                                                 \tag{L-12}
\]

The \(n\log(2\pi)/2\) term from the Gaussian prior cancels the corresponding
Laplace integration term. Since the latent-mode routine returns
\(\boldsymbol C=\boldsymbol H^{-1}\), the implementation evaluates

\[
-\frac12\log|\boldsymbol H|=\frac12\log|\boldsymbol C|.          \tag{L-13}
\]

All determinants and linear solves use Cholesky factors; the implementation
does not form a determinant by multiplying eigenvalues.

## Parameter optimisation

The regression coefficients are initialised by an independent generalized
linear-model calculation. Positive covariance parameters are represented on
the log scale:

\[
\boldsymbol\gamma=
(\log\sigma_s^2,\log\rho,\log\sigma_e^2).                        \tag{L-14}
\]

At every trial value \((\boldsymbol\beta,\boldsymbol\gamma)\), SPmodels:

1. constructs \(\boldsymbol\Sigma\);
2. finds \(\widehat{\boldsymbol s}\) by damped Newton iteration;
3. evaluates (L-12); and
4. returns its negative to a bounded optimiser.

If covariance parameters are fixed, only \(\boldsymbol\beta\) is optimised.
Otherwise all components are optimised jointly with fixed smoothness.
If a trial point cannot produce a positive-definite covariance, converged
latent mode or finite likelihood, the objective returns the finite sentinel
\(10^{20}\). A finite sentinel keeps bounded numerical differentiation from
overflowing when one side of a trial step is invalid; it is far above any
accepted objective in the supported data sizes. The final parameters must
still produce a usable state, and non-zero optimizer convergence is reported
rather than hidden.

## Reported uncertainty and fitted values

When the numerical Hessian of the approximated negative log likelihood is
positive definite, its inverse supplies a local parameter covariance matrix:

\[
\widehat{\operatorname{Var}}
\begin{pmatrix}\widehat{\boldsymbol\beta}\\
\widehat{\boldsymbol\gamma}\end{pmatrix}
=\left\{\nabla^2[-\widetilde\ell]\right\}^{-1}.                 \tag{L-15}
\]

The delta method transforms the covariance-parameter block back from log
scale, as in the Gaussian derivation. If curvature is not usable, standard
errors are reported as missing rather than repaired silently.

Fitted values are conditional response means at the latent mode:

\[
\widehat y_i=
\begin{cases}
m_i\operatorname{logit}^{-1}
(\boldsymbol{x}_i^{\mathsf T}\widehat{\boldsymbol\beta}+o_i+
\widehat s_i), & \text{binomial},\\
\exp(\boldsymbol{x}_i^{\mathsf T}\widehat{\boldsymbol\beta}+o_i+
\widehat s_i), & \text{Poisson}.
\end{cases}                                                       \tag{L-16}
\]

They are not marginal means integrated over latent uncertainty.

## Point prediction under the Laplace approximation

Let \(u_*\) be a latent effect at a new point and let
\(\boldsymbol c_*\) contain its spatial covariances with the fitted sites.
The independent latent component contributes no cross-covariance at a new
site. Conditional Gaussian identities give

\[
u_*\mid\boldsymbol s\sim N\left(
\boldsymbol c_*^{\mathsf T}\boldsymbol\Sigma^{-1}\boldsymbol s,
v_0-\boldsymbol c_*^{\mathsf T}\boldsymbol\Sigma^{-1}
\boldsymbol c_*
\right),                                                         \tag{L-17}
\]

where \(v_0=\sigma_s^2\) for the smooth latent target and
\(v_0=\sigma_s^2+\sigma_e^2\) for a new observation-site target. Combining
(L-17) with the Laplace posterior
\(\boldsymbol s\mid\boldsymbol y\approx
N(\widehat{\boldsymbol s},\boldsymbol C)\) yields an approximate normal
predictive linear predictor with mean

\[
m_* = \boldsymbol x_*^{\mathsf T}\widehat{\boldsymbol\beta}+o_*
+\boldsymbol c_*^{\mathsf T}\boldsymbol\Sigma^{-1}
\widehat{\boldsymbol s}                                          \tag{L-18}
\]

and variance

\[
v_* = v_0-\boldsymbol c_*^{\mathsf T}\boldsymbol\Sigma^{-1}
\boldsymbol c_*
+\boldsymbol c_*^{\mathsf T}\boldsymbol\Sigma^{-1}
\boldsymbol C\boldsymbol\Sigma^{-1}\boldsymbol c_*.             \tag{L-19}
\]

The last term propagates local latent-mode uncertainty. Covariance-parameter
and regression-parameter uncertainty are not included.

For a Poisson model, lognormal moments give the marginal intensity

\[
E(\mu_*\mid\boldsymbol y)=\exp(m_*+v_*/2),                        \tag{L-20}
\]

and its variance

\[
\operatorname{Var}(\mu_*\mid\boldsymbol y)=
E(\mu_*\mid\boldsymbol y)^2\{\exp(v_*)-1\}.                      \tag{L-21}
\]

For a future Poisson count, conditional sampling variation adds
\(E(\mu_*\mid\boldsymbol y)\) to (L-21).

For a binomial model, there is no corresponding elementary logistic-normal
mean. SPmodels evaluates

\[
\overline p_* = \int_{-\infty}^{\infty}
\operatorname{logit}^{-1}(m_*+\sqrt{v_*}z)\phi(z)\,dz             \tag{L-22}
\]

and

\[
E(p_*^2\mid\boldsymbol y) = \int_{-\infty}^{\infty}
\{\operatorname{logit}^{-1}(m_*+\sqrt{v_*}z)\}^2\phi(z)\,dz      \tag{L-23}
\]

by numerical quadrature. Their difference gives latent probability variance.
For a future Bernoulli outcome the predictive variance simplifies to

\[
\operatorname{Var}(Y_*\mid\boldsymbol y)=
\overline p_*(1-\overline p_*).                                  \tag{L-24}
\]

Thus response prediction integrates through the nonlinear link; it does not
merely apply the inverse link to (L-18).

## Limitations and validation map

Laplace accuracy can be poor for sparse binary data, small cluster sizes,
large latent variances or weakly identified covariance parameters. Optimiser
convergence does not establish approximation accuracy. Later chapters derive
fitted-model simulation, conditional residual diagnostics, MCML correction and
optional Hessian-based parameter uncertainty in joint prediction. Those
extensions do not alter the preliminary Laplace equations in this chapter.

| Equations | Implementation | Independent checks |
|---|---|---|
| L-01--L-07 | `.sp_response_data()` and model-frame construction | Equivalent binomial encodings, Poisson offset and invalid-input tests |
| L-08--L-13 | `.sp_laplace_state()` | Direct reconstruction of all four terms in (L-12) |
| L-14--L-15 | `.sp_fit_nongaussian()` | Positivity, convergence reporting and explicit missing-curvature policy |
| L-16 | `fitted.spmodels_fit()` | Family-specific range checks |
| L-17--L-24 | `.sp_predict_nongaussian()` | Direct Poisson lognormal moments and binomial numerical integration |

Row-permutation invariance is checked with a simultaneous permutation of
responses, predictors, and coordinates. Reference comparison and stochastic
full reference validation remains separate future work.
