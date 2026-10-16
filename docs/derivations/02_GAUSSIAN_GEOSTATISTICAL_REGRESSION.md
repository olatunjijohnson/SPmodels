# Gaussian geostatistical regression: derivation and algorithm

## Purpose and prerequisites

This chapter derives the first complete SPmodels fitting path. It assumes
introductory linear models, multivariate normal distributions, matrix
multiplication and logarithms. The implementation is independent and follows
the equations below directly.

The current path answers three questions:

1. How is a linear mean combined with spatially correlated residual variation?
2. How are regression and covariance parameters fitted by maximum likelihood?
3. How are predictions and their conditional uncertainties computed?

It does not yet cover non-Gaussian responses, anisotropy, spatio-temporal
dependence, repeated observations, restricted maximum likelihood or covariance
parameter uncertainty.

## Data, assumptions and notation

For locations \(s_1,\ldots,s_n\), let \(Y_i\) be a continuous response and let
the row vector \(x_i^T\) contain its mean-model predictors. Stack these into
the response vector \(y\) and the full-rank \(n\times p\) design matrix \(X\).
An optional known offset is the vector \(o\).

The model is

\[
Y = X\beta + o + S + \epsilon,                                      \tag{G-01}
\]

where

\[
S \sim N(0,\sigma_s^2 R),\qquad
\epsilon \sim N(0,\tau^2 I),\qquad S\perp\epsilon.                  \tag{G-02}
\]

The assumptions are:

- the conditional response is Gaussian with identity link;
- both supplied coordinate axes are Euclidean and use the explicitly recorded
  common unit;
- dependence is stationary and isotropic, so it depends only on distance;
- the Matérn smoothness \(\nu>0\) is fixed;
- \(\sigma_s^2>0\), \(\rho>0\), and \(\tau^2>0\);
- the observation-level errors are mutually independent and independent of
  the spatial process;
- the design matrix has full column rank.

Here \(\sigma_s^2\) is the spatial variance, \(\rho\) is the Matérn scale and
\(\tau^2\) is the nugget variance. The distinction matters: the latent process
has variance \(\sigma_s^2\), while an observation has marginal variance
\(\sigma_s^2+\tau^2\).

## From coordinates to a covariance matrix

For two coordinate vectors \(s_i=(s_{i1},s_{i2})^T\) and
\(s_j=(s_{j1},s_{j2})^T\), the Euclidean distance is

\[
d_{ij}=\sqrt{(s_{i1}-s_{j1})^2+(s_{i2}-s_{j2})^2}.                  \tag{G-03}
\]

Let \(M_\nu(d/\rho)\) be the unit-variance Matérn correlation derived in
`01_MATERN_CORRELATION.md`. Then

\[
R_{ij}=M_\nu(d_{ij}/\rho),                                         \tag{G-04}
\]

and, after adding the independent nugget,

\[
V=\operatorname{Var}(Y)=\sigma_s^2R+\tau^2I.                       \tag{G-05}
\]

Because the nugget is strictly positive, \(V\) is positive definite in exact
arithmetic even when locations are duplicated. The initial implementation
still requires at least two distinct locations because a spatial scale cannot
be learned or meaningfully checked when every distance is zero.

## Gaussian log-likelihood

Define the offset-adjusted response \(y_o=y-o\). From (G-01)--(G-05),

\[
Y\sim N(X\beta+o,V).                                                \tag{G-06}
\]

The multivariate normal density gives

\[
\ell(\beta,\theta)=
-\frac12\left\{n\log(2\pi)+\log|V|+
(y_o-X\beta)^TV^{-1}(y_o-X\beta)\right\},                         \tag{G-07}
\]

where \(\theta=(\sigma_s^2,\rho,\tau^2)^T\), with \(\nu\) fixed.

### Independent limiting model

With `spatial = FALSE`, the spatial process is absent and

\[
Y=X\beta+o+\varepsilon,\qquad
\varepsilon\sim N(0,\tau^2I).                                    \tag{G-06a}
\]

Equations G-07--G-09 then reduce to ordinary least squares, while maximum
likelihood estimates

\[
\widehat\tau^2_{\mathrm{ML}}
=\frac{1}{n}(y_o-X\widehat\beta)^T(y_o-X\widehat\beta).           \tag{G-06b}
\]

The divisor is \(n\), not \(n-p\), because this is the maximum-likelihood
variance estimator. This route supplies an independent limiting check and a
non-spatial comparison without pretending that a spatial fit is an ordinary
linear model.

For any fixed \(\theta\), differentiate (G-07) with respect to \(\beta\):

\[
\frac{\partial\ell}{\partial\beta}
=X^TV^{-1}(y_o-X\beta).                                            \tag{G-08}
\]

Setting (G-08) to zero produces the generalized least-squares estimator

\[
\widehat\beta(\theta)=
(X^TV^{-1}X)^{-1}X^TV^{-1}y_o.                                    \tag{G-09}
\]

Substituting (G-09) into (G-07) produces the profile log-likelihood
\(\ell_p(\theta)=\ell\{\widehat\beta(\theta),\theta\}\). SPmodels maximises
this ordinary maximum likelihood criterion over the logarithms of the three
positive covariance parameters. Log parameters enforce positivity because
\(\exp(a)>0\) for every finite \(a\). This is ML, not REML.

The default starting scale is the median positive pairwise distance. An
ordinary least-squares residual variance supplies starting variance values,
split 80% spatial and 20% nugget. User-supplied positive starting values take
precedence. Bounds are deliberately broad safeguards against numerical
overflow, not inferential constraints.

## Stable computation without an explicit inverse

Writing \(V^{-1}\) makes the algebra readable, but explicitly computing an
inverse is unnecessary and can amplify rounding error. Since \(V\) is positive
definite, its Cholesky factorisation is

\[
V=L L^T                                                        \tag{G-10}
\]

(or equivalently \(V=U^TU\) for the upper-triangular convention used by R).
To obtain \(V^{-1}b\), solve two triangular systems. Also,

\[
\log|V|=2\sum_{i=1}^n\log L_{ii}.                                  \tag{G-11}
\]

The same approach is used for \(X^TV^{-1}X\). A failed Cholesky factorisation
means the proposed numerical system is not usable; the optimiser rejects that
proposal. The implementation does not silently add jitter because doing so
would change the stated statistical model.
Rejection uses the finite objective sentinel \(10^{20}\), avoiding overflow in
bounded numerical differentiation near an invalid trial point. The final
parameter vector must still produce valid Cholesky factors.

Conditional on fitted covariance parameters,

\[
\operatorname{Var}(\widehat\beta\mid\widehat\theta)
=(X^TV^{-1}X)^{-1}.                                                \tag{G-12}
\]

The `vcov()` and coefficient standard errors report (G-12). They do not
propagate uncertainty in \(\widehat\theta\).

### Approximate covariance-parameter uncertainty

The optimiser works with
\(\eta=(\log\sigma_s^2,\log\rho,\log\tau^2)^T\). Let \(H_\eta\) be
the numerical Hessian of the negative profile log-likelihood at the optimum.
When it is finite and positive definite, the observed-information
approximation is

\[
\operatorname{Var}(\widehat\eta)\approx H_\eta^{-1}.               \tag{G-12a}
\]

For \(\theta=\exp(\eta)\), the delta-method Jacobian is
\(J=\operatorname{diag}(\widehat\theta)\). Therefore

\[
\boxed{
\operatorname{Var}(\widehat\theta)\approx JH_\eta^{-1}J^T.
}                                                                  \tag{G-12b}
\]

The square roots of the diagonal of (G-12b) appear in `summary()`. If
\(H_\eta\) is not positive definite, they are reported as missing rather than
created by silently altering the Hessian. These are local large-sample
approximations and do not include boundary corrections.

## Conditional spatial prediction

Consider a new location \(s_*\), predictor row \(x_*^T\), and offset \(o_*\).
Let

\[
c_*=\operatorname{Cov}(Y,S_*)
=\sigma_s^2\{M_\nu(d(s_1,s_*)/\rho),\ldots,
M_\nu(d(s_n,s_*)/\rho)\}^T.                                      \tag{G-13}
\]

The nugget is absent from (G-13): a future independent observation error is
not shared with the observed errors. The Gaussian conditioning identity gives
the universal-kriging predictor

\[
\widehat Y_*=x_*^T\widehat\beta+o_*+
c_*^TV^{-1}(y_o-X\widehat\beta).                                  \tag{G-14}
\]

For known \(\beta\), conditioning reduces the latent variance from
\(\sigma_s^2\) by \(c_*^TV^{-1}c_*\). Estimating \(\beta\) adds the
universal-kriging correction. Define

\[
q_*=x_*-X^TV^{-1}c_* .                                             \tag{G-15}
\]

Then the reported latent prediction variance is

\[
v_{*,\mathrm{latent}}=
\sigma_s^2-c_*^TV^{-1}c_*+
q_*^T(X^TV^{-1}X)^{-1}q_* .                                       \tag{G-16}
\]

For a future noisy observation, its independent nugget is added:

\[
v_{*,\mathrm{observation}}=v_{*,\mathrm{latent}}+\tau^2.          \tag{G-17}
\]

`predict(..., se.fit = TRUE, target = "latent")` uses (G-16), while
`target = "observation"` uses (G-17). Both condition on fitted covariance
parameters. Therefore these standard errors are not complete parameter-
uncertainty intervals and must not be described as such.

## Algorithm implemented

1. Construct the model frame using ordinary R formula rules, apply the subset
   and missing-value policy, and build \(X\).
2. Match those rows to the two numeric coordinate columns and verify finite
   values, explicit units, full-rank \(X\), and at least two distinct sites.
3. Compute the Euclidean distance matrix once using (G-03).
4. For fixed covariance values, proceed directly to step 6. Otherwise construct
   starting values and optimise the negative profile log-likelihood on the log
   parameter scale.
5. At every objective evaluation, build (G-05), factor it, calculate (G-09),
   and evaluate (G-07) using (G-10)--(G-11).
6. Store the final triangular factors and matrix products needed by standard
   methods.
7. For prediction, build cross-distances and (G-13), then evaluate
   (G-14)--(G-17) by triangular solves.

The direct dense factorisation requires approximately \(O(n^3)\) computation
and \(O(n^2)\) memory. Large-data approximations are future capabilities and
must be derived and validated separately.

## Independent checks and interpretation

`tests/test-gaussian-fit.R` checks this implementation against matrix equations
written directly with base R:

- generalized least squares (G-09);
- log-likelihood (G-07);
- conditional mean (G-14);
- latent and observation variances (G-16)--(G-17);
- invariance to row permutation;
- positivity and convergence of an estimated covariance fit;
- explicit failures for unsupported or ambiguous input.

The half-integer identity \(M_{1/2}(d/\rho)=\exp(-d/\rho)\) makes the test
covariance independently transparent. These are local mathematical checks, not
evidence of agreement with the separately managed reference implementation.

## Code map

| Equations or task | Implementation |
| --- | --- |
| G-03 | `.sp_dist()` in `R/fit.R` |
| G-04--G-05 | `.sp_covariance()` in `R/fit.R` and `matern_cor()` |
| G-07--G-12 | `.sp_gaussian_state()` and `sp_fit()` in `R/fit.R` |
| G-13--G-17 | `predict.spmodels_fit()` in `R/methods.R` |
| Input and optimiser controls | `sp_control()` in `R/control.R` |
| Direct checks | `tests/test-gaussian-fit.R` |

## Practical cautions

- Coordinates such as longitude and latitude in degrees are not automatically
  projected. Supply appropriate Euclidean coordinates and state their units.
- An estimated scale is expressed in the supplied coordinate unit and is the
  parameter in \(M_\nu(d/\rho)\), not a named practical range.
- A small dataset may weakly identify spatial variance, scale and nugget.
- A convergence code of zero is numerical evidence that the optimiser met its
  stopping rule, not proof of global optimality or a suitable scientific model.
- Residual diagnostics, covariance-parameter intervals and reference numerical
  comparisons remain required before research conclusions rely on this path.
