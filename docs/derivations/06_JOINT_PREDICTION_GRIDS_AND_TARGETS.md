# Joint prediction, grids and draw-based targets

## Why marginal standard errors are insufficient

A set of pointwise prediction means and standard errors does not describe how
prediction errors move together. Regional averages, maxima and exceedance
summaries depend on the joint distribution. SPmodels therefore constructs one
predictive mean vector and covariance matrix, draws all target locations
together, and applies transformations or aggregation separately to each joint
draw.

By default, results in this chapter are conditional on fitted covariance
parameters.  Section J-15--J-18 adds the optional local Hessian approximation
used when `parameter_uncertainty = TRUE`.  The non-Gaussian results additionally
inherit the Laplace approximation.

## Gaussian joint prediction

Let \(\boldsymbol y\) be observed Gaussian responses with covariance
\(\boldsymbol V\). At \(m\) new locations, let
\(\boldsymbol K_{**}\) be the \(m\times m\) spatial covariance matrix and
\(\boldsymbol C\) the \(n\times m\) covariance between observed and new
locations. Define

\[
\widehat{\boldsymbol\beta}=
(\boldsymbol X^{\mathsf T}\boldsymbol V^{-1}\boldsymbol X)^{-1}
\boldsymbol X^{\mathsf T}\boldsymbol V^{-1}
(\boldsymbol y-\boldsymbol o),                                   \tag{J-01}
\]

and

\[
\boldsymbol Q=\boldsymbol X_*
-\boldsymbol C^{\mathsf T}\boldsymbol V^{-1}\boldsymbol X.      \tag{J-02}
\]

The joint universal-kriging mean is

\[
\boxed{
\boldsymbol m_*=
\boldsymbol X_*\widehat{\boldsymbol\beta}+\boldsymbol o_*
+\boldsymbol C^{\mathsf T}\boldsymbol V^{-1}
(\boldsymbol y-\boldsymbol o-
\boldsymbol X\widehat{\boldsymbol\beta}).
}                                                                 \tag{J-03}
\]

The covariance for the smooth latent target is

\[
\boxed{
\boldsymbol S_*=
\boldsymbol K_{**}
-\boldsymbol C^{\mathsf T}\boldsymbol V^{-1}\boldsymbol C
+\boldsymbol Q
(\boldsymbol X^{\mathsf T}\boldsymbol V^{-1}\boldsymbol X)^{-1}
\boldsymbol Q^{\mathsf T}.
}                                                                 \tag{J-04}
\]

The first subtraction is uncertainty reduction from the observations. The
last term accounts for estimating the regression mean while treating
covariance parameters as fixed. For a future noisy observation, add nugget
variance to every diagonal entry:

\[
\boldsymbol S_{*,\mathrm{obs}}=
\boldsymbol S_*+\tau^2\boldsymbol I_m.                            \tag{J-05}
\]

For a non-spatial Gaussian fit, (J-03) reduces to the regression mean and
(J-04) reduces to
\(\boldsymbol X_*\widehat{\operatorname{Var}}(\widehat{\boldsymbol\beta})
\boldsymbol X_*^{\mathsf T}\).

## Non-Gaussian joint prediction

For binomial and Poisson fits, Chapter 5 gives the Laplace posterior

\[
\boldsymbol s\mid\boldsymbol y
\ \dot\sim\ N(\widehat{\boldsymbol s},\boldsymbol C_s).          \tag{J-06}
\]

Let \(\boldsymbol\Sigma\) be the fitted latent covariance at observed sites,
\(\boldsymbol C\) the spatial cross-covariance, and
\(\boldsymbol K_{**}\) the new-site spatial covariance. Applying total
expectation and total covariance gives predictive link mean

\[
\boxed{
\boldsymbol m_*=
\boldsymbol X_*\widehat{\boldsymbol\beta}+\boldsymbol o_*
+\boldsymbol C^{\mathsf T}\boldsymbol\Sigma^{-1}
\widehat{\boldsymbol s}
}                                                                 \tag{J-07}
\]

and covariance

\[
\boxed{
\boldsymbol S_*=
\boldsymbol K_{**}
-\boldsymbol C^{\mathsf T}\boldsymbol\Sigma^{-1}\boldsymbol C
+\boldsymbol C^{\mathsf T}\boldsymbol\Sigma^{-1}
\boldsymbol C_s\boldsymbol\Sigma^{-1}\boldsymbol C.
}                                                                 \tag{J-08}
\]

The last term propagates local latent posterior uncertainty. Adding
\(\sigma_e^2\boldsymbol I_m\) produces the new observation-site link
covariance. Regression- and covariance-parameter uncertainty are not included
in the default conditional calculation; the optional approximation is derived
below.

## Generating joint draws

For either model class, let \(\boldsymbol S_*\) denote the selected predictive
covariance and choose a matrix \(\boldsymbol A\) satisfying
\(\boldsymbol A\boldsymbol A^{\mathsf T}=\boldsymbol S_*\). For draw \(b\),

\[
\boldsymbol z^{(b)}\sim N(\boldsymbol 0,\boldsymbol I_m),\qquad
\boldsymbol\eta_*^{(b)}=
\boldsymbol m_*+\boldsymbol A\boldsymbol z^{(b)}.                 \tag{J-09}
\]

SPmodels uses a Cholesky factor when the covariance is positive definite. If
rounding leaves a positive-semidefinite matrix with tiny zero or negative
eigenvalues, it checks those values against a scale-dependent tolerance and
sets only tolerable negative values to zero. A materially negative eigenvalue
raises an error.

Gaussian response and link draws coincide. For non-Gaussian latent response
targets,

\[
p_*^{(b)}=\operatorname{logit}^{-1}(\eta_*^{(b)})                \tag{J-10}
\]

or

\[
\mu_*^{(b)}=\exp(\eta_*^{(b)}).                                  \tag{J-11}
\]

For future observations, SPmodels adds response sampling:

\[
Y_{*i}^{(b)}\sim\operatorname{Bernoulli}(p_{*i}^{(b)})            \tag{J-12}
\]

or

\[
Y_{*i}^{(b)}\sim\operatorname{Poisson}(\mu_{*i}^{(b)}).          \tag{J-13}
\]

A supplied seed is local to the function and does not alter the caller's
random-number stream.

## Optional fitted-parameter uncertainty

The default draws above condition on fitted regression and covariance
parameters. When a usable observed-Hessian covariance is available, SPmodels
can add a local asymptotic parameter layer. Let

\[
\boldsymbol\psi=
(\boldsymbol\beta^{\mathsf T},
\log\sigma_s^2,\log\rho,\log\sigma_e^2)^{\mathsf T}.              \tag{J-14}
\]

For a non-Gaussian Laplace fit, the joint numerical Hessian supplies
\(\widehat{\boldsymbol V}_\psi\). The approximation is

\[
\boldsymbol\psi^{(b)}\sim N(
\widehat{\boldsymbol\psi},
\widehat{\boldsymbol V}_\psi).                                  \tag{J-15}
\]

For a Gaussian fit, the available profile covariance and conditional
regression covariance do not contain their cross-covariance. SPmodels makes
the explicit approximation

\[
\widehat{\boldsymbol V}_\psi=
\begin{pmatrix}
\widehat{\operatorname{Var}}(\widehat{\boldsymbol\beta})&0\\
0&\widehat{\operatorname{Var}}(\widehat{\boldsymbol\gamma})
\end{pmatrix},                                                   \tag{J-16}
\]

where \(\boldsymbol\gamma\) contains log covariance parameters and its block
is obtained by the inverse delta method. This block-independence assumption is
reported rather than hidden.

For every parameter draw, covariance parameters are transformed back to the
positive scale. SPmodels reconstructs the Gaussian conditioning system or the
non-Gaussian latent mode, then generates one draw using (J-09). Consequently,
parameter and conditional predictive uncertainty are combined in the final
draw collection. An unusable Hessian covariance or a parameter draw that
cannot produce a valid conditioning system raises an error. Draws are not
silently truncated or repaired.

## Weighted and grouped targets

Suppose locations assigned to region \(g\) have non-negative raw weights
\(a_{ig}\). SPmodels requires positive total weight in every region and uses

\[
w_{ig}=\frac{a_{ig}}{\sum_{j\in g}a_{jg}}.                       \tag{J-17}
\]

The regional target for draw \(b\) is calculated before summarisation:

\[
T_g^{(b)}=\sum_{i\in g}w_{ig}D_i^{(b)},                          \tag{J-18}
\]

where \(D_i^{(b)}\) is the requested point-level draw. This order retains
spatial covariance. Combining marginal variances after the fact would not.

The reported Monte Carlo mean and standard deviation are

\[
\overline T_g=\frac1B\sum_{b=1}^B T_g^{(b)},\qquad
s_g=\left\{\frac1{B-1}\sum_{b=1}^B
(T_g^{(b)}-\overline T_g)^2\right\}^{1/2}.                       \tag{J-19}
\]

Quantiles are calculated directly from \(\{T_g^{(b)}\}\). For threshold
\(c\), the exceedance estimate is

\[
\widehat P(T_g>c)=\frac1B\sum_{b=1}^B
\mathbb 1\{T_g^{(b)}>c\}.                                      \tag{J-20}
\]

This is Monte Carlo uncertainty conditional on the fitted model; it is not a
frequentist confidence interval for the fitted parameters.

## Projected grid construction

`sp_grid()` operates only in a known projected CRS so that cell sizes have the
declared coordinate units. Its algorithm is:

1. validate a non-empty projected `sf` object and positive x/y spacing;
2. use the union of polygon inputs, an explicit polygon boundary, or the convex
   hull of point inputs;
3. reject non-polygonal or zero-area boundaries;
4. construct regular cell centres over the boundary bounding box;
5. retain centres intersecting the boundary; and
6. order retained points from top to bottom and then left to right, assigning
   stable `grid_id` values.

The function does not guess a projection, transform longitude/latitude, or
claim that a centre-point grid performs exact polygon integration.

## Point-to-polygon allocation

`sp_polygon_weights()` assigns each projected prediction point to a projected
polygon region using strict point-within-polygon membership. It requires the
same CRS, non-missing region identifiers and finite non-negative point weights.
Points assigned to multiple distinct region identifiers raise an error. Points
outside every polygon either raise an error or are explicitly dropped, as
selected by the caller.

The resulting indices, region identifiers and raw weights feed directly into
(J-17)--(J-20). This is centre-point allocation. It is appropriate for point
support or a grid-centre approximation. Exact fractional cell-area overlay for
coarse cells that cross regional boundaries is derived in derivation 11.

## Code and validation map

| Equations or algorithm | Implementation | Independent checks |
|---|---|---|
| J-01--J-05 | `.sp_joint_gaussian()` | Direct universal-kriging mean, complete covariance and seeded matrix transform |
| J-06--J-08 | `.sp_joint_nongaussian()` | Marginal variances equal pointwise prediction variances |
| J-09--J-13 | `sp_predict_draws()` | Local seed, support checks and exact seeded Gaussian transform |
| J-14--J-16 | parameter layer in `sp_predict_draws()` | Local reproducibility, recorded draws and family-specific fitted-state reconstruction |
| J-17--J-20 | `sp_predict_target()` | Direct weighted draws, normalisation, quantiles and exceedance frequencies |
| Grid algorithm | `sp_grid()` | CRS rejection, convex hull, stable ordering and degenerate-boundary failure |
| Polygon allocation | `sp_polygon_weights()` | Direct regional weighted reconstruction and CRS mismatch failure |

Exact bounded-memory pointwise grid prediction is derived in Chapter 18.
Chapter 19 derives the optional bounded-memory nearest-predecessor
approximation for joint draws; ordinary `sp_predict_draws()` remains the exact
route for the declared fitted predictive distribution. Exact fractional
cell-area weights are implemented and derived in derivation 11.
