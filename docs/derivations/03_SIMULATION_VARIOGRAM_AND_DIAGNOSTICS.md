# Gaussian simulation, empirical variograms and residual diagnostics

## Purpose

This chapter derives the three assessment operations added to the Gaussian
geostatistical model: unconditional response simulation, an empirical residual
variogram and decorrelated residuals. It also states what each operation can and
cannot establish.

The notation follows the Gaussian model chapter:

\[
\boldsymbol{Y}=\boldsymbol{X}\boldsymbol{\beta}+\boldsymbol{o}
+\boldsymbol{S}+\boldsymbol{\varepsilon},
\qquad
\boldsymbol{Y}\sim N(\boldsymbol{\mu},\boldsymbol{V}),
\]

where

\[
\boldsymbol{\mu}=\boldsymbol{X}\boldsymbol{\beta}+\boldsymbol{o},
\qquad
\boldsymbol{V}=\sigma_s^2\boldsymbol{R}+\tau^2\boldsymbol{I}.
\]

## Unconditional simulation

Let \(\boldsymbol{V}=\boldsymbol{L}\boldsymbol{L}^{\mathsf T}\) be a
Cholesky factorisation and draw
\(\boldsymbol{z}\sim N(\boldsymbol{0},\boldsymbol{I})\). Define

\[
\boxed{
\boldsymbol{Y}^{(b)}=\widehat{\boldsymbol{\mu}}
+\widehat{\boldsymbol{L}}\boldsymbol{z}^{(b)}.
}
\tag{A-01}
\]

Because \(E(\boldsymbol{z})=\boldsymbol{0}\),

\[
E\{\boldsymbol{Y}^{(b)}\}=\widehat{\boldsymbol{\mu}}.              \tag{A-02}
\]

Using \(\operatorname{Var}(\boldsymbol{A}\boldsymbol{z})=
\boldsymbol{A}\operatorname{Var}(\boldsymbol{z})\boldsymbol{A}^{\mathsf T}\),

\[
\operatorname{Var}\{\boldsymbol{Y}^{(b)}\}
=\widehat{\boldsymbol{L}}\boldsymbol{I}
\widehat{\boldsymbol{L}}^{\mathsf T}
=\widehat{\boldsymbol{V}}.                                        \tag{A-03}
\]

Thus (A-01) has exactly the fitted model's response mean and covariance. Each
simulation column uses an independent standard-normal vector. The current
`simulate()` method generates outcomes at the observed sites only. It is
unconditional on the observed responses and conditional on fitted parameter
values; it does not propagate parameter uncertainty.

When a seed is supplied, SPmodels temporarily sets it and restores the caller's
previous random-number state. Without a supplied seed, ordinary R random-number
state advances normally.

## Empirical semivariance

For residuals \(e_i\) and \(e_j\) at sites separated by distance \(d_{ij}\),
the pairwise empirical semivariance is

\[
\boxed{
g_{ij}=\frac{1}{2}(e_i-e_j)^2.
}
\tag{A-04}
\]

Split the interval from zero to a chosen maximum distance into \(K\)
equal-width bins. Let \(P_k\) be the set of unordered pairs in bin \(k\), and
let \(n_k=|P_k|\). The binned empirical variogram is

\[
\boxed{
\widehat\gamma_k=\frac{1}{n_k}
\sum_{(i,j)\in P_k}g_{ij}.
}
\tag{A-05}
\]

SPmodels always returns \(n_k\). A semivariance based on few pairs is unstable,
so the pair count is part of the result rather than merely plotting metadata.
Empty bins have a missing empirical value.

The default residual is the mean-model residual

\[
\boldsymbol{e}_{\mathrm{mean}}
=\boldsymbol{y}-\boldsymbol{o}-
\boldsymbol{X}\widehat{\boldsymbol{\beta}}.                        \tag{A-06}
\]

This retains the spatial signal that the variogram is intended to reveal. The
optional conditional residual subtracts the fitted conditional spatial effect
as well. Its variogram answers a different diagnostic question: how much
dependence remains after spatial prediction?

## Fitted model variogram

For two distinct response observations under the stationary model,

\[
\operatorname{Var}(Y_i-Y_j)
=2\{\sigma_s^2+\tau^2-\sigma_s^2R(d_{ij})\}.
\]

The theoretical semivariogram is half this quantity:

\[
\boxed{
\gamma(d)=\tau^2+\sigma_s^2\{1-R(d)\},\qquad d>0.
}
\tag{A-07}
\]

At exactly zero lag the conventional value is \(\gamma(0)=0\), with a jump of
size \(\tau^2\) immediately to the right when a nugget is present. SPmodels
evaluates (A-07) at positive bin midpoints and places it beside (A-05). This is
a descriptive diagnostic comparison, not a second fitting criterion.

## Decorrelated residuals

Under the fitted Gaussian model, the mean-model residual vector should have
approximately covariance \(\boldsymbol{V}\). If
\(\boldsymbol{V}=\boldsymbol{L}\boldsymbol{L}^{\mathsf T}\), define

\[
\boxed{
\boldsymbol{z}_{\mathrm{res}}=
\boldsymbol{L}^{-1}
(\boldsymbol{y}-\boldsymbol{o}-
\boldsymbol{X}\widehat{\boldsymbol{\beta}}).
}
\tag{A-08}
\]

If the mean and covariance parameters were known and the model were correct,
then

\[
\operatorname{Var}(\boldsymbol{z}_{\mathrm{res}})
=\boldsymbol{L}^{-1}\boldsymbol{V}\boldsymbol{L}^{-\mathsf T}
=\boldsymbol{I}.                                                    \tag{A-09}
\]

The components would therefore be independent standard normal variables. In
practice parameters are estimated, so (A-09) is a diagnostic reference rather
than an exact finite-sample distribution. A normal quantile plot can reveal
strong skewness, heavy tails or outliers, but it cannot by itself prove model
adequacy.

## Projected geometry policy

For point `sf` input, SPmodels extracts the two geometry coordinates only after
checking that:

1. every geometry is a non-empty two-dimensional point;
2. a coordinate reference system is known;
3. the CRS is projected rather than longitude/latitude;
4. prediction data use the same CRS as the fitting data.

The package does not guess a projection or silently transform prediction data.
For ordinary data frames, the user supplies two coordinate columns and states
their shared unit. This policy makes every distance in (A-04)--(A-07) explicit.

## Algorithms and code map

| Operation | Algorithm | Implementation |
| --- | --- | --- |
| Response simulation | A-01 using the stored Cholesky factor | `simulate.spmodels_fit()` |
| Pair semivariance | A-04 for the lower triangle | `sp_variogram()` |
| Bin estimate | A-05 with equal-width bins and counts | `sp_variogram()` |
| Fitted variogram | A-07 at bin midpoints | `sp_variogram()` |
| Decorrelated residual | A-08 by triangular solve | `sp_diagnostics()` |
| Geometry checks | Explicit projected-point policy | `.sp_coordinates()` |

## Independent validation

The analytic tests reconstruct simulations from the identical seeded normal
draws and stored Cholesky factor, verify that every unordered pair is counted by
the variogram, and reconstruct (A-08) directly. Geometry tests fit and predict
with projected point data and reject longitude/latitude data, duplicate
coordinate specifications and mismatched prediction CRSs.

These tests validate the stated algorithms. Simulation assessment across many
replicates and broader reference comparisons remain future work. Formal
permutation and fitted-model variogram envelopes are derived separately in
derivation 11.
