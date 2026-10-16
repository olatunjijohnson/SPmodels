# Matérn correlation: derivation and implementation

## Question, scope and prerequisites

This chapter derives the unit-variance Matérn correlation implemented by
`matern_cor()`. It answers: given two point-supported observations separated by
a non-negative distance, how strongly are their latent spatial effects
correlated under a stationary, isotropic Matérn model?

The chapter assumes introductory probability, calculus and linear algebra.
Readers need to know that a covariance matrix must be symmetric and positive
semidefinite. The modified Bessel function of the second kind, written
`K_nu(x)`, is introduced where it is used.

This first slice implements correlation evaluation only. It does not yet
implement a variance multiplier, nugget or measurement-error term, covariance
derivatives, model fitting, prediction, or a practical-range conversion.

## Notation and assumptions

Let `s_i` and `s_j` be locations in a Euclidean coordinate system whose distance
unit is known. Define

- `d_ij = ||s_i - s_j|| >= 0`, the distance between the locations;
- `a > 0`, the Matérn scale in the same units as `d_ij`;
- `nu > 0`, the dimensionless smoothness parameter; and
- `x_ij = d_ij / a`, a dimensionless distance.

We assume stationarity, so covariance does not depend on absolute location, and
isotropy, so it depends on separation only through distance. Coordinates in
angular degrees must not be passed as though they were planar distances.

## Correlation definition

For positive distance, SPmodels uses the scale parameterisation

**(M-01)**

`rho(d; a, nu) = [2^(1-nu) / Gamma(nu)] x^nu K_nu(x)`,

where `x = d/a`, `Gamma` is the gamma function, and `K_nu` is the modified
Bessel function of the second kind. Calling `a` a scale rather than a range is
deliberate: different fields use incompatible range conventions.

### Value at zero

For `nu > 0`, the small-argument behaviour of the Bessel function is

`K_nu(x) ~ 2^(nu-1) Gamma(nu) x^(-nu)` as `x` tends to zero from above.

Substitution into (M-01) gives

`rho(d; a, nu)
 ~ [2^(1-nu) / Gamma(nu)] x^nu
    [2^(nu-1) Gamma(nu) x^(-nu)]
 = 1`.

Therefore the continuous definition is

**(M-02)** `rho(0; a, nu) = 1`.

The implementation assigns this value directly. It does not evaluate a
singular Bessel expression at zero.

### Units and rescaling

Because distance and scale enter only through their ratio,

**(M-03)** `rho(c d; c a, nu) = rho(d; a, nu)` for every `c > 0`.

Thus converting both metres to kilometres leaves the correlation unchanged.
Converting only one of them changes the model and is an input error that a
future geometry layer must prevent.

## Half-integer checks

Half-integer Bessel functions have elementary forms. For `nu = 1/2`,

`K_(1/2)(x) = sqrt(pi / (2x)) exp(-x)` and `Gamma(1/2) = sqrt(pi)`.

Putting these identities into (M-01),

`rho = [2^(1-1/2)/sqrt(pi)] x^(1/2)
       sqrt(pi/(2x)) exp(-x)
     = exp(-x)`.

Hence

**(M-04)** `rho(d; a, 1/2) = exp(-d/a)`.

The recurrence relations for half-integer Bessel functions similarly give

**(M-05)** `rho(d; a, 3/2) = (1 + x) exp(-x)`, and

**(M-06)** `rho(d; a, 5/2) = (1 + x + x^2/3) exp(-x)`.

These three identities provide calculations independent of the general Bessel
formula and are used in deterministic tests.

## From a correlation to a matrix

For locations `s_1, ..., s_n`, define `R_ij = rho(d_ij; a, nu)`. Euclidean
distance is symmetric, so `d_ij = d_ji` and consequently `R_ij = R_ji`.
Equation (M-02) gives `R_ii = 1`.

The Matérn family has a non-negative spectral density; by the spectral
characterisation of stationary covariance functions, its inverse Fourier
transform is positive definite. Therefore, for every real vector `z`,

**(M-07)** `z' R z >= 0`.

Finite-precision eigenvalues can be slightly negative, so tests allow only a
small predeclared numerical tolerance. This mathematical property applies when
the supplied distances come from a valid Euclidean configuration; an arbitrary
non-negative matrix need not be a valid distance matrix.

## Numerical algorithm

Directly evaluating (M-01) can overflow near zero because `K_nu(x)` diverges,
and can underflow at large `x` because it decays exponentially. R provides the
exponentially scaled value

**(M-08)** `K_scaled(x, nu) = exp(x) K_nu(x)`.

For positive `x`, take logarithms of (M-01) and use (M-08):

**(M-09)**

`log rho = (1-nu) log 2 - log Gamma(nu) + nu log x
           + log K_scaled(x, nu) - x`.

The implementation evaluates `log Gamma(nu)` with `lgamma()`, evaluates the
scaled Bessel function with `besselK(..., expon.scaled = TRUE)`, and exponentiates
only at the end. Underflow to zero at very large distance is a valid
floating-point representation of a correlation smaller than the available
precision. At extremely small positive arguments, a Bessel overflow is replaced
by the limiting value one only below machine epsilon. Any other non-finite
result is an explicit error. Values displaced a few floating-point units beyond
the theoretical interval are clipped to `[0, 1]`.

For the exact half-integer values `nu = 1/2, 3/2, 5/2`, implementation uses
(M-04)--(M-06) directly. This is algebraically identical to the generic
definition and avoids repeated Bessel evaluation during covariance
optimization. Every other positive smoothness retains the scaled-Bessel log
path.

Algorithm:

1. Validate that all distances are finite and non-negative.
2. Validate that `a` and `nu` are finite positive scalars.
3. Form `x = d/a`.
4. Assign one where `x = 0` using (M-02).
5. For positive `x`, use (M-04)--(M-06) at their exact half-integer
   smoothness values; otherwise evaluate (M-09).
6. Apply the documented extreme-small-argument recovery and mathematical
   bounds, preserving the input dimensions and names.

The algorithm is deterministic and uses no random numbers.

## Interpretation and limitations

Larger scale makes a fixed physical distance smaller on the dimensionless
scale and therefore usually increases correlation. Smoothness controls both
short-distance behaviour and differentiability of the associated Gaussian
field. It is not a variance parameter. A complete covariance model will later
combine this correlation with separately named spatial variance, nugget and
measurement-error components.

This implementation assumes distances have already been computed correctly.
It does not inspect coordinate reference systems, convert units, validate that
a supplied matrix is Euclidean, or determine whether points and polygons have
compatible support.

## Implementation and validation map

| Equations | Implementation | Test evidence |
| --- | --- | --- |
| M-01, M-02, M-08, M-09 | `R/matern.R`, `matern_cor()` | zero value and general evaluation |
| M-03 | same | joint distance/scale rescaling |
| M-04--M-06 | same | three half-integer closed forms |
| M-07 | same, applied elementwise to a distance matrix | symmetry, unit diagonal and eigenvalue bound |

The tests also cover dimension/name preservation and invalid distances, scales
and smoothness values. Reference-implementation comparison remains planned in
the separate pinned validation environment; no agreement claim is made here.

The public `matern_cor()` function performs the complete contract validation
above. During model fitting, distances and scalar covariance parameters have
already passed those checks before numerical optimization begins. The internal
covariance builder therefore uses the identical closed forms for
\(\nu=1/2,3/2,5/2\), and the public generic Bessel calculation otherwise,
without repeatedly scanning the whole distance matrix for the same input
conditions. This changes validation placement rather than the mathematical
function. Tests compare the internal and public results element by element for
specialized and generic smoothness values and require preservation of matrix
dimensions and labels.

## Worked calculation

Take two locations 500 metres apart, scale `a = 1000` metres and `nu = 1/2`.
Then `x = 500/1000 = 0.5`. By (M-04),

`rho = exp(-0.5) = 0.6065307` (rounded).

Changing both quantities to kilometres gives `d = 0.5`, `a = 1` and the same
answer, illustrating (M-03).

## Exercises

1. Use (M-04) to find the distance at which correlation is `0.05` when
   `a = 2 km`. Answer: solve `exp(-d/a) = 0.05`, giving
   `d = -2 log(0.05) km`.
2. Verify algebraically that (M-05) equals one at zero.
3. Explain why applying (M-01) to longitude differences measured in degrees and
   a scale measured in kilometres violates (M-03)'s unit requirement.
4. Construct three one-dimensional locations, calculate their distance matrix,
   and verify numerically that the resulting correlation matrix has no
   materially negative eigenvalue.

## References

- M. L. Stein (1999), *Interpolation of Spatial Data: Some Theory for
  Kriging*, Springer, sections on the Matérn class.
- C. E. Rasmussen and C. K. I. Williams (2006), *Gaussian Processes for
  Machine Learning*, MIT Press, covariance functions chapter.
- NIST Digital Library of Mathematical Functions, chapters 5 and 10, gamma and
  modified Bessel functions: <https://dlmf.nist.gov/>.

