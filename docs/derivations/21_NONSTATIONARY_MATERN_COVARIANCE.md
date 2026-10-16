# Non-stationary Matérn covariance with local range and anisotropy

## Question and scope

A stationary process has covariance that depends on the separation between two
locations, not on where the pair occurs. This chapter defines a valid
non-stationary covariance whose local range and directionality can change over
space. SPmodels currently implements the covariance construction only. It does
not estimate the local matrices, fit regression coefficients under this
covariance, or expose fitted-model prediction, diagnostics or cross-validation.

The construction follows Paciorek and Schervish (2006). Its validity is
important: replacing Euclidean distance by an arbitrary spatially varying
distance does not generally produce a positive-semidefinite covariance.

## Prerequisites, notation and assumptions

Assume multivariate normal distributions, matrix determinants and Cholesky
factorisation. There are `n` locations \(s_i\in\mathbb R^2\), and a symmetric
positive-definite matrix \(\Sigma_i\) is specified at each location. Its
entries have squared coordinate units. For an isotropic local range
\(\rho_i>0\), \(\Sigma_i=\rho_i^2 I_2\). More generally, its eigenvectors
specify local axes and the square roots of its eigenvalues specify local
length-scales.

The global Matérn smoothness \(\nu>0\) is fixed. The field is zero-mean and
unit variance for the correlation construction. A process variance and
independent measurement nugget are separate model components and are not
included in the returned matrix.

## Deriving the local distance

For a stationary anisotropic Gaussian kernel with scale matrix \(\Sigma\),
the squared distance between locations \(s_i,s_j\) is the Mahalanobis form

\[
q_{ij}(\Sigma)=(s_i-s_j)^\mathsf T\Sigma^{-1}(s_i-s_j).
\tag{NS-01}
\]

When local matrices differ, use the symmetric pair-average

\[
\bar\Sigma_{ij}=\frac{\Sigma_i+\Sigma_j}{2},\qquad
q_{ij}=(s_i-s_j)^\mathsf T\bar\Sigma_{ij}^{-1}(s_i-s_j).
\tag{NS-02}
\]

Both matrices in the average are positive definite, so their average is
positive definite and the quadratic form is non-negative. Numerically, write
\(\bar\Sigma_{ij}=LL^\mathsf T\) using its Cholesky factor. Solving
\(Lz=s_i-s_j\) gives \(q_{ij}=z^\mathsf Tz\) without explicitly forming a
matrix inverse.

## Determinant adjustment and covariance

The local Matérn term alone does not in general define a valid covariance
matrix. The determinant adjustment is

\[
a_{ij}=\frac{|\Sigma_i|^{1/4}|\Sigma_j|^{1/4}}
{|\bar\Sigma_{ij}|^{1/2}}.
\tag{NS-03}
\]

The non-stationary correlation is

\[
R_{ij}=a_{ij}\,M_\nu(\sqrt{q_{ij}}),
\tag{NS-04}
\]

where \(M_\nu(r)=2^{1-\nu}r^\nu K_\nu(r)/\Gamma(\nu)\), with
\(M_\nu(0)=1\). Because determinant is positive on positive-definite
matrices, \(a_{ij}\) is well-defined and symmetric. For \(i=j\),
\(\bar\Sigma_{ii}=\Sigma_i\), \(q_{ii}=0\), and \(a_{ii}=1\); hence
\(R_{ii}=1\).

Here is the validity argument in steps. For \(q\geq0\), the Matérn
correlation has the positive mixture representation

\[
M_\nu(\sqrt q)=\frac{1}{\Gamma(\nu)}
\int_0^\infty t^{\nu-1}
\exp\!\left(-t-\frac{q}{4t}\right)\,dt.
\tag{NS-05}
\]

The mixing density \(t^{\nu-1}e^{-t}/\Gamma(\nu)\) is non-negative and
integrates to one. It is therefore enough to show that, for each \(t>0\),
the determinant-adjusted Gaussian kernel is positive semidefinite.

Write \(h_{ij}=s_i-s_j\), \(A_{ij}=\Sigma_i+\Sigma_j\), and
\(\bar\Sigma_{ij}=A_{ij}/2\). Then
\(q_{ij}=2h_{ij}^\mathsf T A_{ij}^{-1}h_{ij}\), so the exponential in
(NS-05) is \(\exp[-h_{ij}^\mathsf T A_{ij}^{-1}h_{ij}/(2t)]\).
For \(u\in\mathbb R^d\), define the local Gaussian kernel

\[
g_{i,t}(u)=(\pi t)^{-d/4}|\Sigma_i|^{-1/4}
\exp\!\left[-\frac{1}{2t}(u-s_i)^\mathsf T
\Sigma_i^{-1}(u-s_i)\right],\qquad d=2.
\tag{NS-06}
\]

Completing the square in the exponent gives

\[
\begin{aligned}
&(u-s_i)^\mathsf T\Sigma_i^{-1}(u-s_i)
+(u-s_j)^\mathsf T\Sigma_j^{-1}(u-s_j)\\
&=(u-m_{ij})^\mathsf T
(\Sigma_i^{-1}+\Sigma_j^{-1})(u-m_{ij})
+h_{ij}^\mathsf T A_{ij}^{-1}h_{ij},
\end{aligned}
\]

where
\(m_{ij}=(\Sigma_i^{-1}+\Sigma_j^{-1})^{-1}
(\Sigma_i^{-1}s_i+\Sigma_j^{-1}s_j)\). Integrating the Gaussian term and
using
\(|\Sigma_i^{-1}+\Sigma_j^{-1}|=|A_{ij}|/(|\Sigma_i||\Sigma_j|)\)
then yields

\[
\int_{\mathbb R^d}g_{i,t}(u)g_{j,t}(u)\,du
=\frac{|\Sigma_i|^{1/4}|\Sigma_j|^{1/4}}
{|\bar\Sigma_{ij}|^{1/2}}
\exp\!\left[-\frac{1}{2t}
h_{ij}^\mathsf T A_{ij}^{-1}h_{ij}\right].
\tag{NS-07}
\]

This is an inner product, so its finite Gram matrix is positive semidefinite:
for any real \(c_1,\ldots,c_n\),
\(\sum_{i,j}c_ic_j\int g_{i,t}g_{j,t}
=\int(\sum_i c_i g_{i,t})^2\,du\geq0\). Equation (NS-05) expresses
(NS-04) as a positive mixture of these Gram matrices. Thus (NS-04) is
positive semidefinite for every finite set of positive-definite local
matrices. The construction itself, not an empirical eigenvalue check,
guarantees validity.

## Stationary limit and interpretation

If \(\Sigma_i=\Sigma\) for all locations, then \(a_{ij}=1\) and
\(\bar\Sigma_{ij}=\Sigma\). Thus

\[
R_{ij}=M_\nu\!\left(\sqrt{(s_i-s_j)^\mathsf T\Sigma^{-1}
(s_i-s_j)}\right),
\tag{NS-08}
\]

which is stationary anisotropic Matérn correlation. With
\(\Sigma=\rho^2I_2\), the argument reduces to Euclidean distance divided by
\(\rho\), matching `matern_cor(distance, scale = rho, smoothness = nu)`.

Where \(\Sigma_i\) changes slowly, the covariance has approximately the local
range and directional behaviour encoded by that matrix. The model remains
globally valid even when local matrices change abruptly, but abrupt changes
may be scientifically implausible. Covariance non-stationarity here concerns
local dependence shape; the process variance is constant. A spatially varying
variance model would require a separate, explicit variance function.

## Computation and numerical checks

`sp_nonstationary_cor(coords, local_covariance, smoothness)` requires finite
two-dimensional coordinates and an `n x 2 x 2` array of symmetric positive-
definite matrices. It validates every matrix, factors each pair-average by
Cholesky, evaluates (NS-02)--(NS-04), and returns a symmetric unit-diagonal
matrix. It does not silently project an invalid matrix to positive definite.
The pair loop uses quadratic storage and cubic pairwise work in the number of
locations for fixed dimension two. Large-data sparse or low-rank computation
is not implemented.

Tests compare the constant isotropic case directly with the stationary Matérn
function, verify symmetry, unit diagonal and positive-semidefinite eigenvalues
in varying-range/directional fixtures, and reject malformed, asymmetric or
non-positive-definite local matrices. These finite tests supplement but do not
replace the construction's validity argument. No numerical comparison with an
external reference implementation has yet been run.

## Reference

Paciorek, C. J. and Schervish, M. J. (2006). Spatial modelling using a new
class of nonstationary covariance functions. *Environmetrics*, 17, 483--506.
https://doi.org/10.1002/env.785
