# Matérn scale derivatives and non-Gaussian latent modes

## Why these operations are needed

Binomial and Poisson geostatistical regression introduce a latent Gaussian
spatial vector inside a nonlinear response model. Later parameter estimation
needs covariance derivatives, a stable conditional latent mode and the
curvature at that mode. This chapter derives those deterministic building
blocks. It does not claim that non-Gaussian model fitting is complete.

## Matérn notation

Write the SPmodels Matérn correlation as

\[
R(d;\rho,\nu)=C_\nu x^\nu K_\nu(x),
\qquad
C_\nu=\frac{2^{1-\nu}}{\Gamma(\nu)},
\qquad
x=\frac{d}{\rho}.                                                  \tag{N-01}
\]

Distance \(d\) and smoothness \(\nu\) are held fixed while differentiating
with respect to scale \(\rho\). Since

\[
\frac{dx}{d\rho}=-\frac{x}{\rho},                                 \tag{N-02}
\]

the scale derivative follows from the Bessel identity

\[
\frac{d}{dx}\{x^\nu K_\nu(x)\}=-x^\nu K_{\nu-1}(x).               \tag{N-03}
\]

Combining (N-01)--(N-03),

\[
\boxed{
\frac{\partial R}{\partial\rho}
=\frac{C_\nu}{\rho}x^{\nu+1}K_{\nu-1}(x)
=\frac{R}{\rho}x\frac{K_{\nu-1}(x)}{K_\nu(x)}.
}                                                                  \tag{N-04}
\]

At \(d=0\), correlation is exactly one for every scale, so both scale
derivatives are zero.

## Second scale derivative

Differentiate (N-04), using

\[
\frac{dK_a(x)}{dx}=-\frac12\{K_{a-1}(x)+K_{a+1}(x)\}              \tag{N-05}
\]

and the Bessel recurrence relation. After collecting terms,

\[
\boxed{
\frac{\partial^2R}{\partial\rho^2}
=\frac{R}{\rho^2}
\left[
x^2-(2\nu+1)x\frac{K_{\nu-1}(x)}{K_\nu(x)}
\right].
}                                                                  \tag{N-06}
\]

For \(\nu=1/2\), \(R=\exp(-d/\rho)\), so the independently derived closed
forms are

\[
\frac{\partial R}{\partial\rho}
=R\frac{d}{\rho^2},
\qquad
\frac{\partial^2R}{\partial\rho^2}
=R\left(\frac{d^2}{\rho^4}-\frac{2d}{\rho^3}\right).              \tag{N-07}
\]

These provide sign-sensitive analytic tests. General smoothness values are
also compared with centred finite differences using predeclared tolerances.
The implementation evaluates ratios of exponentially scaled Bessel functions,
so their common exponential factor cancels.

## Non-Gaussian latent model

Let \(s\) be a latent Gaussian vector with

\[
s\sim N(0,\Sigma),\qquad \eta=a+s,                                \tag{N-08}
\]

where \(a=X\beta+o\) is supplied to the mode calculation. Conditional on
\(s\), observations are independent.

For binomial-logit data,

\[
Y_i\mid s_i\sim\operatorname{Binomial}(m_i,p_i),
\qquad p_i=\frac{1}{1+\exp(-\eta_i)}.                              \tag{N-09}
\]

For Poisson-log data,

\[
Y_i\mid s_i\sim\operatorname{Poisson}(\mu_i),
\qquad \mu_i=\exp(\eta_i).                                       \tag{N-10}
\]

Ignoring constants that do not depend on \(s\), the latent log posterior is

\[
h(s)=\sum_i\ell_i(\eta_i)-\frac12s^T\Sigma^{-1}s.                 \tag{N-11}
\]

## Score and curvature

For both supported families, define the conditional mean vector \(\mu\). The
likelihood score with respect to \(\eta\) is \(y-\mu\). Therefore

\[
\boxed{
g(s)=\frac{\partial h}{\partial s}
=y-\mu-\Sigma^{-1}s.
}                                                                  \tag{N-12}
\]

Let \(W\) be diagonal. Its entries are

\[
W_{ii}=m_ip_i(1-p_i)\quad\text{for binomial data},                 \tag{N-13}
\]

and

\[
W_{ii}=\mu_i\quad\text{for Poisson data}.                         \tag{N-14}
\]

The negative Hessian is

\[
\boxed{
H(s)=-\frac{\partial^2h}{\partial s\partial s^T}
=W+\Sigma^{-1}.
}                                                                  \tag{N-15}
\]

Because \(\Sigma^{-1}\) is positive definite and \(W\) is non-negative,
\(H(s)\) is positive definite. The log posterior is consequently strictly
concave and has one mode.

## Damped Newton algorithm

Starting from zero or a supplied vector, solve

\[
H(s^{(k)})\Delta^{(k)}=g(s^{(k)})                                 \tag{N-16}
\]

by Cholesky factorisation, then propose

\[
s^{(k+1)}=s^{(k)}+\alpha\Delta^{(k)},                              \tag{N-17}
\]

where \(\alpha\) begins at one. If (N-17) decreases (N-11), halve \(\alpha\)
until the objective does not decrease. This damping protects against an
over-large Newton step when the starting point is far from the mode.

Convergence is declared when

\[
\max_i|\alpha\Delta_i|
\leq \epsilon\{1+\max_i|s_i|\}.                                  \tag{N-18}
\]

At the final point, the returned local Gaussian covariance is

\[
\boxed{C_{\mathrm{mode}}=H(\widehat s)^{-1}.}                     \tag{N-19}
\]

This covariance describes local curvature. Returning it does not turn the
calculation into Bayesian parameter estimation and does not supply MCML by
itself.

## Metropolis-adjusted Langevin sampling

Let \(C=H(\widehat s)^{-1}\) be the fixed mode-based preconditioner from
(N-19), let \(g(s)\) be (N-12), and choose step size \(\delta>0\). From current
state \(s\), MALA proposes

\[
\boxed{
s' = s + \frac{\delta^2}{2}C g(s) + \delta C^{1/2}z,
\qquad z\sim N(0,I).
}                                                                  \tag{N-20}
\]

Thus the forward proposal density is

\[
q(s'\mid s)=N\!\left(
s';\ s+\frac{\delta^2}{2}Cg(s),\ \delta^2C
\right).                                                           \tag{N-21}
\]

Because the drift depends on the current state, this proposal is asymmetric.
The Metropolis correction is essential. Accept \(s'\) with probability

\[
\boxed{
\alpha(s,s')=\min\left\{1,
\frac{\exp\{h(s')\}\,q(s\mid s')}
     {\exp\{h(s)\}\,q(s'\mid s)}
\right\}.
}                                                                  \tag{N-22}
\]

SPmodels evaluates (N-22) on the log scale. The common proposal determinant
cancels because \(C\) and \(\delta\) are fixed. For proposal mean
\(m(s)=s+\delta^2Cg(s)/2\), the remaining forward quadratic is

\[
Q(s',s)=\frac{1}{\delta^2}
\{s'-m(s)\}^TC^{-1}\{s'-m(s)\}.                                   \tag{N-23}
\]

The log acceptance ratio is consequently

\[
h(s')-h(s)-\frac12Q(s,s')+\frac12Q(s',s).                          \tag{N-24}
\]

The sampler starts at the latent mode, discards the declared burn-in and keeps
every declared thinning interval. A supplied random seed is local to the call.
The result reports the realised acceptance rate; it does not automatically
retune a poorly chosen step size or certify convergence.

## Autocorrelation and effective sample size

Retained MCMC draws are dependent. For one latent component, let
\(\rho_k\) denote the sample autocorrelation at lag \(k\). If \(N\) draws are
retained, the ideal infinite-series definition of integrated autocorrelation
time is

\[
\tau_{\mathrm{int}}=1+2\sum_{k=1}^{\infty}\rho_k.                \tag{N-25}
\]

The corresponding effective sample size is

\[
\boxed{N_{\mathrm{eff}}=\frac{N}{\tau_{\mathrm{int}}}.}          \tag{N-26}
\]

An infinite noisy sum is not usable directly. SPmodels groups the estimated
autocorrelations into adjacent pairs,

\[
P_j=\widehat\rho_{2j-1}+\widehat\rho_{2j},                       \tag{N-27}
\]

and stops immediately before the first non-positive pair. If \(J\) positive
pairs remain, it uses

\[
\widehat\tau_{\mathrm{int}}
=\max\left\{1,1+2\sum_{k=1}^{2J}\widehat\rho_k\right\},
\qquad
\widehat N_{\mathrm{eff}}
=\min\left\{N,\frac{N}{\widehat\tau_{\mathrm{int}}}\right\}.    \tag{N-28}
\]

For one chain, `sp_sampler_diagnostics()` reports this estimate for every latent
component and trace and autocorrelation values for one selected component. A
satisfactory trace, acceptance rate or effective sample size is not proof of
convergence, which motivates the multiple-chain extension below.

## Multiple chains and dispersed initial states

Let $M$ be the requested number of chains. Each chain uses the same posterior
and MALA transition but evolves with its own sequence of random proposals and
acceptance uniforms. When the caller does not supply starting values, chain one
starts at the latent mode $\widehat{\boldsymbol s}$. For chains
$m=2,\ldots,M$,
SPmodels draws

\[
\boldsymbol s_0^{(m)}=\widehat{\boldsymbol s}
+c_0\boldsymbol L\boldsymbol z_m,
\qquad \boldsymbol z_m\sim N(\boldsymbol0,\boldsymbol I),          \tag{N-29}
\]

where $\boldsymbol L\boldsymbol L^{\mathsf T}=\boldsymbol P$ and
$c_0$ is `initial_scale`. Thus $c_0=0$ starts every chain at the mode,
whereas a positive value deliberately tests whether chains started in
different posterior regions mix towards the same distribution. A supplied
matrix instead defines every initial state explicitly.

The returned three-dimensional array keeps retained draw, latent component and
chain as separate dimensions. A pooled matrix is also returned for backwards
compatibility, but between-chain diagnostics always use the chain array.

## Rank-normalized split R-hat

For one scalar component, split each of the $M$ chains into its first and last
halves, discarding the middle draw when the retained length is odd. This gives
$2M$ sequences of length $n$. Rank all $2Mn$ values jointly. For pooled
rank $r$, use the normal score

\[
z=\Phi^{-1}\left(\frac{r-3/8}{2Mn+1/4}\right).                  \tag{N-30}
\]

For the normalized split sequences, let $\bar z_m$ and $s_m^2$ be the
chain mean and sample variance. The average within-chain variance and scaled
between-chain variance are

\[
W=\frac{1}{2M}\sum_{m=1}^{2M}s_m^2,
\qquad
B=n\operatorname{Var}(\bar z_1,\ldots,\bar z_{2M}).             \tag{N-31}
\]

Estimate the marginal variance by

\[
\widehat V^+=\frac{n-1}{n}W+\frac{B}{n},
\qquad
\widehat R=\sqrt{\frac{\widehat V^+}{W}}.                      \tag{N-32}
\]

Rank normalization makes this comparison less sensitive to a skewed scale.
Location agreement alone can miss chains with different spreads, so SPmodels
also folds the original draws around their pooled median, rank-normalizes the
absolute deviations, and repeats equations N-31--N-32. The reported value is

\[
\widehat R_{\mathrm{reported}}
=\max(\widehat R_{\mathrm{rank}},\widehat R_{\mathrm{folded}}). \tag{N-33}
\]

A value near one is necessary but not sufficient for convergence. SPmodels
records whether every component is at most 1.01 but does not discard or alter
draws to force that result. R-hat is unavailable for a single selected chain.

## Multi-chain ESS and Monte Carlo standard error

For retained chain $m$, centre draws at their chain mean and estimate the lag
$t$ autocovariance using divisor $n$. Average it over chains to obtain
$\widehat\gamma_t$. Using $W$ and $\widehat V^+$ on the original draw
scale, estimate pooled autocorrelation by

\[
\widehat\rho_t
=1-\frac{W-\widehat\gamma_t}{\widehat V^+}.                    \tag{N-34}
\]

SPmodels groups adjacent autocorrelations as in N-27, truncates before the
first non-positive pair, and makes the retained pair sums non-increasing using
cumulative minima. With integrated autocorrelation estimate
$\widehat\tau_{\mathrm{int}}$, the multi-chain effective sample size is

\[
\widehat N_{\mathrm{eff}}
=\min\left\{Mn,\frac{Mn}{\widehat\tau_{\mathrm{int}}}\right\}. \tag{N-35}
\]

Finally, if $s$ is the standard deviation of all pooled retained draws, the
Monte Carlo standard error of their mean is

\[
\boxed{\operatorname{MCSE}(\bar s)
=\frac{s}{\sqrt{\widehat N_{\mathrm{eff}}}}.}                  \tag{N-36}
\]

MCSE describes simulation error in a posterior mean conditional on the
supplied model parameters. It is not the posterior standard deviation and not
a standard error for estimated regression or covariance parameters.

This is conditional latent sampling for supplied parameters. Sampling latent
effects does not make regression or covariance parameter estimation Bayesian,
and it is not yet the Monte Carlo maximum-likelihood estimator required by the
full non-Gaussian model.

## Input and failure conditions

The covariance must be finite, symmetric and positive definite. Binomial
successes and trials must be non-negative integers with \(y_i\leq m_i\).
Poisson outcomes must be non-negative integers. Only logit and log links are
accepted. Overflowing proposals are rejected during step halving; unusable
curvature raises an error rather than being silently modified.

## Code and validation map

| Equations | Implementation | Independent checks |
| --- | --- | --- |
| N-04, N-06 | `matern_scale_derivative()` | Exponential identities and centred finite differences |
| N-09--N-15 | `.sp_latent_components()` | Direct scores and curvature matrices |
| N-16--N-19 | `sp_latent_mode()` | Near-zero final score, inverse curvature and row permutation |
| N-20--N-24 | `sp_latent_sample()` | Seed reproducibility and moments against numerical integration |
| N-25--N-28 | `sp_sampler_diagnostics()` | Single-chain ESS bounds and requested ACF length |
| N-29 | `sp_latent_sample()` | Local seeded reproducibility, dispersed starts and retained chain identity |
| N-30--N-36 | `sp_sampler_diagnostics()` | Direct MCSE identity, multi-chain ESS bounds, chain selection and rank-normalized split R-hat |

The relevant tests are in `tests/test-nongaussian-foundation.R`. The
one-dimensional binomial posterior mean and standard deviation are checked
against numerical integration using predeclared absolute margins. Four seeded
chains test reproducibility, chain identities, dispersed initialization,
multi-chain ESS/MCSE, chain selection and a converged rank-normalized split
R-hat. External reference comparisons remain pending.
