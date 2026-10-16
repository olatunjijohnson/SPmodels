# Cross-validation and predictive-draw assessment

## Purpose and separation of tasks

Cross-validation measures prediction after withholding observations and
refitting. Predictive-draw assessment compares an already generated ensemble
with known simulated truth or held-out values. These are different operations:
the first guards against training-data leakage, while the second evaluates
calibration and accuracy of a predictive distribution.

## Fixed-fold refitting

Let \(f_i\in\{1,\ldots,K\}\) identify the fold for observation \(i\). For
fold \(k\), define training and test index sets

\[
\mathcal T_k=\{i:f_i\ne k\},\qquad
\mathcal H_k=\{i:f_i=k\}.                                       \tag{V-01}
\]

SPmodels refits the complete declared model using only \(\mathcal T_k\), then
predicts observations in \(\mathcal H_k\). Covariance parameters are therefore
re-estimated when the original control requested estimation; they are not
carried over from the full-data fit. Fixed covariance controls remain fixed.

User-supplied fold identifiers are preserved. If only a fold count is given,
SPmodels randomly permutes a balanced repetition of \(1,\ldots,K\). A supplied
seed is local to fold construction. Random folds are a reproducibility
convenience, not a substitute for scientifically justified spatial blocking.

## Held-out response scale

For Gaussian and Poisson models, the held-out value is \(y_i\). For a binomial
response with \(y_i\) successes from \(m_i\) trials, response-scale error uses

\[
r_i=\frac{y_i}{m_i}-\widehat p_i.                                \tag{V-02}
\]

For other supported families,

\[
r_i=y_i-\widehat\mu_i.                                           \tag{V-03}
\]

Across \(N\) held-out predictions, root mean squared error and mean absolute
error are

\[
\operatorname{RMSE}=
\left(\frac1N\sum_{i=1}^N r_i^2\right)^{1/2},\qquad
\operatorname{MAE}=\frac1N\sum_{i=1}^N|r_i|.                    \tag{V-04}
\]

These are reported overall and separately by fold.

## Negative log score

The negative log score penalises low probability assigned to what occurred.
For a Gaussian held-out predictive mean \(\widehat\mu_i\) and standard
deviation \(\widehat\sigma_i\),

\[
S_i=-\log\phi(y_i;\widehat\mu_i,\widehat\sigma_i^2).              \tag{V-05}
\]

For a non-Gaussian held-out row, the local predictive approximation on the
link scale is

\[
\eta_i\mid Y_{\mathcal T_k}\ \dot\sim\
N(m_i^*,v_i^*).                                                \tag{V-06}
\]

For a binomial observation with \(n_i\) trials, let \(a_i\) be any fitted
intervention multiplier (one when absent). The predictive mass is

\[
q_i(y_i)=\int
\Pr\{Y_i=y_i\mid n_i,a_i\operatorname{logit}^{-1}(\eta)\}
\phi(\eta;m_i^*,v_i^*)\,d\eta.                                \tag{V-07}
\]

For Poisson data it is

\[
q_i(y_i)=\int
\Pr\{Y_i=y_i\mid \exp(\eta)\}
\phi(\eta;m_i^*,v_i^*)\,d\eta.                                \tag{V-08}
\]

SPmodels evaluates each one-dimensional integral adaptively after writing
\(\eta=m_i^*+\sqrt{v_i^*}z\), \(z\sim N(0,1)\), and reports
\(S_i=-\log q_i(y_i)\). A failed or underflowed integral is an explicit error,
not a clipped finite score. The reported mean is \(N^{-1}\sum_iS_i\); smaller
is better. Integration is conditional on fitted parameters and the local
Gaussian predictive approximation.

## Assessment from predictive draws

Suppose \(D_i^{(1)},\ldots,D_i^{(B)}\) are predictive draws for location \(i\)
and \(t_i\) is known truth. The draw mean and error are

\[
\overline D_i=\frac1B\sum_{b=1}^B D_i^{(b)},\qquad
e_i=\overline D_i-t_i.                                          \tag{V-09}
\]

Bias, RMSE and MAE aggregate these errors across locations:

\[
\operatorname{Bias}=\frac1M\sum_i e_i,\qquad
\operatorname{RMSE}=\left(\frac1M\sum_i e_i^2\right)^{1/2},\qquad
\operatorname{MAE}=\frac1M\sum_i|e_i|.                          \tag{V-10}
\]

For central probability \(1-\alpha\), empirical draw quantiles define
\([L_i,U_i]\). Coverage and mean width are

\[
\operatorname{Coverage}=\frac1M\sum_i\mathbb 1(L_i\le t_i\le U_i),
\qquad
\operatorname{Width}=\frac1M\sum_i(U_i-L_i).                    \tag{V-11}
\]

## Empirical continuous ranked probability score

For one location, the empirical CRPS is

\[
\boxed{
\operatorname{CRPS}_i=
\frac1B\sum_{b=1}^B|D_i^{(b)}-t_i|
-\frac{1}{2B^2}\sum_{b=1}^B\sum_{c=1}^B
|D_i^{(b)}-D_i^{(c)}|.
}                                                                 \tag{V-12}
\]

It rewards both calibration and concentration; smaller is better. Directly
forming all \(B^2\) differences is unnecessary. If
\(D_{i(1)}\le\cdots\le D_{i(B)}\) are sorted draws, then

\[
\frac{1}{2B^2}\sum_{b,c}|D_i^{(b)}-D_i^{(c)}|
=\frac1{B^2}\sum_{b=1}^B(2b-B-1)D_{i(b)}.                       \tag{V-13}
\]

SPmodels uses (V-12) and reports the mean CRPS across locations.

For optional threshold \(q\), it also reports

\[
\widehat P(D_i>q)=\frac1B\sum_b\mathbb 1(D_i^{(b)}>q)            \tag{V-14}
\]

alongside \(\mathbb 1(t_i>q)\). No arbitrary classification cutoff is imposed.

## Calibration curve and finite-location bounds

One coverage value at one chosen interval level cannot show whether predictive
intervals are systematically too narrow or too wide. For each requested
nominal central coverage \(p\), SPmodels calculates the type-8 draw quantiles
at \((1-p)/2\) and \((1+p)/2\). If these are \(L_i(p)\) and \(U_i(p)\), the
empirical calibration curve is

\[
\widehat c(p)=\frac1M\sum_{i=1}^M
  \mathbb 1\{L_i(p)\le t_i\le U_i(p)\}.                         \tag{V-15}
\]

A calibrated ensemble should lie near the identity line
\(\widehat c(p)=p\). This is a diagnostic relationship, not a hypothesis test.
To expose the sampling uncertainty caused by having only \(M\) assessed
locations, let \(z\) be the \((1+\gamma)/2\) standard-normal quantile and set

\[
d=1+z^2/M,\qquad
a=\frac{\widehat c+z^2/(2M)}{d},\qquad
h=\frac{z}{d}\left\{
\frac{\widehat c(1-\widehat c)}M+\frac{z^2}{4M^2}
\right\}^{1/2}.                                                \tag{V-16}
\]

The reported Wilson interval is

\[
[\max(0,a-h),\min(1,a+h)].                                     \tag{V-17}
\]

These are pointwise finite-location bounds. They are not simultaneous bands,
and spatial dependence between locations can reduce their effective
information, so the bounds must not be read as a spatially adjusted confidence
statement.

## Display algorithms

The validation and simulation plot methods are deliberately views over stored,
unrounded results. Cross-validation displays observed against predicted values,
error against prediction, or error by fold. Predictive assessment displays the
stored mean and interval against truth, or error against truth. Parameter
recovery displays one named parameter across replicates, its stored truth and,
when requested, Wald intervals formed from the already reported standard
errors and recovery level. Plotting does not refit a model or alter any metric.

## Code and validation map

| Equations | Implementation | Independent checks |
|---|---|---|
| V-01 | `sp_cross_validate()` | One fold reproduced by a separate manual refit; supplied fold IDs preserved |
| V-02--V-08 | `sp_cross_validate()`, `sp_predictive_log_score()` | Gaussian direct density; independent binomial/Poisson Gaussian-link integration; finite fold scores |
| V-09--V-11 | `sp_simulation_assess()` | Direct draw means, zero-RMSE fixture and interval coverage |
| V-12--V-13 | `.sp_empirical_crps()` | Agreement with brute-force pairwise absolute differences |
| V-14 | `sp_simulation_assess()` | Direct exceedance frequencies |
| V-15--V-17 | `sp_calibration()` | Direct type-8 interval coverage and independent Wilson-bound reconstruction |
| Display algorithms | validation and simulation `plot()` methods | Every supported view renders from fixed synthetic stored values |

Buffered, regular-block and environmental fold constructors, Cartesian
parameter-grid generation, scenario-level parallel recovery and structured
reporting are now implemented and derived in the later chapters. Remaining
work includes further domain-specific fold constructors only when a consumer
use case requires them, and consumer-specific integration where authorized.
