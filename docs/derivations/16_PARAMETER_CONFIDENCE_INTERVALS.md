# Transformed local parameter confidence intervals

## Purpose and interpretation

`confint()` summarizes local fitted-parameter uncertainty without allowing a
Wald interval to cross a parameter's mathematical support. These are
large-sample intervals derived from observed-Hessian standard errors, not exact
confidence sets, likelihood-profile intervals, boundary corrections or Monte
Carlo error intervals.

Let \(\widehat\theta\) have standard error \(s_\theta\), and let
\(z=z_{(1+c)/2}\) for confidence level \(c\).

## Unconstrained coefficients

Mean coefficients have real-line support, so their interval is

\[
\widehat\theta\pm z s_\theta.                                 \tag{H-01}
\]

## Positive parameters

Variances, spatial and temporal scales, random-effect variances and decay are
positive. With \(g(\theta)=\log\theta\), the delta method gives

\[
\operatorname{se}\{g(\widehat\theta)\}
\approx s_\theta/\widehat\theta.                              \tag{H-02}
\]

Transforming a symmetric log-scale interval back gives

\[
\left[\exp\{\log\widehat\theta-zs_\theta/\widehat\theta\},
\exp\{\log\widehat\theta+zs_\theta/\widehat\theta\}\right]. \tag{H-03}
\]

Both limits are positive; ordinary-scale asymmetry is intentional.

## Intervention impact constrained to zero--one

For \(0<\alpha<1\), use \(g(\alpha)=\operatorname{logit}(\alpha)\). Since

\[
g'(\alpha)=\{\alpha(1-\alpha)\}^{-1},                         \tag{H-04}
\]

the transformed standard error is

\[
s_g=s_\alpha/\{\widehat\alpha(1-\widehat\alpha)\}.            \tag{H-05}
\]

The ordinary-scale interval is

\[
\left[\operatorname{logit}^{-1}\{g(\widehat\alpha)-zs_g\},
\operatorname{logit}^{-1}\{g(\widehat\alpha)+zs_g\}\right], \tag{H-06}
\]

which lies strictly between zero and one.

## Selection, missing curvature and checks

Only estimated parameters appear. Qualified identifiers such as `mean::x`,
`covariance::scale`, `intervention::alpha` and `temporal::scale` avoid
ambiguous names. Selection is exact.

If an observed Hessian cannot provide a finite standard error, both limits
remain missing. SPmodels does not silently substitute curvature. Weak
identification, boundary estimates and small samples can make the local normal
approximation unreliable.

| Equations | Implementation | Independent checks |
|---|---|---|
| H-01 | `confint.spmodels_fit()` | Direct coefficient interval |
| H-02--H-03 | positive branch | Direct nugget, decay and temporal-scale transforms; positivity |
| H-04--H-06 | impact branch | Direct logit reconstruction and zero--one support |
| selection/failure policy | qualified rows | Fixed omission, ordering and invalid-selection errors |

The checks are in `tests/test-confint.R`.
