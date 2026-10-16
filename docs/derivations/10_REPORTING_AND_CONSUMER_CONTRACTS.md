# Structured reporting and consumer contracts

## Purpose

Reporting functions do not introduce a new fitted probability model. They
provide stable transformations of already defined SPmodels results so that an
application, lesson or report need not inspect internal object fields.

## Model table

For mean parameter \(\hat\beta_j\) with reported covariance matrix
\(V_{\hat\beta}\), the table records

\[
\operatorname{SE}(\hat\beta_j)=\sqrt{(V_{\hat\beta})_{jj}},\qquad
z_j=\frac{\hat\beta_j}{\operatorname{SE}(\hat\beta_j)},\qquad
p_j=2\Phi(-|z_j|).
\]

These are large-sample Wald summaries, not exact finite-sample tests. Covariance
parameters retain their fitted standard errors but receive no automatic Wald
p-value because boundary constraints make that interpretation unreliable.
Every row states whether a value was estimated or fixed. The Matérn scale row
also carries the fitted coordinate distance unit.

Implementation: `R/reporting.R`, function `sp_model_table()`.

## Prediction table

For stored joint draws \(t_j^{(1)},\ldots,t_j^{(B)}\) at target \(j\), the
reported quantities are

\[
\bar t_j=B^{-1}\sum_b t_j^{(b)},\qquad
s_j=\left\{(B-1)^{-1}\sum_b(t_j^{(b)}-\bar t_j)^2\right\}^{1/2},
\]

the type-8 empirical quantiles used throughout SPmodels, and, for threshold
\(c\),

\[
\widehat{P}(T_j>c)=B^{-1}\sum_b I(t_j^{(b)}>c).
\]

Grouped targets are first calculated inside every joint draw using the weights
defined in derivation 06. This preserves spatial dependence and is different
from aggregating marginal means after prediction.

Implementation: `R/reporting.R`, function `sp_prediction_table()`, delegating
the calculations to `sp_predict_target()`.

## Validation table

`sp_validation_table()` changes representation only. It reshapes previously
calculated validation metrics into the common tuple

\[
(\text{assessment},\text{scope},\text{id},\text{metric},\text{value}).
\]

The definitions of RMSE, MAE, log score, coverage, CRPS, recovery bias and
recovery RMSE remain those in derivations 07 and 09. No rounding occurs in the
machine-readable table.

## Presentation formatting

`sp_format_table()` first obtains one of the three stable tables above, or
accepts an already constructed data frame. Let \(v\) be a finite numerical
value and \(d\) the requested number of decimal places. Its displayed value is
the fixed-decimal representation

\[
F_d(v)=\operatorname{format}(v,\text{fixed},d).                 \tag{R-01}
\]

For a p-value column, values below \(10^{-d_p}\) are displayed as
`<0.00...1` at precision \(d_p\), rather than as zero. Missing numerical and
character values use the declared missing-value marker. This operation changes
representation only: the original unrounded data frame is retained verbatim in
the `raw_table` attribute of the formatted view.

## Assembled report

`sp_report()` constructs its statistical sections only by calling
`sp_model_table()`, `sp_prediction_table()` and `sp_validation_table()`. The
result stores two parallel collections:

1. the unrounded, machine-readable source tables; and
2. formatted character views produced from those same tables.

The report printer adds a title, basic analysis metadata and section labels. It
does not refit, predict, resample, recalculate metrics or replace raw numbers
with rounded ones. The resulting object can therefore support a readable
console report while downstream applications continue to consume exact values.

## Code and validation map

| Operation | Implementation | Independent checks |
|---|---|---|
| Model estimates | `sp_model_table()` | Coefficients, covariance values, units and status reconstructed directly |
| Prediction summaries | `sp_prediction_table()` | Draw means, SD, quantiles and exceedance frequencies reconstructed directly |
| Validation schema | `sp_validation_table()` | Overall, fold and parameter scopes checked against source summaries |
| R-01 and p-value threshold | `sp_format_table()` | Exact fixed-decimal strings, threshold label and missing marker |
| Report assembly | `sp_report()` | Every raw section is identical to its stable public table; printed sections are present |

## Tests and limitations

`tests/test-reporting.R` checks values against direct calculations, validates
the fixed schemas and exercises Gaussian prediction and cross-validation.
`tests/test-report-presentation.R` additionally checks formatting and assembled
reports. The current report is a dependency-free console/report object, not a
typeset journal template, HTML site or office-document exporter. Formatting
remains a view over the unrounded tables.
