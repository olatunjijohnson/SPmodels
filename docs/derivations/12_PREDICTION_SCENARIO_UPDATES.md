# Prediction scenario updates

## Purpose and boundary

A prediction scenario asks what the fitted model predicts if selected covariate
or exposure values are changed while the fitted parameters, observation data,
locations and every unmentioned prediction-data value remain fixed. It is not
a refit, a causal-effect estimator or an intervention model. The calculation
is useful for transparent sensitivity exercises and for constructing matched
baseline/scenario maps.

`sp_prediction_update()` implements this contract. It accepts only variables
used on the right-hand side of the fitted formula, including variables nested
inside `offset()`. It rejects separately supplied fitting offsets because a
safe prediction-data mapping for those values does not exist.

## Row mapping and replacement

Let the prediction data contain rows \(i=1,\ldots,m\), with unique stored row
identifiers \(r_i\). The caller selects a non-empty index set

\[
\mathcal I\subseteq\{1,\ldots,m\}.                              \tag{P-01}
\]

Logical selection must have length \(m\). Integer selection must contain valid,
unique indices. Character selection is matched exactly to the row identifiers;
partial matching is never used. If replacement values are named, their names
must be exactly the selected identifier set and the values are reordered to
the selected data rows before assignment. This prevents a changed value from
being silently attached to the wrong prediction location.

For changed variable \(q\), scalar input is repeated across \(\mathcal I\), or
one replacement \(q_i^*\) is supplied for every selected row. The scenario is

\[
q_i^{\mathrm{new}}=
\begin{cases}
q_i^*, & i\in\mathcal I,\\
q_i,   & i\notin\mathcal I.
\end{cases}                                                     \tag{P-02}
\]

Existing factor levels and their ordering are retained. New or missing factor
levels fail explicitly rather than becoming a missing design value.

## Rebuilding the fitted design and offset

Let \(T\) denote the response-free fitted terms object, including its factor
levels and contrast definitions. SPmodels independently constructs

\[
X=T(\mathcal D),\qquad X^*=T(\mathcal D^*),                    \tag{P-03}
\]

from the baseline and scenario data. Formula transformations are therefore
re-evaluated rather than edited inside a design matrix. Likewise, an
in-formula offset expression produces \(o\) and \(o^*\) from the two model
frames. For example, changing exposure \(E_i\) in
`offset(log(exposure))` changes the offset from \(\log E_i\) to
\(\log E_i^*\) exactly once.

The ordinary link-scale baseline and scenario means have the common fitted
parameter vector and common conditional spatial contribution \(c_i\):

\[
\eta_i=x_i^{\mathsf T}\widehat\beta+o_i+c_i,\qquad
\eta_i^*={x_i^*}^{\mathsf T}\widehat\beta+o_i^*+c_i.           \tag{P-04}
\]

The spatial contribution is common because coordinates are immutable in this
operation. SPmodels calls the package prediction method for both designs, so
Gaussian response predictions use \(\text{P-04}\), while binomial and Poisson response
predictions apply the already derived uncertainty integration rather than only
transforming the link mean. Requested standard errors are recalculated for
each design; they are not copied from the baseline.

The reported scenario contrast is descriptive:

\[
\Delta_i=\widehat y_i^*-\widehat y_i.                           \tag{P-05}
\]

It is conditional on the same fitted model and must not be interpreted as a
causal effect without a separate identification argument and appropriate model.

## Preservation contract

Before prediction, the implementation verifies all of the following:

1. row identifiers and row order are identical;
2. every unmentioned data column is exactly unchanged;
3. coordinate matrices are exactly unchanged;
4. `sf` geometry is exactly unchanged when present; and
5. the fitted object is never modified or refitted.

Coordinate variables cannot be included in `changes`. A location-change
exercise is a new prediction design, not a covariate scenario, and should be
supplied explicitly as separate prediction data.

## Code and validation map

| Contract | Implementation | Independent checks |
|---|---|---|
| P-01--P-02 | `.sp_update_rows()`, `.sp_update_values()` | Numeric, logical and exact named-row mapping; reordered named values; invalid and empty selections |
| P-03 | `sp_prediction_update()` and `.sp_prediction_design()` | Direct manual reconstruction of updated data, design and log-exposure offset; factor-level preservation |
| P-04--P-05 | baseline and scenario calls to `predict()` | Direct equality with two independent prediction calls, including standard errors and differences |
| Preservation contract | explicit returned audit | Unchanged extra columns, coordinates, row order, `sf` geometry and fitted parameter values |

The implementation and tests are in `R/prediction-update.R` and
`tests/test-prediction-update.R`.
