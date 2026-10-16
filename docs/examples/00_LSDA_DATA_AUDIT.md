# lsda teaching-data suitability audit

## Purpose

This record separates datasets that can support the current SPmodels point
geostatistical interface from datasets whose support or metadata require a
different model or additional provenance work. It prevents a convenient but
scientifically invalid example from being presented as evidence of coverage.
The audit is read-only: no lsda source or data file was changed.

## Inspected source

- local package source: `lsda` 0.2.0;
- package licence: CC BY-NC-SA 4.0;
- inspection date: 15 October 2026;
- data remain external to SPmodels and are not runtime dependencies.

## Decisions

| Object | Response/support | Current decision | Reason |
|---|---|---|---|
| `galicia` | Continuous response at 132 points | Suitable; checked Gaussian example exists | Coordinates and response/predictor fields are explicit; points can be projected before fitting |
| `loaloa` | Positives among examined people at 197 villages | Suitable; checked binomial example exists | The response has a direct binomial denominator and documented point coordinates |
| `liberia` | Continuous cluster summary and an incidence *rate* at 267 points | Not a count example | A rate without its event count and exposure cannot be reconstructed as a Poisson response |
| `scotland_lip` | Cases and expected counts on 56 polygons | Defer until an areal count model exists | Replacing polygons by centroids would change the support and imply a point process that was not observed |
| `nw_obesity` | Obese among sampled children on 33 polygons | Defer until an areal/binomial model exists | It is polygon-supported; centroid point fitting would discard the stated support |
| `gm_no2` | Continuous pollution measurements at 25 points | Suitable for a second Gaussian example if needed | It does not fill the remaining count/intervention teaching gaps |
| `galicia_pred` plus `galicia_grid` | Prediction covariates and grid coordinates | Alignment unresolved | The inspected metadata do not state an exact coordinate transformation or identifier join between the two 1,508-row objects |

No inspected object contains an auditable intervention history with unit,
survey time, rollout time and delivery status. An intervention example should
therefore use a documented synthetic fixture or await a suitable independently
provenanced dataset; a cross-sectional field must not be relabelled as history.

## Consequences for the completion checklist

The Gaussian and binomial real-data teaching requirements have suitable
inputs. The remaining count, grid-target and intervention examples are still
required, but their absence is now an explicit data/model-support boundary.
They do not justify copying data, inventing exposure denominators, assuming a
row-wise grid join, or treating areal counts as point observations.

