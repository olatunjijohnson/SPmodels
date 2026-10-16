# SPmodels feature inventory

This inventory specifies the complete initial model collection and its supporting operations. Some capabilities are now implemented and tested; evidence and remaining boundaries are recorded in [the implementation register](11_IMPLEMENTATION_REGISTER.md). Future additions are tracked in [Teaching and future models](08_TEACHING_AND_FUTURE_MODELS.md). Preserve the complete agreed initial scope when deciding implementation priorities.

## Modelling and supporting capabilities

These are functional requirements, not a promise to reproduce an external package's function names or object classes. Each capability needs its own SPmodels contract and comparison cases.

| ID | SPmodels capability | Validation focus |
| --- | --- | --- |
| F01 | Latent sampling interface | Conditional moments, dependence, MCSE, diagnostics |
| F02 | Spatial validation | Fixed folds, refitting policy, scores and calibration |
| F03 | Simulation assessment | Target definitions, errors and classification summaries |
| F04 | Sampler diagnostics | Selected chains, trace/ACF values and ESS conventions |
| F05 | Coordinate identity | Duplicates, row order and unique-location mapping |
| F06 | Geometry helper | Hull geometry, CRS and degenerate inputs |
| F07 | Grid construction | Resolution, point ordering, boundaries and CRS |
| F08 | Intervention fitting | History, survey/GP times, constraints and impact parameters |
| F09 | Distance summaries | Distance units, summary definitions and duplicates |
| F10 | Main fitting interface | Families, likelihood, estimates, uncertainty and controls |
| F11 | Model simulation | Outcome moments, spatial dependence and returned schema |
| F12 | Covariance specification | Argument meaning, smoothness and variance-component mapping |
| F13 | Covariance derivative | Scale convention, numerical values and finite differences |
| F14 | Covariance curvature | Scale convention, numerical values and finite differences |
| F15 | Matérn correlation | Zero distance, smoothness, limits and matrix properties |
| F16 | Latent-mode calculation | Objective, mode, curvature and solver convergence |
| F17 | Calibration display | Underlying calibration values, bounds and labels |
| F18 | Intervention display | Time alignment, effect values and interval meanings |
| F19 | Variogram display | Bins, distances, fitted values and units |
| F20 | Validation display | Model/metric selection, score values and labels |
| F21 | Simulation display | Selected replicate, spatial values and scale |
| F22 | Spatial prediction | Marginal/joint moments, draws and included components |
| F23 | Grid targets | Transformations, thresholds, uncertainty and component flags |
| F24 | Polygon targets | Joint aggregation, weights, region IDs and empty regions |
| F25 | CRS helper | Zone/hemisphere boundaries and invalid inputs |
| F26 | Grouping specification | Levels, multiple effects and variance parameter meaning |
| F27 | Empirical variogram | Residual definition, lag bins and pair counts |
| F28 | Simulation controls | Resolved defaults, validation and effective iteration counts |
| F29 | Surface simulation | Spatial process, response generation and replicates |
| F30 | Table output | Estimates, uncertainty, labels and generated table content |
| F31 | Prediction update | New design, offsets, row alignment and preserved quantities |

## Result methods and reporting

| ID | Operation | Validation focus |
| --- | --- | --- |
| M01 | Coefficient extraction | Names, transformations, ordering and values |
| M02 | Grid target plot | Values, geometry, scale, labels and intervals |
| M03 | Polygon target plot | Regional values, geometry, scale and intervals |
| M04 | Printed fit summary | Numbers, rounding, labels, units and readability |
| M05 | Printed simulation summary | Aggregate values and uncertainty |
| M06 | Printed cross-validation summary | Metrics, model labels and uncertainty |
| M07 | Fit summary | Estimates, uncertainty, diagnostics and warnings |
| M08 | Simulation summary | Recovery, prediction and classification results |
| M09 | Cross-validation summary | Scores, calibration and fold information |

## Data and example workflows

Account for all data workflows needed to exercise the complete agreed scope and the owner-confirmed consumers. Inventory schemas, units, spatial support, CRS, missingness, response encodings and target outputs. Cover continuous, prevalence, count, intervention/time, covariate-grid and regional aggregation examples.

Use small independently generated fixtures by default. Record original sources and permitted uses before distributing real data. Never imply that externally authored datasets belong to SPmodels or silently replace real data with synthetic observations under an existing name.

## Additional simulation requirement

| ID | Capability | Work before implementation |
| --- | --- | --- |
| W01 | Intervention model simulation | Specify generation from fitted or supplied parameters, supported histories and times, output schema and reference comparability |

## Expand each capability at kickoff

For every ID, add its proposed API, resolved defaults, supported configurations, input/output schema, failure conditions, owning module, mathematical source, derivation chapter and equation IDs, MSc prerequisite notes, worked-example dataset, consumer workflows, fixture IDs, numerical criteria and status. Statistical and numerical methods require full derivations; presentation and data-handling utilities require an explicit algorithm and rationale rather than artificial mathematics. Use the validation protocol for the comparison method.

Priorities follow MBGapp, lsda and course needs. They do not remove less frequently used capabilities. Low-level numerical operations remain in scope where they support reproducible research or consumers.

## Coverage completion rule

A function-name checklist is insufficient. Cover supported families, links, grouping and temporal structures, fixed/estimated components, prediction targets, geometry policies and controls. Reconcile additional requirements discovered in the consumer audit before declaring completeness. The omitted comparator branding does not narrow the agreed numerical validation requirement.
