# SPmodels application and course adoption

## Initial consumers

| Consumer | Owner confirmed role | Work to plan |
| --- | --- | --- |
| MBGapp | Application adoption target | Audit fitting, predictions, maps and reports; update to SPmodels after validation |
| lsda | Teaching package adoption target | Audit package, application, data and tutorial workflows |
| Course | Eventual use of SPmodels for the affected models | Update examples, learning materials and installation instructions after student workflows pass |
| SDALGCP2 | Preferred interface style | Review conventions; do not assume a code migration is required |

The known local lsda directory is `/Users/user/Documents/Manchester Job/Teaching/Courses/Longitudinal and Spatial Data Analysis/LSDA_2026/pkg/lsda`. Repository references are `olatunjijohnson/MBGapp` and `olatunjijohnson/lsda`; confirm branches and immutable revisions at kickoff.

## Consumer audit

Record each consumer's version, repository, existing modelling calls, arguments, data support, object fields, class expectations, diagnostics, tests and example workflows. Include application and course files outside the package R directory. Inspect explicit imports, dynamic loading and indirect dependencies without reproducing external package names in project documentation.

Map scientific capabilities to SPmodels feature IDs and its owned API. Keep an explicit account of changes in parameterisation, default behaviour, outputs and uncertainty. Prefer updating consumer code to the new interface over exporting a compatibility interface from SPmodels.

## Adoption record

| Field | Value to establish |
| --- | --- |
| Consumer and immutable revision | Exact code version |
| Workflow and call site | Scientific task and affected files |
| Feature and configuration IDs | Required SPmodels capabilities |
| Input and output contract | Shapes, data support, units and statistical meanings |
| Proposed new call | SPmodels API and required data preparation |
| Validation case and criteria | Numerical tolerances or stochastic margins set in advance |
| Installation check | Clean environment without the comparator |
| Release and rollback | New release and prior stable consumer revision |
| Status | Planned; no migration performed |

## Sequence

1. Freeze the consumer environment and record baseline scientific outputs.
2. Prioritise its required capabilities without shrinking the agreed overall scope.
3. Implement and validate each capability independently before changing consumer calls.
4. Make the smallest consumer change that uses the SPmodels API.
5. Compare estimates, uncertainty, predictions, threshold decisions, plots, downloads and error handling using unchanged inputs.
6. Test installation, loading and complete workflows without the comparator present.
7. Record release notes, student-facing API changes and a rollback version.

## MBGapp scenarios

Plan Gaussian, prevalence and count workflows covering data upload, formula construction, fitting, grid prediction, uncertainty and threshold maps, and generated reports. Establish the actual estimator and interval meaning from the statistical specification; do not rely on a broad Bayesian or frequentist label.

## Course scenarios

For each affected lesson, record the learning objective, model, data, current student output and setup instructions. Prepare an equivalent SPmodels example that preserves the statistical meaning and explains interface changes. Check related lsda tutorials and teaching applications, including undeclared dependencies and files outside the package directory.

Validate the student workflow in a clean environment before switching course materials. Keep lecture notes, exercises and large assets in the course repository/site; provide concise reusable model examples in SPmodels. Link affected lessons to the full MSc-level method derivations, and use suitable existing lsda datasets to support continuity between teaching examples and package workflows. Avoid a circular runtime dependency when lsda adopts SPmodels.

No consumer code or course files have been changed by this preparation.
