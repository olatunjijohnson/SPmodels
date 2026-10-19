# SPmodels project charter

## Purpose

Create an independently implemented collection of spatial models for research and teaching. Students and researchers should be able to move between models using familiar formulas, consistent result methods and clear explanations of assumptions and uncertainty.

MBGapp, lsda and the course are the first adoption targets. SPmodels should supply their required modelling functionality while supporting a broader collection of models over time.

## Requirements

1. Implement statistical methods independently; do not copy or translate external package code.
2. Provide a concise formula, data and family interface with standard R methods, inspired by the owner's SDALGCP2.
3. Preserve the complete agreed initial functional scope, including fitting, simulation, prediction, validation, spatial utilities, reporting and example-data workflows.
4. Check every equivalent initial capability against the selected reference implementation using identical inputs, aligned parameterisations and explicit numerical criteria. Keep comparator identity out of project materials.
5. Validate future models against independent derivations, analytic or limiting cases, simulation studies and appropriate numerical references where available.
6. Make SPmodels, downstream applications and teaching workflows usable without the comparator package installed, directly or indirectly.
7. Supply reproducible small examples, explicit model assumptions, understandable summaries and a progression from simpler to richer models.
8. Define each model's response and data support, likelihood, dependence structure, inference method, prediction targets and uncertainty.
9. Prepare for CRAN submission and easy installation, with optional R-universe distribution. Keep extensive course materials separate from the modelling package.
10. Provide full mathematical derivations for every implemented statistical and numerical method, accessible to MSc students. Include prerequisites, assumptions, intermediate steps, exact versus approximate results, computation and interpretation; citations and final formulas alone are insufficient.
11. Use appropriate lsda datasets for worked examples, alongside synthetic data for validation. Record exact dataset versions, variables, transformations, units, spatial support and original sources.

## Initial model scope

The initial collection includes Gaussian, binomial and Poisson geostatistical regression, Matérn dependence, unstructured random effects, spatial prediction and uncertainty. It also includes decay-adjusted spatio-temporal intervention models, simulation, model assessment, mapping and geometry utilities. Confirm detailed configurations before implementation.

Specify Gaussian maximum likelihood and Monte Carlo maximum likelihood for the initial non-Gaussian fitting workflows where those match the agreed model specification. Sampling latent effects does not by itself make parameter estimation fully Bayesian. Define the estimand and uncertainty explicitly.

A glm-like interface does not remove spatial dependence from the model. Additional response families and model classes belong in the future-model roadmap, with their own specifications and validation.

## Architecture

Separate model-frame construction, geometry, covariance calculations, fitting, simulation, prediction, diagnostics and result methods. Expose reusable low-level numerical capabilities where consumers need them. Keep optional plotting dependencies separate where practical.

Begin with a transparent implementation and small analytic checks. Choose compiled kernels or parallelism after correctness and profiling establish a need. Set performance targets using actual research and student workloads.

## Initial programme acceptance criteria

- Every agreed capability and supported configuration has an implementation, documentation and linked validation evidence.
- Deterministic, optimisation, stochastic and end-to-end comparisons meet criteria defined before evaluation.
- MBGapp, lsda and the affected course workflows pass their adoption checks in an environment without the comparator.
- The installation and runtime dependency graph is independent of the comparator.
- Documentation explains assumptions, parameterisation, targets, uncertainty, reproducibility and use in applications.
- Every released model has a small checked teaching example and a clear statement of data support, using suitable lsda data where appropriate.
- Every implemented statistical or numerical method has a complete, reviewed MSc-accessible derivation, with equations linked to implementation routines and validation cases. The companion guide is versioned with the corresponding package release.
- Supported-platform package checks and clean installation trials pass; benchmarks report accuracy and resources together.

An earlier release may implement an explicitly stated subset. It must not claim complete coverage until the full agreed scope passes validation.

## Method provenance

Maintain a derivation and provenance record for each mathematical component, including published methods, independently written equations, implementation revision and validation cases. Use synthetic fixtures by default and verify the provenance and permitted use of any external data. Do not imply ownership of another author's work.

The owner selected the MIT License on 19 October 2026. Continue to maintain accurate authorship, copyright and method citations.
