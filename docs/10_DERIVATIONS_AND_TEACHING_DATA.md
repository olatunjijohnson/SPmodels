# SPmodels derivations and teaching data

## Requirement and audience

Every statistical and numerical method implemented in SPmodels must have a full, independently written mathematical derivation accessible to MSc students. This requirement covers the initial collection and all future additions. A method is not ready for release with only an external citation, algorithm name or final formula.

Assume introductory probability, regression, calculus and linear algebra. State the prerequisites for each chapter and provide short refreshers for unfamiliar material. Use consistent notation, define dimensions and units, show intermediate steps and connect equations to their statistical meaning. Keep optional advanced proofs in appendices while retaining all steps needed to follow and use the implemented method in the main chapter.

This document specifies future deliverables. It does not begin derivation writing, data extraction or model implementation.

## Required structure for each method

1. **Question and intuition.** Explain the scientific problem, observational unit, data support and what the method estimates or predicts.
2. **Prerequisites and notation.** Define every symbol, vector/matrix dimension, parameter constraint and unit; link to prerequisite material.
3. **Model and assumptions.** Give the data-generating model, conditioning structure, independence/dependence assumptions and identifiability conditions.
4. **Full derivation.** Start from the relevant probability model or objective and show each substantive algebraic step. Derive the implemented likelihood, estimating equations, gradients, Hessians, variance expressions and transformations where applicable.
5. **Computation.** Derive the algorithm actually used, including parameter transformations, constraints, initialisation, updates, stopping rules and numerical stability choices. Explain how equations become implementation routines.
6. **Approximation and simulation.** Distinguish exact identities from approximations. Where closed forms are unavailable, derive the approximation or numerical integration/sampling method and explain its assumptions, accuracy, convergence and diagnostics. Do not invent a closed-form solution.
7. **Prediction and uncertainty.** Derive the supported targets and uncertainty, distinguishing latent from observed outcomes, conditional from marginal quantities, and estimation uncertainty from Monte Carlo error.
8. **Worked examples.** Include a small hand-checkable or synthetic example and, where suitable, a realistic lsda example with all preparation steps explained.
9. **Checks and limitations.** Show a special or limiting case, independent numerical checks, likely failure modes and interpretation limits.
10. **Practice and references.** Include short exercises with solutions or guided answers, and accurate citations to mathematical/statistical sources.

## Coverage

The guide must cover likelihood construction, estimation and optimisation, covariance functions and derivatives, random effects, numerical integration and sampling, model simulation, prediction and aggregation, intervention/time components, uncertainty, diagnostics, validation scores and model-comparison procedures that are actually implemented.

For plotting, formatting and data-handling utilities, explain the algorithm, transformation and rationale; no artificial mathematical derivation is required when no statistical method is involved. Reuse shared derivations through explicit cross-references, and show every model-specific adaptation rather than repeating generic text.

## Equations and implementation must agree

Maintain a coverage record linking each method and supported configuration to its chapter, numbered equations, parameter conventions, implementation routines, test cases and worked-example dataset. Check algebra, dimensions, signs, limiting cases and numerical results. Review accessibility separately: an MSc reader should be able to trace the derivation without filling in unexplained key steps.

Version the guide alongside the package. Each released method needs a reviewed chapter and checked example outputs. Changes to algorithms, parameterisation or uncertainty must update the derivation and its validation links in the same release.

## Candidate lsda examples

The owner has authorised using suitable lsda data. These candidates were identified from local lsda dataset documentation on 12 October 2026; the proposed applications below remain choices for implementation.

| Dataset | Documented content | Proposed use |
| --- | --- | --- |
| galicia | Lead measurements at sampling stations with a co-located covariate | Gaussian geostatistics, variograms, covariance estimation and spatial prediction; justify any transformation |
| gm_no2 | Nitrogen dioxide at Greater Manchester monitoring stations | A compact spatial example linking equations, model fitting and maps |
| loaloa | Village survey positives and examined totals with environmental covariates | Binomial spatial regression, denominators, prevalence prediction and uncertainty |
| scotland_lip | Observed and expected counts with county geometries | Future areal modelling and offsets, once that model family is implemented; do not treat polygon support as point observations |

Other lsda datasets can be selected when they fit the model and teaching objective. Their presence in lsda does not by itself extend the implementation scope.

For every selected dataset, record the lsda version or immutable revision, object name, original source, variables, units, CRS, observation support, missing-data handling, transformations and intended teaching question. Preserve appropriate attribution and data terms. Use explicit package-qualified loading in examples and do not silently alter the original objects.

Prefer using lsda data in the companion examples without making SPmodels depend on lsda at runtime, since lsda may itself adopt SPmodels. Decide at release preparation whether selected examples load lsda separately or use appropriately documented data distributed with the teaching materials. No datasets are copied or fitted during this planning update.

## Accessible delivery

Plan a versioned companion guide with searchable web chapters and a downloadable edition for offline study. Package help and concise vignettes should link directly to the relevant method chapter. Keep the complete derivations easy to find from the package documentation and course materials while keeping the package distribution focused on reusable modelling functionality.
