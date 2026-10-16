# SPmodels teaching and future models

## Purpose

Develop a model collection that serves research and teaching through a consistent R interface. The initial collection implements the agreed geostatistical capabilities independently. Later models should be selected for scientific and teaching needs, with their own statistical specifications and validation.

This is a planning document. No additional model implementation or course revision has been authorised by creating this roadmap.

## Teaching requirements

- Start each model with the question it answers, the observational unit and the spatial or temporal support of the data.
- Present its likelihood, mean structure, dependence assumptions and identifiable parameters before interpreting a fit.
- Provide full mathematical derivations for every implemented statistical and numerical method, written for MSc students and following [Derivations and teaching data](10_DERIVATIONS_AND_TEACHING_DATA.md). Explain intermediate steps and approximation error, not only the final formulas.
- Use familiar formulas and standard methods, with concise defaults and accessible advanced controls.
- Use suitable data from the owner's lsda package for worked examples, with exact versions, variables, provenance, units and coordinate systems recorded. Keep small synthetic fixtures for simulation exercises, hand calculations and known-truth checks.
- Explain coefficients, covariance parameters and uncertainty in terms students can connect to the model equations.
- Include diagnostics, model limitations and prediction exercises. Use a comparable example across related models where that comparison is scientifically meaningful.
- Keep classroom examples practical on student computers. Record seeds and runtime expectations, and use clearly labelled saved results only for optional expensive demonstrations.
- Keep tutorial solutions and regression checks tied to a package version so course material remains reproducible.

## Candidate model progression

These are candidate areas for later prioritisation, not a promised release list. Review the actual course syllabus and research needs at kickoff.

| Stage or area | Teaching purpose | Scope position |
| --- | --- | --- |
| Independent linear and generalised linear baselines | Establish regression, links, offsets and uncertainty before adding dependence | Baseline comparisons and teaching examples |
| Gaussian geostatistical models | Connect covariance, likelihood and spatial prediction | Initial model collection |
| Binomial and Poisson spatial models | Introduce latent dependence for prevalence and counts | Initial model collection |
| Supported temporal and intervention models | Explain temporal structure and intervention effects | Initial reference-supported configurations |
| Areal models such as CAR or SAR | Distinguish neighbourhood-based dependence from distance-based covariance | Future candidate; specification and priorities unresolved |
| Point-process and aggregated-count models | Distinguish event intensity, point support and aggregation | Future candidate; review the relationship with the owner's SDALGCP2 |
| Additional spatio-temporal models | Explore temporal dependence, forecasting and different observation supports | Future candidate beyond the agreed reference |
| Further response families and inference methods | Address research and curriculum needs as they emerge | Add individually after specification and validation |

## Requirements for each additional model

Maintain a model record with its learning objective, scientific use, response and data support, likelihood/link, dependence structure, inference method, parameterisation, input validation, prediction targets and uncertainty interpretation. Link it to a mathematical derivation, published method where applicable, independent numerical checks, an appropriate reference implementation where available and a small worked example.

Add a model only when its fit, simulation, diagnostics and prediction behaviour are defined and validated to the stated release scope. Reuse shared infrastructure where the mathematical contract is the same. Do not force different models into a single opaque argument list or imply that matching function names makes their outputs directly comparable.

Reference comparisons remain mandatory for equivalent initial capabilities. For a model outside that reference, use analytic results, limiting cases, simulation experiments and suitable independent implementations. State unsupported operations explicitly.

## Relationship with existing packages

MBGapp and lsda remain the owner-identified consumers. SPmodels owns statistical modelling and its model documentation; lsda can continue to organise course-specific materials, datasets and teaching applications. Confirm those boundaries when auditing the actual repositories. SDALGCP2 is the interface inspiration, and its possible future relationship with the wider model collection is a separate design decision.

The owner confirmed that the course should eventually use SPmodels for the affected modelling workflows. Record affected lessons, datasets, teaching applications and setup instructions, then migrate only after the student workflows are validated without the comparator installed. Keep the modelling package suitable for CRAN and place extensive course notes on a separate course site or repository.

No current teaching files or dependent package code are changed by this planning update.
