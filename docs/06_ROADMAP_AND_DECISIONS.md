# SPmodels roadmap and decisions

## Current status

The owner instructed implementation to begin on 12 October 2026. Implemented
slices now include Matérn operations, Gaussian fitting and prediction,
binomial/Poisson Laplace fitting and point prediction, conditional latent
sampling, importance-sampling MCML correction, joint predictive draws,
projected grids and draw-based targets.
Reference comparison, consumer/course revision, advanced model structures and
publication remain pending; current evidence is in `11_IMPLEMENTATION_REGISTER.md`.

## Future phases

| Phase | Work | Completion evidence |
| --- | --- | --- |
| 1 Scope and requirements | Confirm consumers, model configurations, reference comparison setup and curriculum needs | Agreed feature/configuration inventory and consumer workflow records |
| 2 Mathematical specification and derivations | Independently derive all implemented methods and develop a complete guide accessible to MSc students, with suitable lsda examples | Reviewed step-by-step chapters, equation-to-code mappings, analytic cases and predefined comparison criteria |
| 3 Core and Gaussian workflow | Build model frames, geometry, covariance, Gaussian fitting/prediction and the owned fit class | Deterministic checks and a full Gaussian validation report |
| 4 Discrete outcomes | Implement binomial/Poisson fitting, grouping effects, integration and diagnostics | Convergence and stochastic equivalence evidence |
| 5 Full initial collection | Complete temporal/intervention configurations, utilities, target aggregation, assessment and reporting | Configuration-level evidence for all agreed initial capabilities |
| 6 Application and course adoption | Update MBGapp, lsda and affected student workflows | End-to-end checks with the comparator absent |
| 7 Distribution | Prepare CRAN submission, installation guidance and optional R-universe distribution | Package checks, fresh-install trials and reviewed release materials |
| 8 Further models | Add prioritised research and teaching models in validated increments | Model-specific specifications, results and teaching examples |

Phase order may change after the consumer audit. A partial release must state its supported scope. Do not set a delivery date before measuring representative workloads.

## Decisions

| ID | Topic | Current position |
| --- | --- | --- |
| D01 | Name | SPmodels is the owner's working name; review name uniqueness and suitability before public submission |
| D02 | Project identity | Independent spatial modelling package for research and teaching; no comparator branding or replacement framing |
| D03 | Reference validation | Preserve the complete numerical comparison requirement; use neutral project wording and a separate validation environment |
| D04 | Consumers | MBGapp, lsda and the course confirmed; exact revisions and workflow inventory remain to establish |
| D05 | Interface | Formula/data/family and standard methods, inspired by SDALGCP2; provisional function names remain editable |
| D06 | Intervention front end | Unified intervention argument versus dedicated front end remains open |
| D07 | Initial scope | Preserve the full agreed functional inventory and supported configurations |
| D08 | Numerical criteria | Set deterministic tolerances and application-specific stochastic margins before evaluation |
| D09 | Data | Owner authorised suitable lsda datasets for examples; retain synthetic fixtures for known-truth checks and record dataset versions, support and provenance |
| D10 | Discrepancies | Resolve through mathematical checks; never hide failures by changing tolerances or scope |
| D11 | Implementation | Independent runtime; choose compilation and parallelism after correctness and profiling |
| D12 | Authorship and licence | Owner decision at kickoff; accurate method citations remain required |
| D13 | Expansion | Prioritise future models using the syllabus and research needs |
| D14 | Teaching | Clear assumptions, consistent methods, practical examples and model progression |
| D15 | Distribution | Plan CRAN and easy installation, with optional R-universe and a versioned source repository |
| D16 | Full derivations | Required for all implemented statistical and numerical methods, including future additions; written for MSc students with complete intermediate steps and worked examples |
| D17 | Companion guide | Versioned accessible web chapters and a downloadable edition proposed; concise package help links to the matching derivations |

## Provisional effort and compute budget

These are planning estimates, not timings from an implementation or benchmark. Assume one sustained development effort with the owner's statistical review at milestones. They describe the initial agreed model collection; future model expansion remains an ongoing programme. The owner has now explicitly required full MSc-accessible derivations for all methods. Reassess authoring and review effort after the chapter inventory is complete; the earlier ranges are provisional and do not establish a measured budget for this expanded documentation requirement.

| Milestone | Indicative elapsed development time | Main uncertainty |
| --- | --- | --- |
| First usable model and end-to-end example | 2–4 weeks | Mathematical and interface decisions, reference agreement |
| Initial collection ready for CRAN submission | 3–6 months overall | Non-Gaussian fitting, intervention models, full validation, consumer/course adoption and documentation |
| Further model families | Prioritised after the initial collection | Scope and complexity of each addition |

CRAN review and acceptance are external to the development estimate. Calendar estimates are not continuous assistant execution time, and no reliable assistant-hour or token estimate is available yet.

Budget approximately **100–500 CPU-hours** for validation as an initial placeholder, then revise it from representative timing runs once implementation is authorised. A CPU-hour is one CPU core working for an hour. For example, an illustrative suite of 40 configurations, four seeds and two implementations, each taking 20 CPU-minutes, totals about 107 CPU-hours before debugging reruns and additional checks. These are assumed inputs, not measured case counts or runtimes.

Independent cases may run concurrently, but memory limits, serial code and convergence affect elapsed time. Do not promise completion time by simply dividing CPU-hours by advertised core count. Benchmark small and representative moderate cases before fixing the full compute budget or requesting additional hardware. Model fitting and validation compute is separate from assistant usage limits and service costs.

## First implementation session

Confirm the configurations and consumer requirements, freeze the separate validation environment and define scientific equivalence margins. Then write independent mathematical specifications and start the first justified implementation milestone.

## Restart prompt

> Begin implementation of SPmodels in `/Users/user/Documents/Manchester Job/Research/SPmodels`. Read AGENTS.md, README.md and the planning documents. Develop it as my independent spatial modelling package for research and teaching, with a glm-like interface inspired by my SDALGCP2. Preserve the full agreed initial functional scope and validate each equivalent capability against the reference implementation identified in our conversation, using a separate environment and the predefined numerical and stochastic criteria. Do not copy external implementation code, name comparator packages in project materials or frame SPmodels as a replacement package. Provide full mathematical derivations of every implemented statistical and numerical method, accessible to MSc students, with assumptions, intermediate steps, approximation details and equation-to-code links. Use suitable data from my lsda package for worked examples. Plan adoption in MBGapp, lsda and my course, future model expansion, CRAN submission and easy installation. Report actual evidence and remaining gaps; do not claim completeness before all required workflows pass.

The owner subsequently gave the instruction to start; implementation status is
tracked in `11_IMPLEMENTATION_REGISTER.md`.
