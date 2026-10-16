# Instructions for SPmodels

## Current scope

This directory is a planning pack, not an R package. The owner requested the folder and kickoff documents only. Do not create package code, a package skeleton, numerical fixtures, installations, benchmarks, migrations, a repository or releases until the owner explicitly asks to begin implementation. Documentation corrections remain within scope. A later explicit instruction to start is sufficient; do not request redundant permission.

## Project identity and independence

Present SPmodels as the owner's independent spatial modelling package for research and teaching. Do not name external comparator packages, reproduce their branding or API catalogues, or frame this work as a clone, replacement or compatibility version in package or project materials. The owner does not claim ownership of external packages.

Implement methods independently from mathematical specifications and appropriately cited research. Do not copy, translate, adapt, decompile or use another package's implementation as a template. Do not copy its tests or documentation prose. Retain accurate citations to statistical methods and any attribution that is actually required by material used.

Preserve the agreed requirement to validate all initial capabilities against the reference implementation identified in the conversation. Refer to it neutrally in these documents. Do not create a named comparator record inside the project. Future comparisons should run in a separate pinned validation environment; SPmodels and its downstream workflows must work without the comparator installed. Do not replace computation with a hidden call to a comparator.

## Research and teaching

Maintain a consistent glm-like interface inspired by the owner's SDALGCP2, ordinary R methods, explicit assumptions, small reproducible examples and transparent uncertainty. The package is intended to expand to additional models; candidate models in the roadmap are not implementation instructions. Every implemented statistical and numerical method must have a full independent derivation accessible to MSc students, including assumptions, notation, intermediate steps, approximations, algorithms and links to the implementation. This applies to future methods as well as the initial collection. Follow docs/10_DERIVATIONS_AND_TEACHING_DATA.md; do not replace derivations with citations or final formulas alone.

The owner has authorised use of suitable datasets from lsda. Select examples that match each model's response and data support, record dataset versions and provenance, and keep synthetic fixtures for analytic and known-truth checks. Do not create a circular runtime dependency between SPmodels and lsda. This authorisation does not start model fitting, data copying or chapter writing in the current planning-only task.

Plan adoption in MBGapp, lsda and the course, including independent installation and student workflows. Plan for CRAN submission and easy installation, with optional R-universe distribution. No repository creation or publication is authorised by the current planning-only task.

## Evidence and scope

Read README.md and the numbered documents before implementation. Preserve the full agreed functional scope; prioritisation does not remove lower-priority capabilities. Reconcile the inventory and consumer needs before making completeness claims.

Do not claim numerical agreement without recorded results. Distinguish deterministic tolerances from stochastic equivalence, and never relax criteria merely to pass a failed comparison. Record discrepancies and independent mathematical checks. Additional models need their own validation where no equivalent reference model exists.

Use the Research/SPmodels directory as the working master. The Codex outputs copy is a delivery snapshot. Follow the owner's later decisions and applicable parent instructions.
