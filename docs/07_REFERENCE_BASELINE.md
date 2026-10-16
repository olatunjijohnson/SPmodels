# SPmodels reference validation preparation

## Purpose

Prepare reproducible numerical validation without using comparator branding in SPmodels or its planning documents. The reference implementation is the one identified by the owner in the conversation. No named comparator manifest, source archive or comparison outputs are included in this project.

**Status: planning only. No numerical comparisons have been run.**

## Reference setup for a future implementation session

Confirm the selected reference build and configuration when implementation begins. Use a separate validation environment and record the immutable build identity, artifact checksums, R and dependency versions, operating system, numerical and spatial libraries, threading, locale and RNG settings there. Keep comparator identity out of package and project materials, as requested by the owner.

Pin each comparison baseline. Do not assume that a live documentation site and a release with the same version label have identical capabilities or defaults. Record the exact arguments and input checksums for each case, and preserve earlier comparison results when evaluating a newer build.

## Matters to resolve before comparison

- Reconcile the complete initial feature scope with actual callable capabilities and supported configurations.
- Check callable defaults when help descriptions and usage signatures differ.
- Define spatial variance, nugget, measurement error and grouping variance separately.
- Distinguish survey time from the Gaussian-process temporal index in intervention models.
- Identify the estimator, likelihood constants, parameter transformations and uncertainty type for every model.
- Check consumer outputs against these definitions rather than relying on broad inference labels.

The detailed test layers, precision requirements and equivalence rules are in the validation protocol. These requirements remain mandatory even though comparator names are omitted from the project.

## Interface and distribution references

The installed SDALGCP2 0.1.1 help was used to confirm the owner's compact formula/data workflow. Its [maintainer README](https://github.com/olatunjijohnson/SDALGCP2) provides the public interface reference. Use [R's glm documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/glm.html) for ordinary formula, family, weight and offset semantics.

The [CRAN repository policy](https://cran.r-project.org/web/packages/policies.html) and [R-universe setup documentation](https://docs.r-universe.dev/publish/set-up.html) inform the future release plan. Recheck current requirements before submission. No release route is currently available for SPmodels.

## Naming

SPmodels is the owner's working project name. An existing spatial-model package
has a similar singular name; review the public identity and repository name
availability before submission. No name has been reserved.
