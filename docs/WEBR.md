# Using SPmodels with WebR

SPmodels can be built for WebR as a WebAssembly R package. WebR cannot install
an ordinary source checkout: package binaries must first be built for
WebAssembly. The workflow in
[`build-webr.yml`](../.github/workflows/build-webr.yml) builds SPmodels and its
required package dependencies into a CRAN-like WebR repository, runs a WebR
smoke test, then publishes the repository through GitHub Pages.

## Build the package repository

The workflow runs after relevant changes are pushed to `main`, or can be
started manually from the repository's **Actions** tab using **Build SPmodels
for WebR**. A successful run deploys the repository to
`https://olatunjijohnson.github.io/SPmodels/`. This is a build and core runtime
smoke test, not evidence that every optional workflow works in a browser.

## Make it available to a course

The built repository is served from GitHub Pages. Install SPmodels from a WebR session:

```r
webr::install("SPmodels", repos = "https://olatunjijohnson.github.io/SPmodels/")
library(SPmodels)
```

Or from JavaScript using the WebR API:

```js
import { WebR } from "webr";

const webR = new WebR();
await webR.init();
await webR.installPackages("SPmodels", {
  repos: "https://olatunjijohnson.github.io/SPmodels/"
});
await webR.evalR("library(SPmodels)");
```

The URL becomes active after the first successful Pages deployment. WebR
package installation downloads prebuilt WebAssembly package binaries from a
repository; it does not compile SPmodels from its private source checkout in
the student's browser. See the [WebR package installation
guide](https://docs.r-wasm.org/webr/latest/packages.html) and [package build
guide](https://docs.r-wasm.org/webr/latest/building.html).

## Current package boundaries

- The core package workflows use R code and numeric coordinate data frames;
  these are the initial target for browser teaching examples.
- `sf` is a suggested, optional dependency, not an installation requirement
  for core workflows. Geometry-based features require a compatible WebR build
  of `sf` and its system dependencies; this build workflow does not promise or
  bundle that optional stack.
- The simulation-recovery utilities can use fork-based parallel workers on
  supported desktop systems. WebAssembly/browser sessions must use
  `workers = 1` for those utilities.
- Dense spatial covariance and fitting use memory and CPU that grow quickly
  with the number of locations. Browser lessons should use small, deliberately
  chosen examples. A successful package build alone does not establish
  practical browser performance or scientific validation.

## Status and next release gates

The reproducible CI build and WebR runtime smoke check passed on 18 October
2026 at SPmodels commit `180a190`. The build compiled SPmodels and its required
dependencies; the WebR 0.6.0 session installed the binary package, loaded it,
constructed a non-stationary Matérn correlation matrix, fit a small fixed-
covariance Gaussian model, predicted at a new point, and confirmed the
WebAssembly sequential-only parallel guard. The builder uses the exact
WebAssembly image digest recorded in the workflow. See the [recorded GitHub
Actions run](https://github.com/olatunjijohnson/SPmodels/actions/runs/39161861749).

An earlier attempt with the `v0.6.0` builder image failed during dependency
metadata resolution before compilation. The passing workflow pins the exact
builder image digest that succeeded; retain this evidence when updating the
build toolchain.

This is evidence for package installation and those core workflows in the
WebR runtime, not a test of every SPmodels feature or the complete browser
course page. The optional `sf` dependency is not part of this build. After
deployment, verify anonymous access, cross-origin loading from the course site,
and the exact tutorial in its target browser. Keep the WebR version/build image
and SPmodels commit with each course offering so results can be reproduced.
