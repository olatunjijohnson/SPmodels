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
`https://olatunjijohnson.github.io/SPmodels/`, then installs the package from
that public URL in WebR and reruns the core smoke workflow. This is not evidence
that every optional workflow works in a browser.

## Make it available to a course

The built repository and a small explanatory landing page are served from
GitHub Pages at `https://olatunjijohnson.github.io/SPmodels/`. The package
index is anonymously accessible and the static host returns
`Access-Control-Allow-Origin: *` for the repository files. Install SPmodels
from a WebR session:

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

WebR package installation downloads prebuilt WebAssembly package binaries from
a repository; it does not compile SPmodels from its GitHub source checkout in
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
WebAssembly image digest recorded in the workflow. The current workflow also
checks installation from the deployed Pages URL after publishing. See the
[recorded successful build run](https://github.com/olatunjijohnson/SPmodels/actions/runs/39239845877).

An earlier attempt with the `v0.6.0` builder image failed during dependency
metadata resolution before compilation. The passing workflow pins the exact
builder image digest that succeeded; retain this evidence when updating the
build toolchain.

This is evidence for package installation and those core workflows in the
WebR runtime, not a test of every SPmodels feature or the complete browser
course page. The smoke test now includes the ordinary-model pre-fit variogram,
the conditional post-fit variogram, chunked point predictions, joint draws on
a deliberately tiny example, and target-to-data-frame conversion. The optional
`sf` dependency is not part of this build. The live
hosted install and core smoke check passed on 19 October 2026 at commit
`4b6422a`; the [successful GitHub Actions run](https://github.com/olatunjijohnson/SPmodels/actions/runs/39240436594)
repeats this check after each deployment. Before course use, still exercise the
exact tutorial in its target browser. Keep the WebR version/build image and
SPmodels commit with each course offering so results can be reproduced.
