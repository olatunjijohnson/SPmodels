# Using SPmodels with WebR

SPmodels can be built for WebR as a WebAssembly R package. WebR cannot install
an ordinary source checkout: package binaries must first be built for
WebAssembly. The workflow in
[`build-webr.yml`](../.github/workflows/build-webr.yml) builds SPmodels and its
required package dependencies into a CRAN-like WebR repository. It saves that
repository as a private GitHub Actions artifact and does not publish the
artifact or enable GitHub Pages.

## Build the package repository

The workflow runs after relevant changes are pushed to `main`, or can be
started manually from the repository's **Actions** tab using **Build SPmodels
for WebR**. A successful run provides an artifact named
`SPmodels-webr-<commit>`. Download and extract it to obtain the WebR package
repository. This is a build step, not evidence that the package has been
tested or that every optional workflow works in a browser.

## Make it available to a course

The extracted repository must be served from a static HTTPS host that the
students' browsers can access. The host must permit cross-origin requests
(CORS) from the tutorial page. Do not point WebR at the private source GitHub
repository or at a GitHub Actions artifact URL. Those are not a stable,
anonymous package repository. Keep the source repository private; publish only
the built repository to a host whose access and exposure have been approved.

Once the package repository is hosted, install SPmodels from a WebR session:

```r
webr::install("SPmodels", repos = "https://YOUR-HOST/wasm-repo/")
library(SPmodels)
```

Or from JavaScript using the WebR API:

```js
import { WebR } from "webr";

const webR = new WebR();
await webR.init();
await webR.installPackages("SPmodels", {
  repos: "https://YOUR-HOST/wasm-repo/"
});
await webR.evalR("library(SPmodels)");
```

The URL is a placeholder until a repository has actually been hosted. WebR
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

The first CI build and WebR runtime smoke check passed on 18 October 2026 at
SPmodels commit `939f03a`. The build compiled SPmodels and its required
dependencies; the WebR session installed the binary package, loaded it,
constructed a non-stationary Matérn correlation matrix, fit a small fixed-
covariance Gaussian model, predicted at a new point, and confirmed the
WebAssembly sequential-only parallel guard. See the [recorded GitHub Actions
run](https://github.com/olatunjijohnson/SPmodels/actions/runs/39161391488).

This is evidence for package installation and those core workflows in the
WebR runtime, not a test of every SPmodels feature or the complete browser
course page. The optional `sf` dependency is not part of this build. No WebR
binary repository is currently hosted for students. Before course use, host
the built repository with suitable CORS access and exercise the exact tutorial
in the target browser. Keep the WebR version/build image and SPmodels commit
with each course offering so results can be reproduced.
