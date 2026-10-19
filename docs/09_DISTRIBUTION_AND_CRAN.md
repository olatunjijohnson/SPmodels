# SPmodels distribution and CRAN plan

## Release objective

The owner intends SPmodels for public use through CRAN and/or another route that makes installation straightforward. Plan for CRAN as the primary stable distribution route, with R-universe as a possible complementary or alternative repository and a versioned source repository for development. A versioned GitHub source repository was created on 16 October 2026. The owner selected the MIT License on 19 October 2026; public repository visibility and the first WebR Pages deployment are being completed. External submissions remain future work.

The release identity should communicate an independently implemented spatial modelling collection for research and teaching, with a consistent API and its own model documentation. Numerical reference comparisons provide evidence for particular statistical capabilities. The wider package must have a clearly documented, useful contribution.

## Release design constraints

The current [CRAN repository policy](https://cran.r-project.org/web/packages/policies.html) requires a non-trivial contribution and says: “Introduction of packages providing back-compatibility versions of already available packages is not allowed.” Plan a distinct modelling API and reusable teaching examples. Prefer migrating MBGapp and lsda to that API instead of publishing a renamed compatibility layer. Independent implementation and numerical agreement alone do not guarantee acceptance; review the actual release against the policy before submitting.

The same policy covers names, dependency availability, portability, documentation size and check runtime. It also says CRAN is not an appropriate venue for distributing course notes. Keep concise function help and model vignettes in SPmodels; host full lectures, exercises and large assets separately. Provide the complete MSc-level derivations as a versioned companion guide, linked from package help and concise vignettes, with an accessible web edition and a downloadable edition proposed. Check the policy again at submission time.

SPmodels remains the working project name. Check current and archived CRAN and
Bioconductor names without relying on capitalisation to distinguish them, and
assess clarity for students before settling the submitted identity. Name
similarity is a release consideration, not proof that a name is unavailable.

## Local CRAN-style audit on 15 October 2026

`R CMD check --as-cran --no-manual --no-build-vignettes` completed every code,
help, example and test stage without an error or warning under R 4.6.1 on
Ubuntu 24.04.4. It reported two expected release notes:

1. incoming feasibility identified a new submission, the development version
   `0.0.0.9000`, and the current non-FOSS placeholder `LICENSE`;
2. `AGENTS.md` and the project `docs` directory are non-standard top-level
   package content.

These are not hidden. The first requires the owner's final licence and release
version decisions. The second requires a release-layout decision: exclude
project-management records from the package tarball and move concise user
guides/derivations to standard package locations or a versioned companion.
The ordinary development check remains independently preserved with
`Status: OK`; a CRAN-ready claim is not made.

## Distribution routes

| Route | Proposed role | User experience to verify |
| --- | --- | --- |
| CRAN | Stable public releases, subject to acceptance | A standard install.packages call using the final accepted package name |
| R-universe | Optional pre-release, complementary or alternative distribution | A documented repository-based install with platform builds checked before recommending it to students |
| Versioned source repository | Development history, issues, release tags and reproducibility | Maintainer and contributor workflow; source installation when appropriate |
| WebR binary repository | Browser-based course use through WebAssembly package binaries | Build SPmodels and required dependencies, host the binary repository with browser-accessible CORS, then install from WebR; see [WebR setup](WEBR.md) |
| Documentation and course sites | Model tutorials, full course materials and support | One clear getting-started page linking to the appropriate installation route |

The WebR build workflow builds and smoke-tests the WebAssembly packages, then
publishes the resulting static repository through GitHub Pages. The owner has
authorized making the GitHub source repository public for this distribution.
Confirm the Pages URL and browser CORS access after the first deployment.

R-universe provides infrastructure for a personal package repository. Account ownership and repository configuration should be verified when publication begins. [R-universe setup documentation](https://docs.r-universe.dev/publish/set-up.html)

Do not show an installation command as currently working until the package is actually available through that route. Use one recommended stable installation path in student instructions and label development builds clearly.

## Proposed release checklist

1. Define the first release's implemented model set and contribution. State accurately whether the full initial validation programme is complete; an earlier partial release cannot claim complete coverage.
2. Finalise the package name, maintainer, authorship, licence, method citations, data provenance and supported R versions.
3. Choose required dependencies and keep optional teaching or plotting extras separate where practical. Test clean installation without the comparator installed.
4. Run R package checks and meaningful numerical tests on the supported platforms. Resolve errors and warnings and explain remaining notes. Keep expensive reference comparisons in a separate reproducible validation workflow.
5. Check small examples and vignettes, links, help, startup behaviour, errors and diagnostics. Ensure routine examples do not require credentials or an active internet connection.
6. Test installation and the first tutorial on fresh Windows, macOS and Linux environments. Record binary availability, system requirements and any source-build instructions; do not assume students have compilers.
7. Complete MBGapp, lsda and course regression checks for the workflows claimed to have migrated. Preserve a versioned teaching environment for each course offering.
8. Verify derivation coverage for every implemented method and publish the matching version of the companion guide. Prepare release notes, submission notes and a model-specific validation report. Review the current repository policies before the owner submits the release.
9. Maintain an issue tracker, a process for responding to check failures and a versioning policy for model/API changes. Plan updates around course stability where practical.

## Boundary between package and course

SPmodels should contain reusable statistical functionality, concise examples, documented small datasets and focused model vignettes. lsda and the course repository/site can retain lesson sequencing, extended lecture notes, exercises and teaching applications. This separation supports an expanding model collection while keeping student installation manageable.
