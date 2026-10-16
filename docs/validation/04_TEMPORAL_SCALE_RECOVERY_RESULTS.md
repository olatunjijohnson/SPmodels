# Non-Gaussian temporal-scale recovery results

## Scope and result

This report records the TSC cases predeclared in
`02_RECOVERY_AND_APPROXIMATION_DESIGN.md`. The study estimates a temporal scale
jointly with spatial covariance and regression parameters in binomial-logit
and Poisson-log latent Gaussian models. It is an internal synthetic recovery
study, not an external reference comparison. Raw fits, MCML repeats, summaries
and the integrity-check script are retained under `data/` and `tools/`.

All four scenarios reached the cap of 400 response attempts. Usability ranged
from 76.8% to 83.0%, below the predeclared 90% requirement. The estimated
natural temporal scales were extremely dispersed and right-skewed; every
scenario missed its temporal-scale bias-MCSE threshold by several orders of
magnitude. Thus all four cases are inconclusive, not recovery successes. No
threshold, start, warning rule or optimizer setting was changed after seeing
results.

| Case | Family | Truth | Attempts / usable | Mean / median estimate | Bias-MCSE / target | At upper bound | Local Wald coverage* | Result |
|---|---|---:|---:|---:|---:|---:|---:|---|
| TSC-01 | Binomial, 10 trials | 0.75 | 400 / 332 (83.0%) | 615.807 / 0.116 | 267.982 / 0.050 | 5 / 332 | 100.0% (178) | Inconclusive |
| TSC-02 | Binomial, 10 trials | 2.00 | 400 / 313 (78.3%) | 1,396.709 / 0.984 | 389.838 / 0.100 | 9 / 313 | 88.2% (211) | Inconclusive |
| TSC-03 | Poisson | 0.75 | 400 / 307 (76.8%) | 915.054 / 0.132 | 341.304 / 0.050 | 7 / 307 | 100.0% (143) | Inconclusive |
| TSC-04 | Poisson | 2.00 | 400 / 307 (76.8%) | 1,830.440 / 0.777 | 452.270 / 0.100 | 11 / 307 | 94.1% (136) | Inconclusive |

*Local Wald coverage is calculated only where finite Hessian-based scale
uncertainty was available. Its unusually broad intervals do not establish
useful or calibrated uncertainty; these values are descriptive.

## Diagnostics

Every unusable fit (341 of 1,600) had optimizer convergence code 52 and an
associated warning; there were no fit errors. The remaining 1,259 fits had
finite target estimates, a positive-definite fitted covariance, converged
latent mode, zero optimizer convergence code and no warning. Hessian
covariance uncertainty was finite in 668 fits. These unsuccessful optimizer
records remain in the attempt denominator and raw file.

Natural-scale temporal estimates hit both optimizer bounds. The registered
bounds for these time points are `0.0001` and `40,000`; the upper bound was hit
by 5, 9, 7 and 11 usable fits, respectively, while the lower bound was hit by
80, 35, 49 and 44. Among usable fits, the 90th percentile was 1.654, 9.225,
2.071 and 10.55, but the respective maxima were all 40,000. The enormous
means and bias-MCSEs therefore reflect a substantial tail and boundary
behavior, not a stable estimate near truth. Median fitted spatial scales were
1.307, 0.949, 1.122 and 0.783 against generating value 2; their 90th
percentiles were 28.3, 10.0, 33,715 and 101.2. This indicates broader
covariance-parameter instability in these cases.

MCML was run on the first 25 usable responses of each case, with three
independent 1,000-draw repeats per response (300 repeat rows total). A response
set passes only when all three repeats satisfy the predeclared ESS, relative
MCSE, optimizer and warning gates.

| Case | Response sets with all three repeats passing | Median ESS fraction | Median maximum relative MCSE | Repeat rows with warnings |
|---|---:|---:|---:|---:|
| TSC-01 | 4 / 25 (16%) | 0.285 | 0.050 | 60 / 75 |
| TSC-02 | 6 / 25 (24%) | 0.330 | 0.045 | 48 / 75 |
| TSC-03 | 1 / 25 (4%) | 0.071 | 0.115 | 57 / 75 |
| TSC-04 | 2 / 25 (8%) | 0.159 | 0.073 | 54 / 75 |

This does not support a claim that the tested MCML correction is operationally
stable for these temporal-scale cases. The MCML gate results do not rescue the
poor Laplace recovery, and neither route is considered validated here.

## Execution and reproducibility

The fits used archive SHA-256
`22d059ff921242770fc1110af93b7ce9670020a9ddee19b97193f0841431b374`, R 4.6.1
on x86_64 Ubuntu (`linux-gnu`) and BLAS
`/usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.0`. The first pass saved the
initial 100 usable refits per case and all 100 MCML response sets, but its
process exited with status 1 after writing the final rows: the study script
had been edited while R was still reading it. No fitted rows were lost or
rerun. A resume-capable runner was then used to continue the saved scenarios
through the registered 200/400 stages, preserving the original archive hash
and all response seeds. That continuation exited successfully. First-pass and
continuation wall times
were 53:32.10 and 48:37.33; combined elapsed time was 1:42:09.43. Their
externally observed peak RSS values were 197,832 and 154,980 KiB, respectively.
These are whole-process measurements on one machine, not per-case memory
claims.

Recreate the summaries and check raw-row integrity with:

```sh
Rscript tools/summarize_temporal_scale_recovery.R \
  docs/validation/data/temporal_scale_recovery.csv \
  docs/validation/data/temporal_scale_mcml_repeats.csv \
  docs/validation/data/temporal_scale_recovery_summary.csv \
  docs/validation/data/temporal_scale_mcml_summary.csv
Rscript tools/validate_temporal_scale_recovery.R
```

The raw fit, raw MCML, fit-summary and MCML-summary SHA-256 values are
`0fb2a12b65712a1121c826cdb106261866a761d54985b51794c4903446f404a3`,
`fd5f89b8afb45c5445f77e60f0bd6f419975a19541da25992979717d090970e3`,
`808d45e984ea6db929cf80016e37ec9d6183a696e79afed4e1d2f9bc6e4f207e` and
`c75db1d40e29439f964aa447c8c33f541830f13adbc702802b4fdc4dda936e66`.

The next scientific step is to diagnose temporal-scale identifiability and
optimizer behavior independently, without changing these acceptance rules.
The separate MCML proposal-policy and nearest-predecessor prediction studies
remain outstanding. No reference comparison was performed.
