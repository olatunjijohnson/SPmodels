# MCML proposal-policy study results

## Scope and status

This report records the synthetic MCP policy matrix predeclared in
`02_RECOVERY_AND_APPROXIMATION_DESIGN.md`. Two fixed response fixtures were
used: a 24-site binomial case and a 48-site Poisson case. Both base Laplace fits
converged with no warnings. Each of four proposal policies was run five times
at 500, 1,000 and 2,000 main draws per case, for 120 independent MCML
attempts. The raw fixture rows and attempts are retained in `data/`.

No policy was operational for either synthetic case under the predeclared
four-of-five ESS-fraction / relative-MCSE gate at 2,000 draws. This is a failed
policy study, not a reason to select the best-looking seed or simply increase
draws. Every optimizer converged, but optimizer convergence did not imply
reliable importance weights.

| Case | Policy | Attempts passing at 2,000 | Median ESS fraction | Median relative MCSE | Optimizer convergence | Attempts with warnings |
|---|---|---:|---:|---:|---:|---:|
| MCP-01 binomial, n=24 | Single fixed | 2 / 5 | 0.0069 | 0.269 | 5 / 5 | 3 / 5 |
| MCP-01 binomial, n=24 | Two-scale fixed | 2 / 5 | 0.1510 | 0.052 | 5 / 5 | 0 / 5 |
| MCP-01 binomial, n=24 | Single adaptive | 1 / 5 | 0.0294 | 0.129 | 5 / 5 | 4 / 5 |
| MCP-01 binomial, n=24 | Two-scale adaptive | 1 / 5 | 0.0559 | 0.092 | 5 / 5 | 4 / 5 |
| MCP-02 Poisson, n=48 | Single fixed | 2 / 5 | 0.0690 | 0.082 | 5 / 5 | 3 / 5 |
| MCP-02 Poisson, n=48 | Two-scale fixed | 0 / 5 | 0.0009 | 0.745 | 5 / 5 | 5 / 5 |
| MCP-02 Poisson, n=48 | Single adaptive | 0 / 5 | 0.0227 | 0.147 | 5 / 5 | 5 / 5 |
| MCP-02 Poisson, n=48 | Two-scale adaptive | 0 / 5 | 0.0009 | 0.765 | 5 / 5 | 5 / 5 |

The fixed two-scale mixture gave the strongest median diagnostics for the
binomial fixture, but it still passed only two attempts. No tested policy was
adequate for Poisson. Pilot adaptation was used in five binomial attempts for
each adaptive policy; for Poisson it was used in five single-scale adaptive
attempts and zero two-scale adaptive attempts. These counts do not establish
that adaptation generally helps or hurts.

The noisy shared-parameter objectives are retained for every attempt. They
were recomputed at the common base-fit parameter vector using each attempt's
retained main sample and its complete normalized proposal density. These are
importance-sampling estimates with their own Monte Carlo variability, not
exact objective values; the base-fit log likelihoods were -40.388 for MCP-01
and -69.587 for MCP-02. The optimized likelihood estimates and parameter
estimates are also present in the raw attempt table.

## Loa loa teaching-data diagnostic

The separately declared diagnostic used the checked `lsda` 0.2.0 `loaloa`
dataset (197 sites) with the same explicit Lambert azimuthal equal-area
projection and binomial model call as the package's teaching example. It is
not a known-truth case and is not used to select a policy. The file remains in
the lsda course package; the study script reads it by an explicit input path,
does not copy it into SPmodels, and adds no runtime dependency.

The four policies used five independent 2,000-draw attempts each. All 20
optimizations converged, but every attempt emitted a low-ESS warning. Median
ESS fractions / maximum weights / relative MCSEs were:

| Policy | Median ESS fraction | Median maximum weight | Median relative MCSE | Adaptation used |
|---|---:|---:|---:|---:|
| Single fixed | 0.00144 | 0.585 | 0.589 | 0 / 5 |
| Two-scale fixed | 0.00061 | 0.904 | 0.905 | 0 / 5 |
| Single adaptive | 0.00133 | 0.602 | 0.613 | 4 / 5 |
| Two-scale adaptive | 0.00058 | 0.926 | 0.927 | 0 / 5 |

One adaptive single-scale attempt had ESS fraction 0.1115, but its other
diagnostics did not make this 197-site case a stable-likelihood result. The
Loa case therefore reinforces the synthetic warning: finite optimizer output
is not evidence of adequate importance-sampling overlap. This remains a
descriptive teaching-data diagnostic and cannot establish known-truth
recovery or select a policy.

## Reproducibility

The synthetic MCML attempts used source archive SHA-256
`dd0c7f5d1532f70c9eeadf53f62558945dfa77eec6b4674c746f10f26608fce4`, R 4.6.1
on x86_64 Ubuntu, and BLAS
`/usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.0`. The 120-attempt run took
4:50.59 wall-clock and peaked at 150,964 KiB RSS under GNU `time`.

Recreate the synthetic summaries and integrity checks with:

```sh
Rscript tools/summarize_mcml_proposal_study.R
Rscript tools/validate_mcml_proposal_study.R
```

The fit, summary, operational-gate, fixture, and base-fit CSV SHA-256 values
are `6ff6738054172c444cbbbb129468202081e3098fbefa440f8a447b260ea2499f`,
`0f1f49d82d3b724e750e3170ebd48a8f0c59c0dcb514d90b6065744e1ea7df68`,
`2c2642f55322a281f5c5f3a6db99dfbaac90d9ce9a6a6586164eb9f694d86b0f`,
`cc627247b65c3351f0a24c56e4d4fe2eac269e6757a6eda642264e3b44297d9f` and
`5c04b9422355023aba292c0b3415098165385892be6ba3085a36be684d34b575`.
The Loa diagnostic used the same SPmodels source archive and the documented
lsda input-file SHA-256 above; it took 27:34.37 and peaked at 249,808 KiB RSS.
Its raw output SHA-256 is
`5b7e89977bec82d5bee8ded0534c2dcd577d3aa8cbb52cd0c9f53f0fcd985727`.

The immediate technical priority is a proposal redesign informed by the raw
ESS and weight diagnostics, followed by a newly predeclared independent study.
No policy is recommended from these results.
