# Reference-validation records

This directory predeclares validation cases and acceptance rules without
identifying or embedding an external implementation. The executable comparison
harness and its pinned external environment must live outside SPmodels. The
installed package, its examples and its tests must not require that environment.

For every run, the separate record must contain:

- case ID and SPmodels source checksum;
- R, platform and dependency versions;
- complete inputs, seed and configuration;
- both raw outputs before rounding;
- the predeclared criterion and measured discrepancy;
- pass, fail or inconclusive status; and
- an independent mathematical check when results disagree.

An inconclusive stochastic comparison is not a pass. Criteria must not be
relaxed after results are seen merely to make a failing case pass.

The local recovery and approximation studies for newly implemented extensions
are predeclared in
[`02_RECOVERY_AND_APPROXIMATION_DESIGN.md`](02_RECOVERY_AND_APPROXIMATION_DESIGN.md).
They do not replace the external reference cases in the initial matrix.
