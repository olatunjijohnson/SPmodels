# SPmodels 0.1.1

- Allow the nugget variance to be fixed independently while other covariance
  parameters are estimated. `nugget_variance = 0` now fits the no-nugget
  boundary model, and `estimate_nugget = FALSE` fixes any supplied
  non-negative nugget.
- Propagate fixed-nugget parameter maps through dense Gaussian, Laplace,
  Monte Carlo maximum-likelihood, nearest-predecessor, prediction-uncertainty,
  simulation-recovery, likelihood and reporting workflows.
- Report fixed covariance parameters without invented standard errors or
  likelihood degrees of freedom.
- Correct the transformed-covariance ordering used by predictive parameter
  draws when grouped random slopes and temporal scale are both estimated.
- Add native and WebAssembly regression coverage for no-nugget models.

# SPmodels 0.1.0

- First tagged course release.
