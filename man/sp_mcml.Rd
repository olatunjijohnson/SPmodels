\name{sp_mcml}
\alias{sp_mcml}
\title{Refine a Non-Gaussian Fit by Monte Carlo Maximum Likelihood}
\usage{
sp_mcml(object, control = sp_mcml_control(), seed = NULL, verbose = FALSE)
}
\arguments{
  \item{object}{A fitted binomial-logit or Poisson-log
  \code{spmodels_fit}.}
  \item{control}{Settings returned by \code{sp_mcml_control()}.}
  \item{seed}{Optional integer seed local to this call.}
  \item{verbose}{Whether to trace optimization.}
}
\value{An updated \code{spmodels_fit}. The \code{mcml} component contains the
proposal, retained importance draws and weights, effective sample size,
maximum normalized weight and an approximate Monte Carlo standard error. When
adaptation is requested it also retains the original moments and complete
pilot record.}
\description{Replaces the Laplace fitting objective with an
importance-sampling estimate of the binomial or Poisson marginal likelihood.}
\details{The importance distribution is either a scaled Gaussian
approximation to the latent posterior at the supplied fit or a normalized
same-centre mixture of declared scales. Mixture draws are balanced across
components and divided by the complete mixture density. The proposal is normalized,
so the returned Monte Carlo log likelihood is absolute rather than known only
up to an omitted constant. Common random numbers make optimization deterministic
within a call. A small effective sample size or dominant normalized weight
indicates that the importance approximation is unreliable. Mixture likelihood
MCSE uses its component-stratified variance. With \code{adapt = TRUE}, an
independent pilot supplies importance-weighted moments that are shrunk toward
the Laplace moments. The pilot is diagnosed, retained for audit and discarded
before a new main sample is drawn. Conditional on the pilot, weighting by the
complete normalized adapted density preserves the same marginal likelihood.
If pilot ESS is below \code{ess_warning}, the original moments are retained and
the recorded adaptation status explains the fallback. Adaptation is not
guaranteed to improve overlap. Downstream latent prediction
continues to use local Gaussian conditional calculations at the MCML parameter
estimates. If intervention impact and decay were declared
estimated, they are re-estimated on logit and log scales and their optional
normalized beta penalty is included in the MCML criterion.}
\examples{
d <- data.frame(y = c(0, 1, 0, 1), x = c(0, 1, 0, 1),
                east = c(0, 1, 0, 1), north = c(0, 0, 1, 1))
fit <- sp_fit(y ~ x, d, family = binomial(), coords = ~ east + north,
  control = sp_control(spatial_variance = 1, scale = 1,
    nugget_variance = 0.2, estimate = FALSE, distance_units = "km"))
corrected <- sp_mcml(fit,
  control = sp_mcml_control(nsim = 200, maxit = 20), seed = 1)
corrected$mcml$effective_sample_size
}
