\name{sp_sampler_diagnostics}
\alias{sp_sampler_diagnostics}
\alias{plot.spmodels_sampler_diagnostics}
\alias{print.spmodels_sampler_diagnostics}
\title{Diagnose Retained Latent Samples}
\usage{
sp_sampler_diagnostics(object, component = 1L, max_lag = NULL, chains = NULL)

\method{plot}{spmodels_sampler_diagnostics}(x, ...)

\method{print}{spmodels_sampler_diagnostics}(x, ...)
}
\arguments{
  \item{object}{An object returned by \code{sp_latent_sample()}.}
  \item{component}{Component index or name selected for trace and ACF output.}
  \item{max_lag}{Positive maximum lag smaller than retained sample size.}
  \item{chains}{Optional distinct chain indices or names. The default uses all
  available chains.}
  \item{x}{A sampler-diagnostics object.}
  \item{...}{Additional graphical arguments for the trace panel.}
}
\value{A \code{spmodels_sampler_diagnostics} object containing per-component
means, standard deviations, multi-chain effective sample sizes, Monte Carlo
standard errors and rank-normalized split R-hat values, plus chain acceptance
rates and selected trace/autocorrelation values. Methods invisibly return the
object.}
\description{Calculates autocorrelation-aware effective sample sizes and Monte
Carlo standard errors. With multiple chains it also calculates the maximum of
rank-normalized split and folded rank-normalized split R-hat.}
