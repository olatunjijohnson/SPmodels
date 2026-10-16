\name{sp_report}
\alias{sp_report}
\alias{print.spmodels_report}
\title{Assemble an SPmodels analysis report}
\usage{
sp_report(
  object,
  prediction = NULL,
  validation = NULL,
  title = "SPmodels analysis report",
  digits = 3L,
  p_digits = 3L,
  validation_detail = c("summary", "fold")
)

\method{print}{spmodels_report}(x, ...)
}
\arguments{
\item{object}{A fitted \code{spmodels_fit} object.}
\item{prediction}{Optional joint-prediction draws or prediction-target result.}
\item{validation}{Optional cross-validation, predictive-assessment or
simulation-recovery result.}
\item{title}{One non-empty report title.}
\item{digits, p_digits}{Presentation precision passed to
\code{sp_format_table()}.}
\item{validation_detail}{Overall or fold-level cross-validation output.}
\item{x}{An assembled \code{spmodels_report}.}
\item{...}{Reserved for print-method compatibility.}
}
\value{An \code{spmodels_report} containing analysis metadata, unrounded
stable tables and separate formatted views.}
\description{Assembles model, optional prediction and optional validation
sections while keeping presentation formatting separate from all unrounded
statistical output.}
