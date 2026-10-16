\name{sp_intervention}
\alias{sp_intervention}
\alias{print.spmodels_intervention}
\title{Specify a Decay-Adjusted Intervention History}
\usage{
sp_intervention(
  history,
  unit,
  time,
  alpha,
  decay,
  time_units,
  estimate = FALSE,
  alpha_penalty = NULL,
  treated = NULL
)

\method{print}{spmodels_intervention}(x, ...)
}
\arguments{
\item{history}{Data frame with one row per unit/intervention time.}
\item{unit}{Name of the unit-identifying column in history and survey data.}
\item{time}{Name of the intervention-time column in history and survey-time column in survey data.}
\item{alpha}{Immediate proportional reduction between zero and one.}
\item{decay}{Positive decay time in the same units as \code{time}.}
\item{time_units}{Character label for the units shared by \code{time} and \code{decay}.}
\item{estimate}{Estimate \code{alpha} and \code{decay} during Laplace fitting?}
\item{alpha_penalty}{Optional two positive beta-shape parameters for a normalized alpha penalty.}
\item{treated}{Optional name of a zero-one history column.}
\item{x}{An intervention specification.}
\item{...}{Unused.}
}
\value{An object of class \code{spmodels_intervention}.}
\description{Defines a fixed intervention history. Every intervention strictly
before a survey multiplies baseline prevalence by
\code{1 - alpha * exp(-(survey_time - intervention_time) / decay)}.
Same-time interventions are excluded, corresponding to a survey immediately
before rollout.}
