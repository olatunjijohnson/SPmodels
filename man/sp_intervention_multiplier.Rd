\name{sp_intervention_multiplier}
\alias{sp_intervention_multiplier}
\title{Calculate Decay-Adjusted Intervention Multipliers}
\usage{sp_intervention_multiplier(intervention, data)}
\arguments{
\item{intervention}{An object returned by \code{sp_intervention()}.}
\item{data}{A data frame containing the declared unit and survey-time columns.}
}
\value{One multiplier in \code{[0, 1]} per row of \code{data}.}
\description{Evaluates the intervention-history multiplier for survey or prediction rows.}
