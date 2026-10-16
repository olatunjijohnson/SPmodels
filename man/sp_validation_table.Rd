\name{sp_validation_table}
\alias{sp_validation_table}
\title{Structured validation table}
\usage{
sp_validation_table(object, detail = c("summary", "fold"))
}
\arguments{
\item{object}{An SPmodels cross-validation, simulation-assessment,
simulation-recovery, or multi-scenario recovery object.}
\item{detail}{Return overall summaries or fold summaries for cross-validation.}
}
\value{A long-form data frame with assessment, scope, identifier, metric and
value columns.}
\description{
Provides a common machine-readable metric schema for SPmodels validation and
simulation-assessment results.
}
