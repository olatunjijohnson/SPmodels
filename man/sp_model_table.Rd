\name{sp_model_table}
\alias{sp_model_table}
\title{Structured model-estimate table}
\usage{
sp_model_table(object, include_smoothness = TRUE)
}
\arguments{
\item{object}{A fitted \code{spmodels_fit} object.}
\item{include_smoothness}{Include the fixed Matérn smoothness in the table?}
}
\value{A stable data frame with one row per mean or covariance parameter.}
\description{
Converts an SPmodels fit into a machine-readable estimate table without
requiring downstream code to inspect internal fit fields.
}
