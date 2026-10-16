\name{sp_format_table}
\alias{sp_format_table}
\alias{print.spmodels_formatted_table}
\title{Format a stable SPmodels table for presentation}
\usage{
sp_format_table(
  object,
  digits = 3L,
  p_digits = 3L,
  missing = "—",
  detail = c("summary", "fold")
)

\method{print}{spmodels_formatted_table}(x, ..., row.names = FALSE)
}
\arguments{
\item{object}{A data frame or an SPmodels fit, prediction, validation or
simulation result supported by the stable table functions.}
\item{digits}{Decimal places for ordinary numeric columns.}
\item{p_digits}{Decimal places and lower display threshold for
\code{p_value} columns.}
\item{missing}{Character representation of missing values.}
\item{detail}{Overall or fold detail for cross-validation.}
\item{x}{A formatted SPmodels table.}
\item{...}{Further arguments passed to \code{print.data.frame()}.}
\item{row.names}{Display row names?}
}
\value{A character data frame of class \code{spmodels_formatted_table}. The
exact unrounded source remains available in the \code{raw_table} attribute.}
\description{Creates a presentation-only view over a stable SPmodels table. It
never changes the machine-readable values from which the view was produced.}
