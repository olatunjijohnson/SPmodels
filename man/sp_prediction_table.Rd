\name{sp_prediction_table}
\alias{sp_prediction_table}
\title{Structured prediction-summary table}
\usage{
sp_prediction_table(object, threshold = NULL, probs = c(0.025, 0.5, 0.975))
}
\arguments{
\item{object}{An object returned by an exact or nearest-predecessor
prediction-draw or prediction-target function.}
\item{threshold}{Optional finite exceedance threshold for prediction draws.}
\item{probs}{Numeric quantile probabilities between zero and one.}
}
\value{A data frame of draw summaries and prediction metadata.}
\description{
Summarises joint prediction draws at each location or returns an existing
grouped-target summary using a stable downstream schema.
}
