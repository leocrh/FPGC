% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SNPformat.R
\name{rmmissingsnp}
\alias{rmmissingsnp}
\title{Remove missing from X}
\usage{
rmmissingsnp(X = NULL, threshold = 20, plot.missing = F)
}
\arguments{
\item{X}{Is a marker matrix with markers in columns and gids in rows}

\item{threshold}{Is a threshold for the \% of the missing values in X}

\item{plot.missing}{Is a Logical to indicate if the the histogram of missing values is plotted. Default is False}
}
\description{
This function removes the missing marker data at certain threshold and plots  a histogram of the missing data
}
