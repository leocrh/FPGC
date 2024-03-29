% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NASAPOWER.R
\name{periodicextraction}
\alias{periodicextraction}
\title{Summarize climdata outputs}
\usage{
periodicextraction(
  df = NULL,
  climate.data = NULL,
  total.days = 120,
  by.periods = 10,
  start.date = "sowing.d",
  loc.id = "Loc_no",
  Date = "Date",
  mv.names = 9:18,
  pp.name = "pp"
)
}
\arguments{
\item{df}{A dataframe of one row (one location) that contains the date from which meteorological variables (MVs) will be extracted.}

\item{climate.data}{A list of of locations. Each element of the list must contain the MVs.}

\item{total.days}{The number of days from which the MVs will be extracted. Is the number of days after the initial date that will be consider to extract the data.}

\item{by.periods}{The number of days which the MVs will be averaged over. i.e, it must be multiple and <= than total.days.}

\item{start.date}{The column name in df that contains the date from which the extraction will start}

\item{loc.id}{The column name in df that contains the location ids}

\item{Date}{The column name in climate.data that contains the date of the recorded MVs}

\item{mv.names}{The index of the columns that contain the MVs in climate.data}

\item{pp.name}{The name of the column that contains the rainfall. This is used since percipitation is summed over the by.periods, instead of averaged as the other MVs}
}
\description{
Function to summarize the climate date by groups of days for a single location. It requires the locs.merged.climate object which is the climate data.
Dates must be in the appropriate format otherwise errors are shown.
}
