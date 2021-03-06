% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_to_srplv.R
\name{data_to_srplv}
\alias{data_to_srplv}
\title{Data format for srplv
\code{data_to_srplv} converts the dataset to the srplv layout.}
\usage{
data_to_srplv(data, sys.col, time.col, status.col, n.comp)
}
\arguments{
\item{data}{a data frame with three columns, recording the system ID, time and status.}

\item{sys.col}{the index containing the system ID.}

\item{time.col}{the index containing the event times and end-of-observation times.}

\item{status.col}{the index containing censoring or failure status.}

\item{n.comp}{the number of components.}
}
\value{
A list with the following components:
\item{n}{the number of systems in the fleet.}
\item{n.comp}{the number of components in the systems.}
\item{r}{a vector containing the number of failures for each system in the fleet.}
\item{time}{a list containing the failure time for each system.}
\item{cens}{a vector that records the end-of-observation time for each system.}
}
\description{
Data format for srplv
\code{data_to_srplv} converts the dataset to the srplv layout.
}
