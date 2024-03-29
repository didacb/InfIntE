% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_absolute_observations.R
\name{get_absolute_observations}
\alias{get_absolute_observations}
\title{Get absolute observations
Generate logical clauses depicting the absolute abundance change between two samples for each OTU.}
\usage{
get_absolute_observations(otu_abundance, comparisons, exclusion = FALSE)
}
\arguments{
\item{otu_abundance}{list joining all the OTU information. Obtained using the join_abundances function.}

\item{comparisons}{list with the position of the pairs of samples to compare.}

\item{exclusion}{logical indicating if exclusion cases have to be considered for generating the logic clauses. Default is FALSE.}
}
\value{
character vector with the abundance clauses describing the absolute abundance change between samples.
}
\description{
Get absolute observations
Generate logical clauses depicting the absolute abundance change between two samples for each OTU.
}
\examples{
get_absolute_observations(otu_abundance, comparisons)
}
