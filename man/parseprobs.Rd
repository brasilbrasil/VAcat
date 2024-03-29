% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/catnetparsing.R
\name{parseprobs}
\alias{parseprobs}
\title{Parse a vector of nodes and probabilities into a structured list}
\usage{
parseprobs(levels, p)
}
\arguments{
\item{levels}{vector of levels in parent nodes}

\item{p}{vector of probabilities}
}
\description{
Take a vector of levels of parent nodes and recursively
parse a vector of probabilities into a structured list suitable
for \code{\link{cnNew}}.
}

