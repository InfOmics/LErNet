% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LErNet.R
\name{enrich}
\alias{enrich}
\title{Functional enrichment}
\usage{
enrich(entrez_genes, organism)
}
\arguments{
\item{entrez_genes}{list of proteins, in Entrez format, for which to compute the enrichment}

\item{oganism}{oganism name (see \code{\link[ReactomePA]{enrichPathway}})}
}
\value{
a ReactomePA result object.
}
\description{
Computes functional enrichment of a given set of proteins via the ReactomePA package.
}
