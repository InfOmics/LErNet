% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/support.R
\name{load_gtf}
\alias{load_gtf}
\title{Creation of a coordinates data.frame from a GTF file}
\usage{
load_gtf(gtf_file)
}
\arguments{
\item{gtf_file}{the path tot he GTF file}
}
\value{
a data.frame with columns: \code{id type seqname start end}
}
\description{
Creates the coordinates data.frame by reading the data from a GTF file
having 9 columns (whihc is the typical format of GTF files from GENCODE).
}
\examples{
gtf_file <- system.file("extdata", "gencode.v34.annotation.gtf.gz", package = "LErNet")
complete_positions <- LErNet::load_gtf(gtf_file)

}
