% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modified_read_vcf.R
\name{modmod_read_vcf}
\alias{modmod_read_vcf}
\title{Read vcf}
\usage{
modmod_read_vcf(x = "inst/example.vcf", pos.range = NULL, id.range = NULL)
}
\arguments{
\item{x}{The path and name of a vcf file}
}
\value{
Returns a data.table object containing only SNPs
}
\description{
Read vcf
}
\examples{
read_vcf(x=system.file(package="balselr", "example.vcf"))
}
