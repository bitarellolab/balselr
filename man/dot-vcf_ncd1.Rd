% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcf_ncd1.R
\name{.vcf_ncd1}
\alias{.vcf_ncd1}
\title{Make input file for NCD1 from VCF}
\usage{
.vcf_ncd1(x, nind = nind, index.col = index.col, verbose = T)
}
\arguments{
\item{x}{An object of class data.table containing typical elements of a vcf file. Output: data.table with
columns: CHR, POS, REF, ALT, tx_1 (number of alternate allele copies), tn_1 (total number of alleles) and (if type is ncd2) tx_2 and tn_2}

\item{nind}{A vector containing the number of diploid individuals from each
population to be used. For ncd1, only one population is used.}

\item{index.col}{First genotype column in VCF file.}

\item{verbose}{Logical. Printed updates. Default is T.}
}
\description{
Make input file for NCD1 from VCF
}
