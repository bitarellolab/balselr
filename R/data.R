#' Input data to run ncd1 example

#' An R object in data.table format created from a VCF (Variant Call Format) input file generated with the SLiM simulator.

#' @format A data.table with 838 rows and 6 columns
#' \describe{
#'        \item{CHR}{Chromsome, single value - 1}
#'        \item{POS}{Position along the chromosome - first is 92, last is 29992}
#'        \item{REF}{Reference allele - A, C, G, T}
#'        \item{ALT}{Alternate allele - A, C, G, T}
#'        \item{tx_1}{Number of alternate allele copies - min=0, max=216}
#'        \item{tn_1}{Total number of alleles, single value - 216 }
#' }

#' @examples
#' data(ncd1_input)
