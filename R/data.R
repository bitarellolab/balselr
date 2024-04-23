#' #' Input data to run ncd2 example

#' An R object in data.table format
#' @name vcf_obj
#' @docType data
#' @format A data.table with 838 rows and 6 columns
#' #' @details
#'      details for vcf_obj
#' @source created from a VCF (Variant Call Format) input file generated with the SLiM simulator by reading it in with read_vcf().
#' @examples
#' data("vcf_obj")
#'
#' #' Input data to run ncd1 example

#' An R object in data.table format
#' @name ncd1_input
#' @docType data
#' @format A data.table with 838 rows and 6 columns
#' #' @details
#'      details for ncd1_input
#' \describe{
#'        \item{CHR}{Chromsome, single value - 1}
#'        \item{POS}{Position along the chromosome - first is 92, last is 29992}
#'        \item{REF}{Reference allele - A, C, G, T}
#'        \item{ALT}{Alternate allele - A, C, G, T}
#'        \item{tx_1}{Number of alternate allele copies - min=0, max=216}
#'        \item{tn_1}{Total number of alleles, single value - 216 }
#' }

#' @source created from a VCF (Variant Call Format) input file generated with the SLiM simulator by reading it in with parse_vcf(type="ncd1")
#' @examples
#' data("ncd1_input")


#' Input data to run ncd2 example

#' An R object in data.table format
#' @name ncd2_input
#' @docType data
#' @format A data.table with 838 rows and 8 columns
#' @details
#'      details for ncd2_input
#' \describe{
#'        \item{CHR}{Chromsome, single value - 1}
#'        \item{POS}{Position along the chromosome - first is 92, last is 29992}
#'        \item{REF}{Reference allele - A, C, G, T}
#'        \item{ALT}{Alternate allele - A, C, G, T}
#'        \item{tx_1}{Number of alternate allele copies in the ingroup - min=0, max=216}
#'        \item{tn_1}{Total number of alleles in the ingroup, single value - 216 }
#'        \item{tx_2}{Number of alternate allele copies in the outgroup - min=0, max=4}
#'        \item{tn_2}{Total number of alleles in the outgroup, single value - 216 }
#' }

#' @source created from a VCF (Variant Call Format) input file generated with the SLiM simulator by reading it in with parse_vcf(type="ncd2")
#' @examples
#' data("ncd2_input")
#'

