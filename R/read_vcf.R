#' Read vcf
#'
#' @param x The path and name of a vcf file
#' @param only.bi Logical. If TRUE, only bi-allelic SNP positions from VCF file
#' are kept.
#' @return Returns a data.table object
#' @export
#'
#' @examples read_vcf(x="inst/example.vcf", only.bi=T)
#' @import data.table
#' @importFrom data.table ":="
read_vcf <- function(x = "inst/example.vcf",
                     only.bi = T) {
        REF <- ALT <- x2 <- NULL
        inp <- data.table::fread(x, skip = "##", header = T)
        data.table::setnames(inp, "#CHROM", "CHR")
        if (only.bi == T) {
                inp <-
                        inp[REF %in% c("A", "C", "T", "G")][ALT %in% c("A", "C", "T", "G")]
                inp
        } else {
                inp
        }
}

