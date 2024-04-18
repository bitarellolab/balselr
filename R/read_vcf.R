#' Read vcf
#'
#' @param x The path and name of a vcf file
#' @return Returns a data.table object containing only SNPs
#' @export
#'
#' @examples read_vcf(x=system.file(package="balselr", "example.vcf"))
#' @import data.table
#' @importFrom data.table ":="
#'
read_vcf <- function(x = "inst/example.vcf") {
        inp <- data.table::fread(x, skip = "##", header = T)
        data.table::setnames(inp, "#CHROM", "CHR")
        #inp <-
        inp[REF %in% c("A", "C", "T", "G")][ALT %in% c("A", "C", "T", "G")]
}
