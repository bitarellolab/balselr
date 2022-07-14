#' Read vcf
#'
#' @param x The path and name of a vcf file
#' @param only.bi Logical. If TRUE, only bi-allelic SNP positions from VCF file
#' are kept.
#' @param inds A vector specifying which individuals (samples) from vcf file
#' should be kept. If "all", all are kept.
#'
#' @return Returns a data.table object
#' @export
#'
#' @examples read_vcf(x="inst/example.vcf", only.bi=T, inds="all")
#' @import data.table
#' @importFrom data.table ":="
read_vcf <- function(x = "inst/example.vcf",
                     only.bi = T,
                     inds = "all") {
        REF <- ALT <- x2 <- NULL
        inp <- data.table::fread(x, skip = "##", header = T)
        data.table::setnames(inp, "#CHROM", "CHR")
        if (only.bi == T) {
                inp <-
                        inp[REF %in% c("A", "C", "T", "G")][ALT %in% c("A", "C", "T", "G")]
        } else {
                return(inp)
        }

        formatcol <- which(colnames(inp) == "FORMAT")

        if (is.numeric(inds)) {
                cbind(inp[, 1:formatcol], x2[, ...inds])
        } else if (inds == 'all') {
                return(inp)
        }
}
