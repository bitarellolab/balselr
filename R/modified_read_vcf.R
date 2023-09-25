
#' Read vcf
#'
#' @param x The path and name of a vcf file
#' @param id.range A numeric vector specifying the IDs or range of IDs to keep.
#' @return Returns a data.table object containing only SNPs
#' @export
#'
#' @examples read_vcf(x=system.file(package="balselr", "example.vcf"))
#' @import data.table
#' @importFrom data.table ":="
#'
mod_read_vcf <- function(x = "inst/example.vcf", id.range = NULL) {
    inp <- data.table::fread(x, skip = "##", header = T)
    data.table::setnames(inp, "#CHROM", "CHR")

    inp <- inp[REF %in% c("A", "C", "T", "G")][ALT %in% c("A", "C", "T", "G")]

    if (!is.null(id.range)) {
        if (length(id.range) == 2) {
            cols_to_keep <- c("CHR", "POS", "ID", "REF", "ALT", (6 + id.range[1]):(5 + id.range[2]))
        } else {
            cols_to_keep <- c("CHR", "POS", "ID", "REF", "ALT", 5 + id.range)
        }
        inp <- inp[, ..cols_to_keep]
    }

    return(inp)
}
