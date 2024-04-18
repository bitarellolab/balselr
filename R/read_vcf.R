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
read_vcf <- function(x = "inst/example.vcf", pos.range = NULL, id.range = NULL) {
        inp <- data.table::fread(x, skip = "##", header = TRUE)
        data.table::setnames(inp, "#CHROM", "CHR")

        inp <- inp[REF %in% c("A", "C", "T", "G") & ALT %in% c("A", "C", "T", "G")]

        # if (!is.null(pos.range)) {
        #         if (length(pos.range) == 2) {
        #                 inp <- inp[POS >= pos.range[1] & POS <= pos.range[2]]
        #         } else {
        #                 inp <- inp[POS %in% pos.range]
        #         }
        # }
        #
        # col.range <- 1:9
        #
        # if (!is.null(id.range)) {
        #         if (length(id.range) == 2) {
        #                 col.range <- c(col.range, (10 + id.range[1] - 1):(10 + id.range[2] - 1))
        #         } else {
        #                 col.range <- c(col.range, 9 + id.range)
        #         }
        #         col.range <- col.range[col.range <= ncol(inp)]
        # }
        #
        # inp <- inp[, .SD, .SDcols = col.range]

        return(inp)
}
