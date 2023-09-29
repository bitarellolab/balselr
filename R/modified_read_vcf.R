
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
        # Read VCF file
        inp <- data.table::fread(x, skip = "##", header = TRUE)
        data.table::setnames(inp, "#CHROM", "CHR")

        inp <- inp[REF %in% c("A", "C", "T", "G") & ALT %in% c("A", "C", "T", "G")]
        if (!is.null(id.range)) {
                if (length(id.range) != 1 && length(id.range) != 2) {
                        stop("id.range must have either 1 or 2 values.")
                }

                if (length(id.range) == 2 && id.range[1] > id.range[2]) {
                        stop("The first value of id.range must be less than or equal to the second value.")
                }

                basic_cols <- c("CHR", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")

                if (length(id.range) == 2) {
                        sample_cols <- paste0("i", (id.range[1]-1):(id.range[2]-1))
                } else {
                        sample_cols <- paste0("i", id.range-1)
                }

                if (!all(sample_cols %in% names(inp))) {
                        stop("Some specified sample columns do not exist in the VCF file.")
                }

                inp <- inp[, c(basic_cols, sample_cols), with = FALSE]
        }

        return(inp)
}
