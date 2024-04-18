#' Parse vcf
#' @param infile The path and name of a vcf file - single chromosome. Can be NULL if vcf_data is provided.
#' @param vcf_data A data.table object containing VCF data. Can be NULL if infile is provided.
#' @param n0 Number of individuals in pop0. Must be >=2 for both ncd1 and ncd2.
#' @param n1 Number of individuals in pop1. Must be >=1 for ncd2. Must be NULL for ncd1. NULL by default.
#' @param type Which input format will be generated. Options: ncd1, ncd2.
#' @param verbose Logical. Printed updates. Default is T.
#' @return Returns a data table object with relevant columns based on type
#' @export
#' @examples parse_vcf(vcf_data=your_preloaded_vcf, n0=108, type="ncd1")
#' @examples parse_vcf(infile=system.file(package="balselr", "example.vcf"), n0=108, n1=1, type="ncd2")
parse_vcf <- function(infile = NULL,
                      vcf_data = NULL,
                      n0 = NULL,
                      n1 = NULL,
                      type = NULL,
                      verbose=F) {
        if (!is.null(vcf_data)) {
                assertthat::assert_that(is(vcf_data, "data.table"), msg = "vcf_data must be a data.table object.\n")
                inp <- vcf_data
        } else {
                assertthat::assert_that(file.exists(infile), msg = glue::glue("VCF file {infile} does not exist.\n"))
                inp <- read_vcf(x = infile)
        }

        nind <- c(n0, n1)
        type <- tolower(type)
        assertthat::assert_that(type %in% c("ncd1", "ncd2"), msg = "Type must be either 'ncd1' or 'ncd2'.\n")
        if (type == "ncd2") {
                assertthat::assert_that(!is.null(n1), msg = "n1 cannot be 'NULL'. 'ncd2' requires an outgroup with at least one individual.")
        }

        index.col <- which(colnames(inp) == "FORMAT") + 1
        if (type == "ncd2") {
                res <-
                        balselr:::.vcf_ncd2(
                                x = inp,
                                nind = nind,
                                index.col = index.col,
                                verbose = verbose
                        )
        } else if (type == "ncd1") {
                res <-
                        balselr:::.vcf_ncd1(
                                x = inp,
                                nind = nind,
                                index.col = index.col,
                                verbose = verbose
                        )
        }
        return(res)
}
