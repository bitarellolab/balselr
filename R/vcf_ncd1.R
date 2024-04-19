#' Make input file for NCD1 from VCF
#'
#' @param x An object of class data.table containing typical elements of a vcf file. Output: data.table with
#' columns: CHR, POS, REF, ALT, tx_1 (number of alternate allele copies), tn_1 (total number of alleles) and (if type is ncd2) tx_2 and tn_2
#' @param nind A vector containing the number of diploid individuals from each
#' population to be used. For ncd1, only one population is used.
#' @param index.col First genotype column in VCF file.
#' @param verbose Logical. Printed updates. Default is T.
#'
#to do: some of the tests suggested in parse_vcf might pertain here.
.vcf_ncd1 <- function(x,
                      nind = nind,
                      index.col = index.col,
                      verbose = T) {
        npop <- length(nind)
        assertthat::assert_that(npop == 1, msg = "NCD1 only uses one species. 'nind' should have length 1. \n")
        # to do: asssert that input file looks the way it should
        # to do: assert that index.col is an integer>1 and that the column it corresponds to is not named anything like CHR, POS, etc (i.e, the typical vcf columns)

        pop0_cols <-
                colnames(x)[index.col:(index.col + (nind[1] - 1))]

        if (verbose == T) {
                cat(glue::glue("Only {npop} population will be considered:"),
                    "\n")
        }
        if (verbose == T) {
                cat(
                        glue::glue(
                                "Pop0: {nind[1]} diploid individuals.Columns {pop0_cols[1]}:{pop0_cols[length(pop0_cols)]}. This is your population of interest."
                        ),
                        "\n"
                )
        }

        x <- data.table::setDT(x)
        tableout <-
                x %>% dplyr::select(CHR, POS, REF, ALT)
        tableout <- dplyr::bind_cols(
                tableout,
                x %>%
                        dplyr::rowwise() %>%
                        dplyr::reframe(across(all_of(pop0_cols), .count_alleles)) %>%
                        dplyr::reframe(tx_1 = Reduce(`+`, .)),
                x %>%
                        dplyr::rowwise() %>%
                        dplyr::reframe(across(all_of(pop0_cols), .count_nonmissing)) %>%
                        dplyr::reframe(tn_1 = Reduce(`+`, .))
        ) %>%
                dplyr::select(CHR, POS, REF, ALT, tx_1, tn_1) %>%
                as.data.table
        return(tableout)
}
