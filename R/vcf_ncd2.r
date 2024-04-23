#' Make input file for NCD2 from VCF
#
#' @param x A data.table object generated from reading a VCF file with read_vcf.
#' @param nind A vector containing the number of diploid individuals from each
#' population to be used. The length of the vector should be equal to the number
#' of populations being sampled.
#' @param index.col First genotype column in VCF file.
#' @param verbose Logical. Printed updates. Default is T.

#to do: some of the tests suggested in parse_vcf might pertain here.
#to do: check that parameters for the .vcf_ncd2 and .vcf_ncd1 are coded and described similarly.
#check that for the example command in parse_vcf using ncd2, all values in column tn_2 are 2 and all values in tx_2 are between 0 and 2.Do a similar test for the ncd1 case.

.vcf_ncd2 <-
        function(x,
                 nind = nind,
                 index.col = index.col,
                 verbose = T) {
                npop <- length(nind)
                assertthat::assert_that(npop == 2 |
                                                npop == 3, msg = "If only one species is represented, you should use ncd1.\n")
                pop0_cols <-
                        colnames(x)[index.col:(index.col + (nind[1] - 1))]
                if (npop > 2) {
                        nind <- nind[c(1, 2)]
                        if (verbose == T) {
                                cat(
                                        glue::glue(
                                                "Only the first two elements of nind are considered:"
                                        ),
                                        "\n"
                                )
                        }
                }
                npop2 <- length(nind)
                if (verbose == T) {
                        cat(
                                glue::glue(
                                        "Your vcf contains {npop} populations. Only {npop2} are considered:"
                                ),
                                "\n"
                        )
                }
                if (verbose == T) {
                        cat(
                                glue::glue(
                                        "Pop0: {nind[1]} diploid individuals.Columns {pop0_cols[1]}:{pop0_cols[length(pop0_cols)]}. This is your population of interest."
                                ),
                                "\n"
                        )
                }
                npop <- npop2
                pop1_cols <-
                        colnames(x)[(index.col + nind[1]):(index.col + (sum(nind) - 1))]
                if (verbose == T) {
                        cat(
                                glue::glue(
                                        "Pop1: {nind[2]} diploid individuals.Columns {pop1_cols[1]}:{pop1_cols[length(pop1_cols)]}. This is your outgroup."
                                ),
                                "\n"
                        )
                }

                x <- data.table::setDT(x)
                tableout <-
                        x %>% dplyr::select(CHR, POS, REF, ALT) %>% as.data.table
                tableout <- dplyr::bind_cols(
                        tableout,
                        x %>%
                                dplyr::rowwise() %>%
                                dplyr::reframe(across(all_of(pop0_cols), .count_alleles)) %>%
                                dplyr::reframe(tx_1 = Reduce(`+`, .)),
                        x %>%
                                dplyr::rowwise() %>%
                                dplyr::reframe(across(all_of(pop0_cols), .count_nonmissing)) %>%
                                dplyr::reframe(tn_1 = Reduce(`+`, .)),
                        x %>%
                                dplyr::rowwise() %>%
                                dplyr::reframe(across(all_of(pop1_cols), .count_alleles)) %>%
                                dplyr::reframe(tx_2 = Reduce(`+`, .)),
                        x %>%
                                dplyr::rowwise() %>%
                                dplyr::reframe(across(all_of(pop1_cols), .count_nonmissing)) %>%
                                dplyr::reframe(tn_2 = Reduce(`+`, .))
                ) %>%
                        dplyr::select(CHR, POS, REF, ALT, tx_1, tn_1, tx_2, tn_2) %>%
                        as.data.table

                return(tableout)
        }
