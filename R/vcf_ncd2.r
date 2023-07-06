#' Make input file for NCD2 from VCF
#
#' @param x A data.table object generated from reading a VCF file with read_vcf.
#' @param outfile The path and name for the outfile. If not provided,
#' this will be  a timestamp file in the current directory.
#' @param nind A vector containing the number of diploid individuals from each
#' population to be used. The length of the vector should be equal to the number
#' of populations being sampled.
#' @param index.col First genotype column in VCF file.

.vcf_ncd2 <-
        function(x,
                 outfile = outfile,
                 nind = nind,
                 index.col = index.col) {
                #
                #CHR <- POS <- REF <- ALT <- tx_1 <- tn_1 <- NULL
                # tx_2 <- tn_2 <- NULL
                #tictoc::tic("Total runtime")
                npop <- length(nind)
                assertthat::assert_that(npop == 2 |
                                                npop == 3, msg = "If only one species is represented, you should use ncd1.\n")
                pop0_cols <- colnames(x)[index.col:(index.col + (nind[1] - 1))]
                if (npop > 2) {
                        nind <- nind[c(1, 2)]
                }
                npop2 <- length(nind)
                #if (verbose == T) {
                cat(
                        glue::glue(
                                "Your vcf contains {npop} populations. Only {npop2} are considered:"
                        ),
                        "\n"
                )
                #}
                # if (verbose == T) {
                cat(
                        glue::glue(
                                "Pop0: {nind[1]} diploid individuals.Columns {pop0_cols[1]}:{pop0_cols[length(pop0_cols)]}. This is your population of interest."
                        ),
                        "\n"
                )
                # }
                pop1_cols <-
                        colnames(x)[(index.col + nind[1]):(index.col + (sum(nind) - 1))]
                #if (verbose == T) {
                cat(
                        glue::glue(
                                "Pop1: {nind[2]} diploid individuals.Columns {pop1_cols[1]}:{pop1_cols[length(pop1_cols)]}. This is your outgroup."
                        ),
                        "\n"
                )
                # }

               # counter <- 0
               # if (verbose == T) {
               #         cat("Parsing vcf lines....\n")
               # }
                #for (l in 1:nrow(x)) {
                #if(verbose==T){cat(l, "\r")}


                x <- data.table::setDT(x)
                tableout <-
                        x %>% dplyr::select(CHR, POS, REF, ALT) %>% as.data.table
                tableout<-dplyr::bind_cols(tableout,
                                           x %>%
                                                   #dplyr::select(all_of(pop0_cols)) %>%
                                                   dplyr::rowwise() %>%
                                                   dplyr::reframe(across(pop0_cols, .count_alleles)) %>%
                                                   dplyr::reframe(tx_1 = Reduce(`+`,.)),
                                           x %>%
                                                   #dplyr::select(all_of(pop0_cols)) %>%
                                                   dplyr::rowwise() %>%
                                                   dplyr::reframe(across(pop0_cols, .count_nonmissing)) %>%
                                                   dplyr::reframe(tn_1 = Reduce(`+`,.)),
                                           x %>%
                                                   #dplyr::select(all_of(pop1_cols)) %>%
                                                   dplyr::rowwise() %>%
                                                   dplyr::reframe(across(pop1_cols, .count_alleles)) %>%
                                                   dplyr::reframe(tx_2 = Reduce(`+`,.)),
                                           x %>%
                                                   #dplyr::select(all_of(pop1_cols)) %>%
                                                   dplyr::rowwise() %>%
                                                   dplyr::reframe(across(pop1_cols, .count_nonmissing)) %>%
                                                   dplyr::reframe(tn_2 = Reduce(`+`,.))
                ) %>%
                       # dplyr::mutate(tn_1 = nind[1] * 2, tn_2 = nind[2] *
                        #                            2) %>%
                        dplyr::select(CHR, POS, REF, ALT, tx_1, tn_1, tx_2, tn_2) %>%
                        as.data.table

                #if (verbose == T) {
                cat(glue::glue("Printing output to {outfile}..."),
                    "\n")
                #}
                if(!is.null(outfile)){
                data.table::fwrite(tableout,
                                   file = outfile,
                                   sep = "\t",
                                   col.names = T)
                }
                #if(verbose==T){tictoc::toc()}
                return(tableout)
        }
