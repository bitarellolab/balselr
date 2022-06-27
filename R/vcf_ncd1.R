#' Make input file for NCD1 from VCF
#'
#' @param x
#' @param outfile
#' @param nind
#' @param index.col
#' @param verbose
#'
#' @export
#'
#' @examples
vcf_ncd1 <-
        function(x = inp,
                 outfile = outfile,
                 nind = nind,
                 index.col = index.col,
                 verbose = verbose) {
                assertthat::assert_that(length(nind) == 1, msg = "Only one species should be used with ncd1.\n")
                pop0_cols <- c(index.col, index.col + nind-1)
                cat(glue::glue("Your vcf contains 1 population:"), "\n")
                if (verbose == T) {
                        cat(
                                glue::glue(
                                        "Pop0: {nind} diploid individuals.Columns {pop0_cols[1]}:{pop0_cols[2]}. If you have an outgroup species, use ncd2 instead."
                                ),
                                "\n"
                        )
                }
                tableout <-
                        data.table::data.table(
                                CHR = NA,
                                POS = NA,
                                REF = NA,
                                ALT = NA,
                                tx_1 = NA,
                                tn_1 = NA
                        )

                counter = 0
                for (l in 1:nrow(x)) {
                        cat(l, "\r")
                        chr <- x[l, 1]
                        pos <- x[l, 2]
                        ref <- x[l, 4]
                        alt <- x[l, 5]

                        if (!(ref %in% c("A", "C", "G", "T") &
                              alt %in% c("A", "C", "G", "T"))) {
                                if (verbose == T) {
                                        cat(
                                                glue::glue(
                                                        "Line {l} {ref} {alt} will be skipped. Only bi-allelic SNPs PASS"
                                                ),
                                                "\n"
                                        )
                                }
                                counter = counter + 1
                                next
                        }
                        x <- data.table::setDT(x)
                        anc = 0
                        drv = 0
                        total = 0
                        for (i in seq(from=pop0_cols[1],to=pop0_cols[2])) {
                                al <-
                                        split_geno(x=x[l, colnames(x)[i], with = F], split="|")
                                if (al[1] != anc) {
                                        drv = drv + 1
                                }
                                if (al[2] != anc) {
                                        drv = drv + 1
                                }
                                total = total + 2
                        }

                        tableout <-
                                rbind(
                                        tableout,
                                        data.table::data.table(
                                                chr,
                                                pos,
                                                ref,
                                                alt,
                                                tx_1 = drv,
                                                tn_1 = total
                                        ),
                                        use.names = FALSE
                                )
                }
                tableout <- tableout[-1,]
                if (verbose == T) {
                        cat(glue::glue("Printing output to {outfile}..."),
                            "\n")
                }
                data.table::fwrite(tableout,
                                   file = outfile,
                                   sep = "\t",
                                   col.names = F)
                if (verbose == T) {
                        cat(glue::glue("Lines skipped: {counter}"), "\n")
                }

                return(tableout)

        }
