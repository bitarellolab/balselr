#' Make input file for NCD2 from VCF
#
#' @param x A data.table object generated from reading a VCF file with read_vcf.
#' @param outfile The path and name for the outfile. If not provided,
#' this will be  a timestamp file in the current directory.
#' @param nind A vector containing the number of diploid individuals from each
#' population to be used. The length of the vector should be equal to the number
#' of populations being sampled.
#' @param index.col First genotype column in VCF file.
#' @param fold Logical. If TRUE, the output will have alternate allele counts.
#' If FALSE, derived allele counts will be used. Default is TRUE.
#' @param verbose Logical. If TRUE, progress reports will be printed as the
#' function runs.
#'
#' @examples
#' vcf_ncd2(x=inp, outfile=outfile_path("inst/example.vcf"),nind=c(108,10),
#' index.col=10, fold=T, verbose=T)
#'
.vcf_ncd2 <-
  function(x = inp,
           outfile = outfile,
           nind = nind,
           index.col = index.col,
           fold = fold,
           verbose = verbose) {
          #
          tictoc::tic("Total runtime")
    npop <- length(nind)
    assertthat::assert_that(npop == 2 |
    npop == 3, msg = "If only one species is represented, you should use ncd1.\n")
    pop0_cols <- c(index.col, index.col + (nind[1] - 1))
    if (npop > 2) {
      nind <- nind[c(1, 2)]
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
          "Pop0: {nind[1]} diploid individuals.Columns {pop0_cols[1]}:{pop0_cols[2]}. This is your population of interest."
        ),
        "\n"
      )
    }
    pop1_cols <-
      c(index.col + nind[1], index.col + (sum(nind) - 1))
    if (verbose == T) {
      cat(
        glue::glue(
          "Pop1: {nind[2]} diploid individuals.Columns {pop1_cols[1]}:{pop1_cols[2]}. This is your outgroup."
        ),
        "\n"
      )
    }
    tableout <-
      data.table::data.table(
        CHR = NA,
        POS = NA,
        ANC = NA,
        REF = NA,
        ALT = NA,
        tx_1 = NA,
        tn_1 = NA,
        tx_2 = NA,
        tn_2 = NA
      )
    counter <- 0
    if(verbose==T){cat("Parsing vcf lines....\n")}
    for (l in 1:nrow(x)) {
    if(verbose==T){cat(l, "\r")}
      chr <- x[l, 1]
      pos <- x[l, 2]
      ref <- x[l, 4]
      alt <- x[l, 5]


      x <- data.table::setDT(x)
      anc <-
        x[l, colnames(x)[index.col + (sum(nind) - 1)], with = F]
      anc <-
        split_geno(x = anc, split = "|")[1]
      drv <- 0
      total <- 0
      drv2 <- 0
      total2 <- 0

      for (i in seq(from = pop0_cols[1], to = pop0_cols[2])) {
        al <-
          split_geno(x = x[l, colnames(x)[i], with = F], split = "|")
        if (al[1] != anc) {
          drv <- drv + 1
        }
        if (al[2] != anc) {
          drv <- drv + 1
        }
        total <- total + 2
      }
      for (i in seq(from = pop1_cols[1], to = pop1_cols[2])) {
        al <-
          split_geno(x = x[l, colnames(x)[i], with = F], split = "|")
        if (al[1] != anc) {
          drv2 <- drv2 + 1
        }
        if (al[2] != anc) {
          drv2 <- drv2 + 1
        }
        total2 <- total2 + 2
      }
      tableout <-
        rbind(
          tableout,
          data.table::data.table(
            chr,
            pos,
            ref,
            ref,
            alt,
            tx_1 = drv,
            tn_1 = total,
            tx_2 = drv2,
            tn_2 = total2
          ),
          use.names = FALSE
        )
    }
    tableout <- tableout[-1, ]
    if (verbose == T) {
      cat(
        glue::glue("Printing output to {outfile}..."),
        "\n"
      )
    }
    data.table::fwrite(tableout,
      file = outfile,
      sep = "\t",
      col.names = T
    )
    if(verbose==T){tictoc::toc()}
    return(tableout)
  }
