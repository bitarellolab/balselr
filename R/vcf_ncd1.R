#' Make input file for NCD1 from VCF
#'
#' @param x An object of class data.table containing
#' @param outfile The path and name for the outfile. If not provided,
#' this will be  a timestamp file in the current directory.
#' @param nind A vector containing the number of diploid individuals from each
#' population to be used. For ncd1, only one population is used.
#' @param index.col First genotype column in VCF file.
#' @param verbose  Logical. If TRUE, progress reports will be printed as the
#' function runs.
#'
#'
#' @export
#' @examples inp = read_vcf("inst/example.vcf")
#' vcf_ncd1(x=inp, outfile=outfile_path("inst/example.vcf"),nind=c(108),
#' index.col=10, verbose=T)
vcf_ncd1 <- function(x = infile,
                     outfile = outfile,
                     nind = nind,
                     index.col = index.col,
                     verbose = verbose) {
  #
  npop <- length(nind)
  assertthat::assert_that(npop == 1, msg = "NCD1 only uses one species.\n")
  pop0_cols <- c(index.col, index.col + (nind[1] - 1))

  if (verbose == T) {
    cat(
      glue::glue("Only {npop} population will be considered:"), "\n"
    )
  }
  if (verbose == T) {
    cat(
      glue::glue("Pop0: {nind[1]} diploid individuals.Columns {pop0_cols[1]}:{pop0_cols[2]}. This is your population of interest."), "\n"
    )
  }

  tableout <- data.table::data.table(
    CHR = NA,
    POS = NA,
    REF = NA,
    ALT = NA,
    tx_1 = NA,
    tn_1 = NA
  )
  counter <- 0
  for (l in 1:nrow(x)) {
    cat(l, "\r")
    chr <- x[l, 1]
    pos <- x[l, 2]
    ref <- x[l, 4]
    alt <- x[l, 5]
    x <- data.table::setDT(x)

    total <- 0
    alt2 <- 0
    for (i in seq(from = pop0_cols[1], to = pop0_cols[2])) {
      al <-
        split_geno(x = x[l, colnames(x)[i], with = F], split = "|")
      if (al[1] == 1) {
        alt2 <- alt2 + 1
      }
      if (al[2] == 1) {
        alt2 <- alt2 + 1
      }
      total <- total + 2
    }
    tableout <-
      rbind(
        tableout,
        data.table::data.table(
          chr,
          pos,
          ref,
          alt,
          tx_1 = alt2,
          tn_1 = total
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
    col.names = F
  )

  return(tableout)
}
