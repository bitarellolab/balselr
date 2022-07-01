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
#' @examples inp = read_vcf("inst/example.vcf")
#' vcf_ncd1(x=inp, outfile=outfile_path("inst/example.vcf"),nind=c(108),
#' index.col=10, verbose=T)
.vcf_ncd1 <- function(x = infile,
                     outfile = outfile,
                     nind = nind,
                     index.col = index.col,
                     verbose = verbose) {
  #
        tictoc::tic("Total runtime")
  npop <- length(nind)
  assertthat::assert_that(npop == 1, msg = "NCD1 only uses one species.\n")
  pop0_cols <- colnames(x)[c(index.col: index.col + (nind[1] - 1))]

  if (verbose == T) {
    cat(
      glue::glue("Only {npop} population will be considered:"), "\n"
    )
  }
  if (verbose == T) {
    cat(
      glue::glue("Pop0: {nind[1]} diploid individuals.Columns {pop0_cols[1]}:{pop0_cols[length(pop0_cols)]}. This is your population of interest."), "\n"
    )
  }

  x <- data.table::setDT(x)
  tableout<-x %>% dplyr::select(CHR, POS, REF, ALT) %>% as.data.table
  tableout<-dplyr::bind_cols(tableout,x %>% dplyr::select(pop0_cols) %>%
                dplyr::rowwise() %>%
                dplyr::summarise(across(pop0_cols, .count_alleles)) %>%
                dplyr::summarise(tx_1 = Reduce(`+`,.))) %>%
                dplyr::mutate(tn_1 = nind[1]*2) %>%
                dplyr::select(CHR, POS, REF, ALT, tx_1, tn_1) %>%
                as.data.table
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
if(verbose==T){tictoc::toc()}
  return(tableout)
}
