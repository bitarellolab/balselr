now <- Sys.time()
timestamp <- function(time) format(time, "%Y-%B-%d_%H-%M-%S")

#' Create output file with timestamp
#' @param infile Name of input file to modify
#' @example
#' @export outfile_path("inst/example.vcf")
outfile_path <- function(infile) {

        sub("\\..*", paste0("_",timestamp(now),".out"), infile)
}
