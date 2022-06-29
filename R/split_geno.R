#' Split a genotype column from a vcf file
#'
#' @param x A string from a GT column.
#' @param split Pattern to use for splitting.
#'
#' @export
#'
#' @examples split_geno("0|1")
split_geno <- function(...) {
  as.integer(stringr::str_split(string = x, pattern = "|", simplify = F)[[1]][c(2, 4)])
}
