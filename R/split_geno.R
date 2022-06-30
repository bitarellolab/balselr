#' Split a genotype column from a vcf file
#'
#' @param x A string from a GT column.
#' @param split Pattern to use for splitting.
#' @param n Which elements to keep
#' @export
#'
#' @examples split_geno("0|1")
split_geno <- function(x=x, split="|",n=c(1,2)) {
  as.integer(stringr::str_split(string = x, pattern = split, simplify = F)[[1]][c(2,4)][n])
}
