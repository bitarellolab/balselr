#' Split a genotype column from a vcf file
#'
#' @param x
#' @param split
#'
#' @return
#' @export
#'
#' @examples
split_geno <- function(x = x, split = "|") {
  res <- as.integer(stringr::str_split(string = x, pattern = "|", simplify = F)[[1]][c(2, 4)])
  return(res)
}
