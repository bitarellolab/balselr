#' Split a genotype column from a vcf file
#'
#' @param x A string from a GT column.
#' @param split Pattern to use for splitting.

.split_geno <- function(x, split="|") {
        #assertthat::assert_that(is.numeric(which.al), err = "Parameter which.al must be numeric")

       # if(length(x) == 1){
        #as.vector(as.integer(stringr::str_split(string = x, pattern = split, simplify = T)[,c(2,4)]))
        suppressWarnings(as.numeric(unlist(lapply(stringr::str_split(string = x, pattern = "|"), function(x) x[c(2,4)]))))
       # }else{

#}
}
