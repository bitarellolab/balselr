#' Split a genotype column from a vcf file
#'
#' @param x A string from a GT column.
#' @param split Pattern to use for splitting.
#' @param index Which elements to keep 1, 2, or 1 & 2
#'
#' @examples split_geno(x = c("0|0","0|1","1|1"))
#' data.table(col1="1|1", col2="1|0", col3="0|0") %>% dplyr::summarise(across(col1:col3, .split_geno))
split_geno <- function(x, split="|",index=c(1,2)) {
        #assertthat::assert_that(is.numeric(which.al), err = "Parameter which.al must be numeric")

        if(length(x) == 1){
        as.vector(as.integer(stringr::str_split(string = x, pattern = split, simplify = T)[,c(2,4)][index]))
        }else{
        as.vector(as.integer(stringr::str_split(string = x, pattern = split, simplify = T)[,c(2,4)][,index]))
}
}
