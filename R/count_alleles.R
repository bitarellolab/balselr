#' Count derived or ancestral alleles for a given vcf position
#' @param x A vector containing 0s and 1s corresponding to allele states
#' @param der Return derived. Logical. If TRUE, returns the counts of
#' derived (alternate)
#' alleles. If FALSE, returns ancestral (ref) counts
#' @param split Pattern to look for. Inherited from stringr::str_split
#' @param index 1 or 2 or c(1,2). Which alleles to consider.
#' @examples data.table::data.table(col1="1|1", col2="1|0", col3="0|0") %>%
#' dplyr::summarise(across(col1:col3, .count_alleles))
.count_alleles<-function(x, split="|",index=c(1,2), der = F){
        x<-.split_geno(x = x, split = split, index = index)
        if(der == T){sum(x, na.rm=T)}else if (der == F){sum(!is.na(x))-sum(x, na.rm=T)}
}
.count_nonmissing<-function(x, split="|",index=c(1,2)){
        x<-.split_geno(x = x, split = split, index = index)
        sum(!is.na(x))
}
