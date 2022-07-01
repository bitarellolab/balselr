#' Count derived or ancestral alleles for a given vcf position
#' @param x A vector containing 0s and 1s corresponding to allele states
#' @param der Return derived. Logical. If TRUE, returns the counts of derived (alternate)
#' alleles. If FALSE, returns ancestral (ref) counts
#' @examples data.table(col1="1|1", col2="1|0", col3="0|0") %>% dplyr::summarise(across(col1:col3, .count_alleles))
.count_alleles<-function(x, split="|",index=c(1,2), der = T){
        x<-.split_geno(x = x, split = split, index = index)

        if(der == T){sum(x)}else if (der == F){length(x)-sum(x)}
}
