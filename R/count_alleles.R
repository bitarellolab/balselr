#' Count derived or ancestral alleles for a given vcf position
#' @param x A vector containing 0s and 1s corresponding to allele states. 0/0: homozygous reference; 0/1: heterozygous; 1/1: homozygous alternate

#' @param split Pattern to look for. Inherited from stringr::str_split


.count_alleles<-function(x, split="|"){
        x<-.split_geno(x = x, split = split)
        #if(der == T){sum(x, na.rm=T)}else if (der == F){sum(!is.na(x))-sum(x, na.rm=T)}
        #sum(!is.na(x))-sum(x, na.rm=T)
        ref<-sum(x==0, na.rm=T) #number of ref alleles
        alt<-sum(x==1, na.rm=T) #number of alt alleles

}
.count_nonmissing<-function(x, split="|"){
        x<-.split_geno(x = x, split = split)
        sum(!is.na(x))
}
