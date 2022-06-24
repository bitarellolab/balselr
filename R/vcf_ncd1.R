# @
#
# This is a function called vcf_ncd1.R
# which makes an input file for NCD1 from a vcf file
#
# You can learn more about package authoring with RStudio at:
#
#   v
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
#' @export
vcf_ncd1<-function(x=inp,outfile="outfile.txt", nind=nind, index.col=index.col, fold=T){
        assertthat::assert_that(length(nind)==2|length(nind)==3)
        if(fold==T){
                nind<-nind[c(1,2)]
        }
        npop<-length(nind)
        cat(glue::glue("Parsing data from {npop} populations..."),"\n")
        tableout<-data.table::data.table(CHR=NA,POS=NA,ANC=NA,REF=NA,ALT=NA,tx_1=NA,tn_1=NA,tx_2=NA,tn_2=NA)
        counter=0
        for(l in 1:nrow(x)){
                cat(l,"\r")
                chr<-x[l,1]
                pos<-x[l,2]
                ref<-x[l,4]
                alt<-x[l,5]

                if(!(ref %in% c("A","C","G","T") & alt %in% c("A","C","G","T"))){
                        cat(glue::glue("Line {l} {ref} {alt} will be skipped. Only bi-allelic SNPs PASS"),"\n")
                        counter=counter+1
                        next
                }
                x<-data.table::setDT(x)
                if(fold==T){
                anc<-x[l,colnames(x)[index.col+nind[1]+1], with=F]
                }else{
                anc<-x[l,colnames(x)[index.col+sum(nind)], with=F]
                }
                anc<-stringr::str_split(anc,"|",simplify=F)[[1]][[2]]
                drv=0 ; total=0;drv2=0; total2=0;
                for(i in index.col:(index.col+(nind[1]-1))){
                        al<-stringr::str_split(x[l,colnames(x)[i], with=F],"|",simplify=F)[[1]]
                        if(al[1]!=anc){drv=drv+1}
                        if(al[2]!=anc){drv=drv+1}
                        total = total + 2
                }
                for(i in (index.col+nind[1]):(index.col+nind[1]+(nind[2]-1))){
                        al<-stringr::str_split(x[l,colnames(x)[i], with=F],"|",simplify=F)[[1]]
                        if(al[1]!=anc){drv2=drv2+1}
                        if(al[2]!=anc){drv2=drv2+1}
                        total2 = total2 + 2
                }
                tableout<-rbind(tableout, data.table::data.table(chr,pos,anc,ref,alt,tx_1=drv,tn_1=total,tx_2=drv2,tn_2=total2),use.names=FALSE)
        }
        tableout<-tableout[-1,]
        cat(glue::glue("Printing output to {outfile}..."),"\n")
        data.table::fwrite(tableout, file=outfile, sep="\t", col.names=F)
        cat(glue::glue("Lines skipped: {counter}"),"\n")
}
#install.packages(c("devtools", "roxygen2", "testthat", "knitr"))
