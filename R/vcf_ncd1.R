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
vcf_ncd1<-function(x=inp,outfile="outfile.txt", nind=NA, index.col=index.col){
        npop<-length(nind)
        tableout<-data.table::data.table(CHR=NA,POS=NA,REF=NA,ALT=NA,tx_1=NA,tn_1=NA,tx_2=NA,tn_2=NA)
        cat(glue::glue("Parsing data from {npop} populations..."),"\n")
        for(l in 1:nrow(x)){
                cat(l,"\n")
                chr<-x[l,1]
                pos<-x[l,2]
                ref<-x[l,4]
                alt<-x[l,5]
                #print(ref)
                if(!(ref %in% c("A","C","G","T") & alt %in% c("A","C","G","T"))){
                cat("Only bi-allelic SNPs...PASS\n")
                next
                }
                x<-data.table::setDT(x)
                #print(index.col)
                print(colnames(x)[index.col+nind[1]+1])
                anc<-x[l,colnames(x)[index.col+nind[1]+1], with=F]
                anc<-stringr::str_split(anc,"|",simplify=F)[[1]][[2]]
                drv=0 ; total=0
                for(i in index.col:(index.col+nind[1])){
                        print(i)
                        al<-stringr::str_split(x[l,colnames(x)[i], with=F],"|",simplify=F)[[1]]
                        if(al[1]!=anc){drv=drv+1}
                        if(al[2]!=anc){drv=drv+1}
                        total = total + 2
                }
                total2=(nind[2]*2)
                tableout<-rbind(tableout, data.table::data.table(chr,pos,ref,alt,tx_1=drv,tn_1=total,tx_2=total2,tn_2=total2),use.names=FALSE)
               # }
        }
        tableout<-tableout[-1,]
        cat(glue::glue("Printing output to {outfile}..."),"\n")
        sink(outfile)
        print(tableout)
        sink()
        cat("Done.","\n")
}
#install.packages(c("devtools", "roxygen2", "testthat", "knitr"))
