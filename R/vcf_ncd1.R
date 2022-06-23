vcf_ncd1<-function(x=inp,out="outfile.txt", nind=c(108,10)){
 
  tableout<-data.table(CHR=NA,POS=NA,REF=NA,ALT=NA,tx_1=NA,tn_1=NA,tx_2=NA,tn_2=NA)  
  for(l in 1:nrow(inp)){
    chr<-inp[l,1]
    pos<-inp[l,2]
    ref<-inp[l,4]
    alt<-inp[l,5]
    if(ref %in% c("A","C","G","T") & alt %in% c("A","C","G","T")){
      anc<-inp[l,get(colnames(inp)[index.col+nind[1]+1])]
      anc<-str_split(anc,"|",simplify=F)[[1]][[2]]
      drv=0 ; total=0
      for(i in index.col:(index.col+nind[1])){
        al<-str_split(inp[l,get(colnames(inp)[i])],"|",simplify=F)[[1]]
        if(al[1]!=anc){drv=drv+1}
        if(al[2]!=anc){drv=drv+1}
        total = total + 2
      }
    }
    cat(l,"\r")
    tableout<-rbind(tableout, data.table(chr,pos,ref,alt,tx_1=drv,tn_1=total,tx_2=(nind[2]*2),tn_2=(nind[2]*2)))
    }
  tableout[,MAF:=ifelse((tx_1/tn_1)>0.5, 1-(tx_1/tn_1), (tx_1/tn_1))]
  tableout<-tableout[-1,]
  sink(out) 
  print(tableout)
  sink()