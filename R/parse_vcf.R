# @
#
# This is a function called parse_vcf.
# which makes an input file for NCD from a vcf file
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
parse_vcf_slim<-function(dir=".",infile="*.vcf", outfile=NA,ncores=2, index.col=NA, nind=c(108,10,1), outseq=NA, type="ncd1"){
  #
  assertthat::assert_that(dir.exists(dir), msg=glue::glue("Dir {dir} does not exist.\n"))
  assertthat::assert_that(file.exists(glue::glue("{dir}{infile}")), msg=glue::glue("VCF file {infile} does not exist.\n"))
  tmp<-glue::glue("{dir}{infile}")
  pref=gsub(".vcf","",tmp)
  inp<-data.table::fread(tmp, skip="##", header=T)
  data.table::setnames(inp, "#CHROM","CHR")
  if(is.na(index.col)){

  index.col<-which(colnames(inp)=="FORMAT")+1
  cat(glue::glue("index.col not provided so we will assume it is column {index.col}..."),"\n")
  }
  #Index No. of the individual to use as ``ancestral'' sequence
  if(is.na(outseq) & !(type %in% c("ncd1","ncd2"))){
    outseq<-(index.col-1)+sum(nind)
  }else if(type %in% c("ncd1","ncd2")){
    remove(outseq)
  }
  if(is.na(outfile)){
    outfile=gsub(".vcf","_",tmp)
  }else{
    outfile=glue("{dir}{outfile}_")
  }
  cat(glue::glue("output will be written into file {outfile}<name of test>"),"\n")
  #
  if(type=="ncd1"){
          cat(glue::glue("Creating input file for ncd1 from {infile}..."),"\n")
          vcf_ncd1(x=inp,outfile=glue::glue("{outfile}ncd1.out"), nind=nind[c(1,2)], index.col=index.col)
          cat("Done.\n")
  }
}
 # out.ncd1.unf<-glue("{outfile}ncd1_unf.out")
 # out.ncd2<-glue("{outfile}ncd2")
 # out.ncd2.unf<-glue("{outfile}ncd2_unf.out")
 # out.betascan<-glue("{outfile}betascan.out")
 #out.betascan.unf<-glue("{outfile}betascan_unf.out")
 # out.mutebass<-glue("{outfile}mutebass.out")
 # out.balmixder<-glue("{outfile}mutebass.out")
 # out.ballet.snp<-glue("{outfile}ballet_snp.out")
 # out.ballet.rec<-glue("{outfile}ballet_rec.out")
  #
#	for(l in 1:nrow(inp)){
	  #chr<-inp[l,1]
	  #pos<-inp[l,2]
	  #ref<-inp[l,4]
	  #alt<-inp[l,5]
	  #anc<-inp[l,get(colnames(inp)[outseq])]
	  #anc<-str_split(anc,"|",simplify=F)[[1]][[2]]
	  #drv=0 ; total=0; drv2=0; total2=0
  #for(i in index.col:(index.col+nind[1])){
  #    al<-str_split(inp[l,get(colnames(inp)[i])],"|",simplify=F)[[1]]
  #    if(al[1]!=anc){drv=drv+1}
  #      if(al[2]!=anc){drv=drv+1}
#		      total = total + 2
#		}
 # for(i in ((index.col+nind[1])+1):(index.col+nind[2])){
  #    al<-str_split(inp[l,get(colnames(inp)[i])],"|",simplify=F)[[1]]
   #   if(al[1]!=anc){drv2=drv2+1}
    #    if(al[2]!=anc){drv2=drv2+1}
     #     total2 = total2 + 2
     # }
#	sink(out.ncd1)
#  cat("CHR\tPOS\tREF\tALT\tx_1\tn_1\tx_2\tn_2")
#  cat(chr,pos,ref,alt,drv,total, drv2, total2)
#	sink()
#		print $1, $2, 1.2e-8*$2,drv,total, drv2, total2 > outmute # recombination rate
#    		print $2, drv,total > betaout
#                if(drv>0) print $2, total-drv, total > balout".snp"
#                if(drv>0){
#                    print $2, ($2-prevpos)*2.75e-4 > balout".rec" # pop-scaled recombination rate (2Nr = 2*11477*1.2e-8 = 2.75448e-4)
#                    prevpos=$2
#                }
#    		if(drv>0) print $2, "NA", drv, total > balmixder #BalLeRMix_manual p.5: When the user does not have recombination maps for reference, the seceond column
#can be NAs as long as the user makes sure to use --physPos.
	#}
#}'

#	sink(out.mutebass)
#	cat("CHR\tPOS\tGenPOS\tx_1\tn_1\tx_2\tn_2")
#	sink()
#	sink(out.balmixder)
#	cat("physPos\tgenPos\tx\tn")
#	sync()
#	sync(out.ballet.snp)
 # cat("position\tx\tn")
  #sync()
  #sync(out.ballet.rec)
  #cat("position\trate")
  #sync()
  #}
  #outseq=128
#}
