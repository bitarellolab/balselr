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

parse_vcf<-function(dir="example/",infile="example.vcf", outfile=NA, outseq=NA,ncores=2){
  require(tictoc)
  require(glue)
  require(testthat)
  require(data.table)
  tic()
  
  tmp<-glue("{dir}{infile}")
  if(is.na(outfile)){
    out=gsub(".vcf",".out",tmp)
  }
  else{
    out=glue("{dir}{outfile}"
  }
  pref=gsub(".vcf","",tmp)
  inp<-fread("tmp")
  #Index No. of the individual to use as ``ancestral'' sequence
  if(is.na(outseq)){

    
  }
  outseq=128
  
  
}