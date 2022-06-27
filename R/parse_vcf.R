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
parse_vcf_slim <-
        function(infile = "*.vcf",
                 outfile = NA,
                 index.col = NA,
                 nind = c(n0, n2, n1),
                 type = "ncd2",
                 fold = NA,
                 intern = T,
                 verbose = T) {
                assertthat::assert_that(file.exists(infile),
                                        msg = glue::glue("VCF file {infile} does not exist.\n"))
                pref = gsub(".vcf", "", infile)
                inp <- data.table::fread(infile, skip = "##", header = T)
                data.table::setnames(inp, "#CHROM", "CHR")
                type <- tolower(type)
                if (is.na(index.col)) {
                        index.col <- which(colnames(inp) == "FORMAT") + 1
                        if (verbose == T) {
                                cat(
                                        glue::glue(
                                                "index.col not provided so we will assume it is column {index.col}..."
                                        ),
                                        "\n"
                                )
                        }
                }
                if (verbose == T) {
                        cat(glue::glue("Creating input file for {type} from {infile}..."),
                            "\n")
                }
                outfile <- glue::glue("{outfile}_{type}.out")
                #Index No. of the individual to use as ``ancestral'' sequence
                if (fold == F & type == "ncd2") {
                        outseq <- (index.col) + (sum(nind) - 1)
                        if (verbose == T) {
                                cat(
                                        glue::glue(
                                                "You asked for the unfolded version. You need a third species. I will use column {outseq} for this."
                                        ),
                                        "\n"
                                )
                        }
                } else if (type == "ncd1") {
                        assertthat::assert_that(fold == T, msg = "Only the folded option is available when Ancestral/Derived states are unknown.\n")
                }
                if (is.na(outfile)) {
                        outfile = tempfile("file", fileext = type)
                        if (verbose == T) {
                                cat(
                                        glue::glue(
                                                "No outfile provided. Will write this into tmp file {outfile}"
                                        ),
                                        "\n"
                                )
                        }
                }
                #
#
if (type == "ncd2") {
        res <-
                vcf_ncd2(
                        x = inp,
                        outfile = outfile,
                        nind = nind,
                        index.col = index.col,
                        fold = fold,
                        verbose = verbose
                )
} else if (type == "ncd1") {
        res <-
                vcf_ncd1(
                        x = inp,
                        outfile = outfile,
                        nind = nind,
                        index.col = index.col,
                        verbose = verbose
                )
}
if (verbose == T) {
        cat(glue::glue("Finished making input file for {type}"),
            "\n")
}else{
        cat(glue::glue("File written to: {outfile}"),"\n")
}

if (intern == T) {
        return(res)
}
}
#
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
