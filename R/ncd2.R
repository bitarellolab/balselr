#' Calculate non-central statistic (NCD2)
#'
#' @param x A data.table object with
#' columns: CHR, POS, REF, ALT, tx_1 (number of alternate allele copies), tn_1 (total number of alleles), tx_2, and tn_2.
#' @param tf Target frequency. Any value between 0 and 0.5. Default is 0.5.
#' @param w Window size in bp. Default is 1000.
#' @param ncores Number of cores. Increasing this can speed things up for you.
#' Default is 2.
#' @param minIS Minimum number of informative sites. Default is 2. Windows with
#'  less informative sites than this threshold are discarded.
#' @param by Define how to scan the genome. "POS" (default) defined sliding windows based on w. "IS" defined windows around each informative site.
#' @return A data.table object with columns:  Win.ID, S (sites), FD (fixed differences),IS (informative sites), tf (target frequency), ncd2
#' @export
#' @examples ncd2(x=ncd2_input, tf=0.5, w=3000, ncores=2, minIS=2)
#' @import data.table
#' @importFrom data.table ":="
#'

ncd2 <- function(x = x,
                 tf = 0.5,
                 w = 3000,
                 ncores = 2,
                 minIS = 2,
                 by= "POS") {
        assertthat::assert_that(length(unique(x[, CHR])) == 1,
                                msg = "Run one chromosome at a time\n")
        x <- copy(x)
        x[, AF := tx_1 / tn_1]  #allele relative frequency
        x[, AF2 := tx_2 / tn_2]  #allele relative frequency
        x[, ID := seq_along(CHR)]
        w1 <- w / 2
        polpos <- x[AF != 1 & AF != 0]$ID #select positions that are polymorphic
        fdpos <- sort(c(x[AF == 1 & AF2 == 0]$ID, #fixed difference
                        x[AF == 0 & AF2 == 1]$ID)) #also fixed difference
        #to do: check that the there is no intersection between polpos and fdpos.
        x[, SNP := ifelse(ID %in% polpos, T, F)] #logical: True if SNP, False if not.
        x[, fd := ifelse(ID %in% fdpos, T, F)] #logical: True if FD, False if not.
        x[, MAF := ifelse(AF > 0.5, 1 - AF, AF)]
        x[, CHR := stringr::str_replace(stringr::str_replace(CHR,"-",""), ":","")]
        ####################################################################################
        if(by=="POS"){
                #windows (sliding)
                #to do: add some checks here
                vec<-data.table(start=seq(from=x$POS[1], to=x$POS[nrow(x)], by=w1))
                vec[,end:=start+w]
                setkey(vec, start, end)
                x[,start:=POS][,end:=POS]
                res_0<-foverlaps(x, vec, type="within")[, Win.ID:=paste0(CHR,"_",start,"_",end)][,.(POS, ID,SNP,fd,MAF, Win.ID)]
                res_0<-res_0[SNP==T | fd==T][, tf:=tf]
                res_1<-unique(res_0[,.(POS=POS,
                                       S = sum(SNP),
                                       FD = sum(fd),
                                       IS = sum(SNP) + sum(fd),
                                       tf = tf),
                                    by= Win.ID])

                setkey(res_0, Win.ID, POS, tf)
                setkey(res_1, Win.ID, POS, tf)
                res_2<- setorder(res_0[res_1], cols="POS")
                #to do: check that each Win.ID has the number of rows given by IS
                res_3<-res_2[, .(ncd2 = sqrt(sum((MAF-tf)^2)/IS)), by=Win.ID]
                res_3<-unique(res_3, by = "Win.ID")
                res_3b<-res_2[,.(POS,Win.ID,S, FD, IS,tf)]
                res_3b<-unique(res_3b, by="Win.ID")
                setkey(res_3, "Win.ID")
                setkey(res_3b, "Win.ID")
                res4 <- setorder(data.table::merge.data.table(res_3, res_3b)[IS>=minIS],cols="POS")[,POS:=NULL]

        }else if (by=="IS"){
                #to do: write this code.
        }
        return(res4)
}
