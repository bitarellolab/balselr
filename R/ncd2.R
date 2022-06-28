#' Calculate non-central statistic (NCD2)
#'
#' @param x
#' @param tf
#' @param fold
#' @param w
#' @param by.snp
#'
#' @return
#' @export
#'
#' @examples
#' @import data.table
#' @importFrom data.table ":="
ncd2 <- function(x = x,
                 tf = 0.5,
                 fold = T,
                 w = 3000,
                 by.snp = T,
                 ncores=4,
                 valormaf="val",
                 minIS=10,
                 label="example_sim") {
        assertthat::assert_that(length(unique(x[, CHR])) == 1, msg = "Run one chromosome at a time\n")
        x[, MAF := tx_1 / tn_1]
        x[, MAF2 := tx_2 / tn_2]
        x[, ID := seq_along(CHR)]
        w1 = 3000 / 2
        polpos <- x[MAF != 1 & MAF != 0]$ID
        fdpos <- sort(c(x[MAF == 1 & MAF2 == 0]$ID, x[MAF == 0 &
                                                              MAF2 == 1]$ID))
        x[, SNP := ifelse(ID %in% polpos, T, F)]
        x[, FD := ifelse(ID %in% fdpos, T, F)]
        #if(by.snp==T){
        mylist <-
                parallel::mclapply(x[polpos, ]$POS, function(y)
                        x[POS >= y - w1 &
                                  POS < y + w1][, .(POS, MAF, MAF2, ID, SNP, FD, tx_1, tx_2)][, Mid := y], mc.cores =
                                ncores)
        mylist <-
                do.call(rbind,
                        parallel::mclapply(seq_along(mylist), function(y)
                                mylist[[y]][, Win.ID := y], mc.cores = ncores))
        mylist <- data.table::setDT(mylist)
        mylist[, MAF := ifelse(MAF > 0.5, 1 - MAF, MAF)]
        #}
        res <-
                mylist[, .(SegSites = sum(SNP), FDs = sum(FD)), by = Win.ID][, IS := SegSites +
                                                                                     FDs]
        res1 <- mylist[, (MidMaf = MAF[which(Mid == POS)]), by = Win.ID]
        res2 <- mylist[, (maf = MAF[sort(c(which(SNP), which(FD)))]), by = Win.ID]
        res3 <- res2[, (MaxMaf = max(V1)), by = Win.ID]
        res4 <- mylist[, (Mid = Mid[1]), by = Win.ID]

        data.table::setkey(res, Win.ID)
        data.table::setkey(res1, Win.ID)
        data.table::setkey(res3, Win.ID)
        data.table::setkey(res4, Win.ID)
        res5 <- res[res3, ][res4, ][res1, ]
        data.table::setnames(res5,
                             c(
                                     "Win.ID",
                                     "SegSites",
                                     "FDs",
                                     "IS",
                                     "MaxMaf",
                                     "Mid",
                                     "MidMaf"
                             ))
        res5[, temp := res2[, sum((V1 - tf) ^ 2), by = Win.ID]$V1]
        res5[, NCD2 := sqrt(temp / IS)]

        res5[, temp:=NULL]
        if(is.null(minIS)=="FALSE"){
         res5<-res5[IS>=minIS]

        }else{return(res5)}

        if(valormaf=="val"){
                res6<-res5[which.min(res5$NCD2), ][, .(Win.ID, SegSites, FDs, IS, MaxMaf, Mid, MidMaf, NCD2)][, tf := tf]
        }else if(valormaf=="maf"){
                res6<-res5[which.max(res5$MaxMaf), ][, .(Win.ID, SegSites, FDs, IS, MaxMaf, Mid, MidMaf, NCD2)][, tf := tf]
        }else if(valormaf=="all"){
                res6<-res5[which.min(res5$NCD1), ][, .(Win.ID, SegSites,FDs,IS, MaxMaf, Mid, MidMaf, NCD2)][, tf := tf]
                res6<-rbind(res6, res5[which.max(res5$MaxMaf), ][, .(Win.ID, SegSites, FDs,IS, MaxMaf, Mid, MidMaf, NCD2)][, tf := tf])
        }
        if(is.null(label)==F){
                res6<-res5[, label:=label]
        }
        print(res6)
}
