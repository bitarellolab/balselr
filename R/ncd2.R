#' Calculate non-central statistic (NCD2)
#'
#' @param x A data.table object
#' @param tf Target frequency

#' @param w Window size in bp. Default is 1000

#' @param ncores Number of cores. Increasing this can spead things up for you.
#' Default is 4.


#' @param minIS Minimum number of informative sites. Default is 2. Windows with
#'  less informative sites than this threshold are discarded.

#' @return A data.table object
#' @export
#'
#' @examples ncd2(x=ncd2_input, tf=0.5, w=3000, ncores=2, minIS=8)
#' @import data.table
#' @importFrom data.table ":="
#'
ncd2 <- function(x = x,
                 tf = 0.5,
                 w = 3000,
                 ncores = 2,
                 minIS = 2) {
        #Win.ID <- SegSites <- POS <- V1 <- temp <- NCD2 <- CHR <- AF <- AF2 <- NULL
        #tx_1 <- tn_1 <- tx_1 <- tx_2 <- tn_2 <- ID <- SNP <- FD <- MAF <- NULL

        assertthat::assert_that(length(unique(x[, CHR])) == 1,
                                msg = "Run one chromosome at a time\n")

        x[, AF := tx_1 / tn_1]
        x[, AF2 := tx_2 / tn_2]
        x[, ID := seq_along(CHR)]
        w1 <- w / 2
        polpos <- x[AF != 1 & AF != 0]$ID
        fdpos <- sort(c(x[AF == 1 & AF2 == 0]$ID, x[AF == 0 & AF2 == 1]$ID))
        x[, SNP := ifelse(ID %in% polpos, T, F)]
        x[, FD := ifelse(ID %in% fdpos, T, F)]
        x[, MAF := ifelse(AF > 0.5, 1 - AF, AF)]
        x2 <- x[SNP == T | FD == T]
        x2[, ID := seq_along(CHR)]
        polpos2<-x2[SNP==T]$ID
        fdpos2<-x2[FD==T]$ID
        # if(by.snp==T){
        mylist <-
                parallel::mclapply(x$POS, function(y) {
                        x[POS >= y - w1 &
                                  POS < y + w1][, .(POS, AF, AF2, FD, ID, SNP, tx_1, MAF)]
                },
                mc.cores =
                        ncores)
        mylist2 <-
                do.call(rbind,
                        parallel::mclapply(1:length(mylist), function(y)
                                mylist[[y]][, Mid := x2[polpos2,]$POS[y]][, Win.ID :=
                                                                                  y],
                                mc.cores = ncores))
        mylist2 <- data.table::setDT(mylist2)

        ####################
                mylist2[, tf := tf]
                res <-
                        mylist2[, .(SegSites = sum(SNP),
                                   FDs = sum(FD),
                                   IS = sum(SNP)+sum(FD)),
                               by = Win.ID]
                res[, tf := tf]
                res1 <- mylist2 %>%
                        dplyr::group_by(Win.ID) %>%
                        dplyr::reframe(
                                MidMaf = MAF[which(Mid == POS)],
                                MidSNP = Mid[1]) %>%
                                #CenMaf = max(abs(MAF - tf)), Win.ID = Win.ID) %>%
                       # dplyr::ungroup() %>%
                        as.data.table %>%
                        unique()

                res2 <- merge(res, res1)
                res3 <-
                        mylist2 %>%
                        dplyr::filter(SNP == T) %>%
                        dplyr::group_by(Win.ID) %>%
                        dplyr::reframe(temp2 = sum((MAF - tf) ^ 2)) %>%
                        as.data.table
        #}

        ####################

        res4 <- merge(res2, res3) %>% as.data.table
        res4 <- res4 %>%
                dplyr::ungroup() %>%
                as.data.table
        mini_fun<-function(x,y){
        sum((rep(0, x)-y)^2)
        }
        res4[, temp4 := mini_fun(unique(FDs),unique(tf)), by=Win.ID]
        res4[, NCD2 := sqrt(((temp2 + temp4) / IS)), by = Win.ID]
        res4[, temp2 := NULL]
        res4[, temp4 := NULL]
        #if (is.null(minIS) == "FALSE") {
        res4 <- res4[IS >= minIS]
        return(print(res4))
}
