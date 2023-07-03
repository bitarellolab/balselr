#' Calculate non-central statistic (NCD1)
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
#' @examples ncd1(x=ncd1_input, tf=0.5, w=3000, ncores=2, minIS=8)
#' @import data.table
#' @importFrom data.table ":="
ncd1 <- function(x = x,
                 tf = 0.5,
                 w = 3000,
                 ncores = 2,
                 minIS = 2) {
 Win.ID <-SegSites <-POS <- V1 <- temp <- NCD1 <- CHR <- AF <- NULL
        tx_1 <-tn_1 <-AF2 <-tx_2 <- tn_2 <- ID <- SNP <- FD <- MAF <- NULL
        tictoc::tic("Total runtime")
        assertthat::assert_that(length(unique(x[, CHR])) == 1, msg = "Run one
                          chromosome at a time\n")

        x[, AF := tx_1 / tn_1]
        x[, ID := seq_along(CHR)]
        w1 <- w / 2
        polpos <- x[AF != 1 & AF != 0]$ID #must be polymorphic
        x[, SNP := ifelse(ID %in% polpos, T, F)]
        x <- x[SNP == T]
        x[, MAF := ifelse(AF > 0.5, 1 - AF, AF)]
        mylist <-
                parallel::mclapply(x$POS, function(y) {
                        x[POS >= y - w1 &
                                  POS < y + w1][, .(POS, AF, ID, SNP, tx_1, MAF)]
                },
                mc.cores =
                        ncores)
        mylist2 <-
                do.call(rbind,
                        parallel::mclapply(1:length(mylist), function(y)
                                mylist[[y]][, Mid := x$POS[y]][, Win.ID :=
                                                                       y],
                                mc.cores = ncores))
        mylist2 <- data.table::setDT(mylist2)
        # }


                mylist2[, tf := tf]
                res <-
                        mylist2[, .(SegSites = sum(SNP),
                                   IS = sum(SNP)),
                               by = Win.ID]
                res[, tf := tf]
                res1 <- mylist2 %>% dplyr::group_by(Win.ID) %>%
                        dplyr::reframe(
                                MidMaf = MAF[which(Mid == POS)],
                                Mid.SNP = Mid[1]) %>%
                               # TopMaf = which(min(abs(MAF - tf)), Win.ID = Win.ID) %>%
                        #dplyr::ungroup() %>%
                        as.data.table %>%
                        unique() #one row per window

                res2 <- merge(res, res1)
                res3 <-
                        mylist2 %>%
                        dplyr::filter(SNP == T) %>%
                       #dplyr::filter(IS>=minIS) %>%
                        dplyr::group_by(Win.ID) %>%
                        dplyr::reframe(temp2 = sum((MAF - tf) ^ 2)) %>%

                        as.data.table
       # }

        res4 <- merge(res2, res3) %>%
                dplyr::filter(IS>=minIS) %>%
                as.data.table
        res4 <- res4 %>% dplyr::ungroup() %>% as.data.table
        res4[, NCD1 := sqrt(temp2 / IS), by = Win.ID]
        res4[, temp2 := NULL]
        print(res4)
}
