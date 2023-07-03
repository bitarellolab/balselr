#' Calculate non-central statistic (NCD2)
#'
#' @param x A data.table object
#' @param tf Target frequency
#' @param fold Logical. If TRUE, NCD2 will use minor allele frequencies.
#' @param w Window size in bp. Default is 1000
#' @param by.snp Logical. If TRUE, windows are defined around each SNP
#' in the input data. Else, slidding windows in the range first pos:last pos
#' will be used.
#' @param mid Logical. If TRUE runs NCD centered on a core SNP frequency instead of a pre-defined tf. Requires targetpos.
#' @param ncores Number of cores. Increasing this can spead things up for you.
#' Default is 4.
#'@param selectwin Select window. "val" returns the window with the lowest test. "maf" returns the window
#' which has the lowest difference between its central MAF and the tf."mid" returns the window
#' centered closest on a value provided by the user under targetpos."all" returns all windows (idea for when the user wants to select windows after scanning a genome)
#' @param targetpos Position of the target SNP. Only needed for mid==T
#' @param minIS Minimum number of informative sites. Default is 2. Windows with
#'  less informative sites than this threshold are discarded.
#' @param label An optional label to include as the last column of the output
#' @param verbose Logical. If TRUE, progress reports will be printed as the
#' function runs.
#' @return A data.table object
#'
#' @examples ncd2(x=ncd2_input, selectwin="mid", targetpos=15000, mid=TRUE)
#' ncd2(x=ncd2_input, selectwin="val", targetpos=NULL, mid=F)
#' @import data.table
#' @importFrom data.table ":="
ncd2 <- function(
        x = x,
        tf = 0.5,
        fold = T,
        w = 3000,
        by.snp = TRUE,
        mid = FALSE,
        ncores = 2,
        selectwin = "val",
        targetpos = NULL,
        minIS = 2,
        label = NULL,
        verbose = T) {
        Win.ID <- SegSites <- POS <- V1 <- temp <- NCD2 <- CHR <- AF <- AF2 <- NULL
        tx_1 <- tn_1 <- tx_1 <- tx_2 <- tn_2 <- ID <- SNP <- FD <- MAF <- NULL
        tictoc::tic("Total runtime")
        assertthat::assert_that(length(unique(x[, CHR])) == 1,
                                msg = "Run one chromosome at a time\n")
        if (mid == TRUE) {
                assertthat::assert_that(is.null(targetpos) == FALSE, msg = "NCD_mid requires a targetpos")
                assertthat::assert_that(selectwin == "mid", msg = "If mid=T, selectwin must be=='mid'.")
        }
        if (selectwin == "mid") {
                assertthat::assert_that(mid == TRUE, msg = "If selectwin=='mid', mid must be TRUE.")
                assertthat::assert_that(mid == TRUE, msg = "NCD_mid requires a targetpos.")
        }

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
                parallel::mclapply(x2[polpos2, ]$POS, function(y) {
                        x2[POS >= y - w1 &
                                  POS < y + w1][, .(POS, AF, AF2, ID, SNP, FD, tx_1, tx_2, MAF)]
                },
                mc.cores =
                        ncores)
        mylist <-
                do.call(rbind,
                        parallel::mclapply(1:length(mylist), function(y)
                                mylist[[y]][, Mid := x2[polpos2,]$POS[y]][, Win.ID :=
                                                                                y],
                                mc.cores = ncores))
        mylist <- data.table::setDT(mylist)
        # }
        ####################
        if (mid == TRUE) {
                mylist[, temp := abs(Mid - targetpos), by = Win.ID]
                mylist <- mylist[temp == min(temp)]
                mylist[, temp := NULL]
                mylist[, tf := round(mylist[POS == Mid]$MAF[1], 2)]

                res <-
                        mylist[POS != Mid][, .(SegSites = sum(SNP),
                                               FDs = sum(FD),
                                               IS = sum(SNP)+sum(FD)),
                                           by = Win.ID]
                res[, tf := round(mylist[1, tf], 2)]
                res1 <- dplyr::bind_cols((
                        mylist %>%
                                dplyr::summarise(
                                        MidMaf = MAF[which(Mid == POS)],
                                        Mid = Mid[1],
                                        Win.ID = Win.ID[1]
                                )
                ),
                (mylist[POS != Mid] %>% dplyr::summarise(CenMaf = max(
                        abs(MAF - tf)
                )))) %>% as.data.table  #target SNP excluded

                res2 <- merge(res, res1)
                res3 <-
                        mylist %>% dplyr::filter(SNP == T) %>%
                        dplyr::filter(POS != Mid) %>% #exclude target SNP
                        dplyr::summarise(temp2 = sum((MAF - tf) ^ 2), Win.ID = Win.ID[1]) %>%
                        as.data.table
        } else{
                mylist[, tf := tf]
                res <-
                        mylist[, .(SegSites = sum(SNP),
                                   FDs = sum(FD),
                                   IS = sum(SNP)+sum(FD)),
                               by = Win.ID]
                res[, tf := tf]
                res1 <- mylist %>% dplyr::group_by(Win.ID) %>%
                        dplyr::summarise(
                                MidMaf = MAF[which(Mid == POS)],
                                Mid = Mid[1],
                                CenMaf = max(abs(MAF - tf)), Win.ID = Win.ID) %>%
                        dplyr::ungroup() %>%
                        as.data.table

                res2 <- merge(res, res1)
                res3 <-
                        mylist %>% dplyr::filter(SNP == T) %>% dplyr::group_by(Win.ID) %>%
                        dplyr::summarise(temp2 = sum((MAF - tf) ^ 2)) %>% dplyr::ungroup() %>%
                        as.data.table
        }

        ####################

        res4 <- merge(res2, res3) %>% as.data.table
        res4 <- res4 %>% dplyr::ungroup() %>% as.data.table
        mini_fun<-function(x,y){
        sum((rep(0, x)-y)^2)
        }
        res4[, temp4 := mini_fun(unique(FDs),unique(tf)), by=Win.ID]
        res4[, NCD2 := sqrt(((temp2 + temp4) / IS)), by = Win.ID]
        res4[, temp2 := NULL]
        res4[, temp4 := NULL]
        if (is.null(minIS) == "FALSE") {
                res4 <- res4[IS >= minIS]
        }
        if (selectwin == "val") {
                res4 <-
                        res4[which.min(res4$NCD2), ][, .(Win.ID,
                                                         SegSites,
                                                         IS,
                                                         FDs,
                                                         CenMaf,
                                                         Mid,
                                                         MidMaf,
                                                         NCD2,
                                                         tf)]

        } else if (selectwin == "maf") {
                res4 <-
                        res4[which.min(res4$CenMaf), ][, .(Win.ID,
                                                           SegSites,
                                                           IS,
                                                           FDs,
                                                           CenMaf,
                                                           Mid,
                                                           MidMaf,
                                                           NCD2,
                                                           tf)]

        } else if (selectwin == "mid") {
                res4 <- res4[which.min(res4$NCD2)][, .(Win.ID,
                                                      SegSites,
                                                      IS,
                                                      FDs,
                                                      CenMaf,
                                                      Mid,
                                                      MidMaf,
                                                      NCD2,
                                                      tf)]

        } else{
                res4 <- res4[, .(Win.ID,
                                 SegSites,
                                 IS,
                                 FDs,
                                 CenMaf,
                                 Mid,
                                 MidMaf,
                                 NCD2,
                                 tf)]
        }

        setnames(res4,"NCD2",
                 ifelse(mid == TRUE, "NCD2_mid", glue::glue("NCD2_{tf}")))



        if (is.null(label) == F) {
                res4 <- res4[, label := label]
        }
        if (verbose == T) {
                tictoc::toc()
        }
        print(res4)
}
