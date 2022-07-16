#' Calculate non-central statistic (NCD1)
#'
#' @param x A data.table object
#' @param tf Target frequency
#' @param fold Logical. If TRUE, NCD1 will use minor allele frequencies.
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
#' @export
#'
#' @examples ncd1(x=ncd1_input, selectwin="mid", targetpos=15000, mid=TRUE)
#' @import data.table
#' @importFrom data.table ":="
ncd1 <- function(x = x,
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
        Win.ID <-
                SegSites <-
                POS <- V1 <- temp <- NCD1 <- CHR <- AF <- NULL
        tx_1 <-
                tn_1 <-
                AF2 <-
                tx_2 <- tn_2 <- ID <- SNP <- FD <- MAF <- NULL
        tictoc::tic("Total runtime")
        assertthat::assert_that(length(unique(x[, CHR])) == 1, msg = "Run one
                          chromosome at a time\n")
        if (mid == TRUE) {
                assertthat::assert_that(is.null(targetpos) == FALSE, msg = "NCD_mid requires a targetpos")
                assertthat::assert_that(selectwin == "mid", msg = "If mid=T, selectwin must be=='mid'.")
        }
        if (selectwin == "mid") {
                assertthat::assert_that(mid == TRUE, msg = "If selectwin=='mid', mid must be TRUE.")
                assertthat::assert_that(mid == TRUE, msg = "NCD_mid requires a targetpos.")
        }

        x[, AF := tx_1 / tn_1]
        x[, ID := seq_along(CHR)]
        w1 <- w / 2
        polpos <- x[AF != 1 & AF != 0]$ID

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
        mylist <-
                do.call(rbind,
                        parallel::mclapply(1:length(mylist), function(y)
                                mylist[[y]][, Mid := x$POS[y]][, Win.ID :=
                                                                       y],
                                mc.cores = ncores))
        mylist <- data.table::setDT(mylist)
        # }

        if (mid == TRUE) {
                mylist[, temp := abs(Mid - targetpos), by = Win.ID]
                mylist <- mylist[temp == min(temp)]
                mylist[, temp := NULL]
                mylist[, tf := round(mylist[POS == Mid]$MAF[1], 2)]
                res <-
                       res<- mylist[POS!=Mid][, .(SegSites = sum(SNP),
                                   IS = sum(SNP)),
                               by = Win.ID]
                res[, tf := round(mylist[1, tf], 2)]
                res1 <- dplyr::bind_cols((
                        mylist %>%
                                dplyr::summarise(MidMaf = MAF[which(Mid == POS)],
                                                 Mid = Mid[1], Win.ID=Win.ID[1])
                ),
                (mylist[POS != Mid] %>% dplyr::summarise(CenMaf = max(
                        abs(MAF - tf)
                )))) %>% as.data.table  #target SNP excluded

                res2 <- merge(res, res1)
                res3 <-
                        mylist %>% dplyr::filter(SNP == T) %>% dplyr::filter(POS!=Mid) %>% #exclude target SNP
                        dplyr::summarise(temp2 = sum((MAF - tf) ^ 2)) %>% dplyr::ungroup() %>%
                        as.data.table


        } else{
                mylist[, tf := tf]
                res <-
                        mylist[, .(SegSites = sum(SNP),
                                   IS = sum(SNP)),
                               by = Win.ID]
                res[,tf:=tf]
                res1 <- mylist %>% dplyr::group_by(Win.ID) %>%
                        dplyr::summarise(
                                MidMaf = MAF[which(Mid == POS)],
                                Mid = Mid[1],
                                CenMaf = max(abs(MAF - tf), Win.ID=Win.ID)
                        ) %>%
                        dplyr::ungroup() %>%
                        as.data.table

                res2 <- merge(res, res1)
                res3 <-
                        mylist %>% dplyr::filter(SNP == T) %>% dplyr::group_by(Win.ID) %>%
                        dplyr::summarise(temp2 = sum((MAF - tf) ^ 2)) %>% dplyr::ungroup() %>%
                        as.data.table
        }

        res4 <- merge(res2, res3) %>% as.data.table
        res4 <- res4 %>% dplyr::ungroup() %>% as.data.table
        res4[, NCD1 := sqrt(temp2 / IS), by = Win.ID]
        res4[, temp2 := NULL]

        if (is.null(minIS) == "FALSE") {
                res4 <- res4[IS >= minIS]
        }
        if (selectwin == "val") {
                res4 <-
                        res4[which.min(res4$NCD1),][, .(Win.ID,
                                                        SegSites,
                                                        IS,
                                                        CenMaf,
                                                        Mid,
                                                        MidMaf,
                                                        NCD1,
                                                        tf)]
        } else if (selectwin == "maf") {
                res4 <-
                        res4[which.min(res4$CenMaf),][, .(Win.ID,
                                                          SegSites,
                                                          IS,
                                                          CenMaf,
                                                          Mid,
                                                          MidMaf,
                                                          NCD1,
                                                          tf)]
        } else if (selectwin == "mid") {
                res4 <- res4[which.min(res4$NCD)][, .(Win.ID,
                                                      SegSites,
                                                      IS,
                                                      CenMaf,
                                                      Mid,
                                                      MidMaf,
                                                      NCD1,
                                                      tf)]
        } else{
                res4 <- res4[, .(Win.ID,
                                 SegSites,
                                 IS,
                                 CenMaf,
                                 Mid,
                                 MidMaf,
                                 NCD1,
                                 tf)]
        }

        setnames(res4,
                 "NCD1",
                 ifelse(mid == TRUE, "NCD1_mid", glue::glue("NCD1_{tf}")))

        if (is.null(label) == F) {
                res4 <- res4[, label := label]
        }
        if (verbose == T) {
                tictoc::toc()
        }
        print(res4)
}
