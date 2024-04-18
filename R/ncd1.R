#' Calculate non-central statistic (NCD1)
#'
#' @param x A data.table object with
#' columns: CHR, POS, REF, ALT, tx_1 (number of alternate allele copies),
#' tn_1 (total number of alleles).
#' @param tf Target frequency.
#' @param w Window size in bp. Default is 3000.
#' @param ncores Number of cores. Increasing this can speed things up for you.
#' Default is 2.
#' @param minIS Minimum number of informative sites. Default is 2. Windows with
#'  less informative sites than this threshold are discarded.
#' @param by Define how to scan the genome. "POS" (default) defined sliding windows based on w. "IS" defined windows around each informative site.
#' @return A data.table object with columns:  Win.ID, S (sites), IS (informative sites), tf, ncd1
#' @export
#'
#' @examples ncd1(x=ncd1_input, tf=0.5, w=3000, ncores=2, minIS=8)
#' @import data.table
#' @importFrom data.table ":="
#'
# to do: check that ncd1_input (internal object) is the same as the output of example command for ncd1 in the parse_vcf function.
# something like ncd1(x=parse_vcf(infile=system.file(package="balselr", "example.vcf"), n0=108, type="ncd1"), tf=0.5, w=3000, ncores=2, minIS=8) == ncd1(x=ncd1_input, tf=0.5, w=3000, ncores=2, minIS=8)
# do the equivalent thing for ncd2
ncd1 <- function(x = x,
                 tf = 0.5,
                 w = 3000,
                 ncores = 2,
                 minIS = 2,
                 by="POS") {
        assertthat::assert_that(length(unique(x[, CHR])) == 1, msg = "Run one
                          chromosome at a time\n")

        x[, AF := tx_1 / tn_1] #allele relative frequency
        x[, ID := seq_along(CHR)]
        w1 <- w / 2
        polpos <-
                x[AF != 1 & AF != 0]$ID #select positions that are polymorphic
        x[, SNP := ifelse(ID %in% polpos, T, F)] #logical: True if SNP, False if not.
        x[, MAF := ifelse(AF > 0.5, 1 - AF, AF)]
        ####################################################################################
        if(by=="POS"){
        #windows (sliding)
                #to do: add some checks here
                vec<-data.table(start=seq(from=x$POS[1], to=x$POS[nrow(x)], by=w1))
                vec[,end:=start+w]
                setkey(vec, start, end)
                x[,POS2:=POS]
                setnames(x, c("POS","POS2"),c("start","end"))
                res_0<-foverlaps(x, vec, type="within")[, Win.ID:=paste0(CHR,"_",start,"_",end)][, POS:=start][,.(POS,AF, ID,SNP,MAF, Win.ID)]
                res_0<-res_0[SNP==T][, tf:=tf]
                res_1<-unique(res_0[,.(S = sum(SNP),
                        IS = sum(SNP),
                        tf = tf),
                by= Win.ID])
                res_2<-res_0[, .(temp2 = sum((MAF - tf) ^ 2)), by=Win.ID]
                res4 <- merge(res_1, res_2) %>%
                        dplyr::filter(IS >= minIS) %>%
                        as.data.table
                res4<-res4[, ncd1:=sqrt((temp2)/IS)][,temp2:=NULL] %>% dplyr::arrange(Win.ID) %>% as.data.table

        } else if(by=="IS"){
                ###################################
                #windows centered on SNPs
                x2 <- x[SNP == T] #select only polymorphic positions
                mylist <-
                        parallel::mclapply(x2$POS, function(y) {
                                x2[,start:=y-w1][,end:=y+w1][POS >= start &
                                          POS < end][, Mid:=y][,Win.ID:=paste0(CHR,"_",start,"_",end)][, tf:=tf][, .(POS, AF, ID, SNP, MAF, Mid, Win.ID, tf)]
                        },
                        mc.cores =
                                ncores) #creates list where each element is a genomic window

                mylist2 <-
                        do.call(rbind, mylist)

                #remove windows where Win.ID includes negative numbers
                mylist2<-mylist2[-grep("-",mylist2$Win.ID),]
                #remove window where mid==last POS in dataset
                mylist2<-mylist2[Mid!=mylist2$POS[nrow(mylist2)]]
                #to do: check that mylist2 is a data table
        res <-
                unique(mylist2[, .(S = sum(SNP),
                            IS = sum(SNP),
                            tf = tf),
                        by = Win.ID])
        #to do: add some checks here
        res1 <- unique(as.data.table(mylist2[, .(MidMaf = MAF[which(Mid == POS)],
                               Mid = Mid), by=Win.ID]))
                # TopMaf = which(min(abs(MAF - tf)), Win.ID = Win.ID) %>%

        #to do: add some checks here
        res2 <- merge(res, res1)
        #to do: add some checks here
        res3 <-
                mylist2[, .(temp2 = sum((MAF - tf) ^ 2)), by=Win.ID]


        res4 <- merge(res2, res3) %>%
                dplyr::mutate(ncd1:=sqrt((temp2)/IS)) %>%
                dplyr::filter(IS >= minIS) %>%
                dplyr::select(Win.ID, S, IS, tf, MidMaf, Mid) %>%
                dplyr::arrange(ncd1) %>%
                as.data.table

        }
        split_ids <- lapply(strsplit(res4$Win.ID, "_"), as.numeric)
        sort_values <- sapply(split_ids, function(x) x[2])
        res5 <- res4[order(sort_values), ]

        return(res5)
}

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
#' @param label An optional label to include as the last column of the output. Defaults to value in selectwin.
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
                        mylist[POS != Mid][, .(SegSites = sum(SNP),
                                               IS = sum(SNP)),
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
                                   IS = sum(SNP)),
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

        res4 <- merge(res2, res3) %>% as.data.table
        res4 <- res4 %>% dplyr::ungroup() %>% as.data.table
        res4[, NCD1 := sqrt(temp2 / IS), by = Win.ID]
        res4[, temp2 := NULL]

        if (is.null(minIS) == "FALSE") {
                res4 <- res4[IS >= minIS]
        }
        if (selectwin == "val") {
                res4 <-
                        res4[which.min(res4$NCD1), ][, .(Win.ID,
                                                         SegSites,
                                                         IS,
                                                         CenMaf,
                                                         Mid,
                                                         MidMaf,
                                                         NCD1,
                                                         tf)]
        } else if (selectwin == "maf") {
                res4 <-
                        res4[which.min(res4$CenMaf), ][, .(Win.ID,
                                                           SegSites,
                                                           IS,
                                                           CenMaf,
                                                           Mid,
                                                           MidMaf,
                                                           NCD1,
                                                           tf)]
        } else if (selectwin == "mid") {
                res4 <- res4[which.min(res4$NCD1)][, .(Win.ID,
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
                res4 <- res4[, Method := label]
        } else{
                res4<-res4[, Method := paste(selectwin)]
        }
        if (verbose == T) {
                tictoc::toc()
        }
        print(res4)
}
