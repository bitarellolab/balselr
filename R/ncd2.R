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
#' @return A data.table object with columns: TO DO
#' @export
#'
#' @examples ncd2(x=ncd2_input, tf=0.5, w=3000, ncores=2, minIS=8)
#' @import data.table
#' @importFrom data.table ":="
#'
# to do: check that ncd2_input (internal object) is the same as the output of example command for ncd2 in the parse_vcf function.
ncd2 <- function(x = x,
                 tf = 0.5,
                 w = 3000,
                 ncores = 2,
                 minIS = 2,
                 by= "POS") {


        assertthat::assert_that(length(unique(x[, CHR])) == 1,
                                msg = "Run one chromosome at a time\n")

        x[, AF := tx_1 / tn_1]  #allele relative frequency
        x[, AF2 := tx_2 / tn_2]  #allele relative frequency
        x[, ID := seq_along(CHR)]
        w1 <- w / 2
        polpos <- x[AF != 1 & AF != 0]$ID #select positions that are polymorphic
        fdpos <- sort(c(x[AF == 1 & AF2 == 0]$ID, #fixed difference
                        x[AF == 0 & AF2 == 1]$ID)) #also fixed difference
        #to do: check that the there is no intersection between polpos and fdpos.
        x[, SNP := ifelse(ID %in% polpos, T, F)] #logical: True if SNP, False if not.
        x[, FD := ifelse(ID %in% fdpos, T, F)] #logical: True if FD, False if not.
        x[, MAF := ifelse(AF > 0.5, 1 - AF, AF)]
        x2 <- x[SNP == T | FD == T] #select informative sites
        x2[, ID := seq_along(CHR)]
        polpos2<-x2[SNP==T]$ID
        fdpos2<-x2[FD==T]$ID
        # if(by.snp==T){
        mylist <-
                parallel::mclapply(x2$POS, function(y) {
                        x2[POS >= y - w1 &
                                  POS < y + w1][, Mid:=y][, .(POS, AF, ID, SNP, FD, MAF, Mid)]
                },
                mc.cores =
                        ncores) #creates list where each element is a genomic window

        #remove first and last windows in a chromosome
        mylist<-mylist[seq(from=2, to=length(mylist)-1)]

        mylist2 <-
                do.call(rbind,
                        parallel::mclapply(1:length(mylist), function(y)
                                mylist[[y]][, Win.ID := y][,tf:=tf],
                                mc.cores = ncores))
        #to do: check that mylist2 is a data table
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
