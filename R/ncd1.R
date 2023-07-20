#' Calculate non-central statistic (NCD1)
#'
#' @param x A data.table object with
#' columns: CHR, POS, REF, ALT, tx_1 (number of alternate allele copies), tn_1 (total number of alleles).
#' @param tf Target frequency.
#' @param w Window size in bp. Default is 1000.
#' @param ncores Number of cores. Increasing this can speed things up for you.
#' Default is 2.
#' @param minIS Minimum number of informative sites. Default is 2. Windows with
#'  less informative sites than this threshold are discarded.
#' @return A data.table object with columns:  Win.ID SegSites IS tf MidMaf Mid ncd1
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
                 minIS = 2) {
        assertthat::assert_that(length(unique(x[, CHR])) == 1, msg = "Run one
                          chromosome at a time\n")

        x[, AF := tx_1 / tn_1] #allele relative frequency
        x[, ID := seq_along(CHR)]
        w1 <- w / 2
        polpos <-
                x[AF != 1 & AF != 0]$ID #select positions that are polymorphic
        x[, SNP := ifelse(ID %in% polpos, T, F)] #logical: True if SNP, False if not.
        x <- x[SNP == T] #select only polymorphic positions
        x[, MAF := ifelse(AF > 0.5, 1 - AF, AF)]
        mylist <-
                parallel::mclapply(x$POS, function(y) {
                        x[POS >= y - w1 &
                                  POS < y + w1][, Mid:=y][, .(POS, AF, ID, SNP, MAF, Mid)]
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
        res <-
                mylist2[, .(SegSites = sum(SNP),
                            IS = sum(SNP),
                            tf = tf),
                        by = Win.ID] %>% unique()
        #to do: add some checks here
        res1 <- mylist2 %>% dplyr::group_by(Win.ID) %>%
                dplyr::reframe(MidMaf = MAF[which(Mid == POS)],
                               Mid = Mid) %>%
                # TopMaf = which(min(abs(MAF - tf)), Win.ID = Win.ID) %>%
                #dplyr::ungroup() %>%
                as.data.table %>%
                unique() #one row per window
        #to do: add some checks here
        res2 <- merge(res, res1)
        #to do: add some checks here
        res3 <-
                mylist2 %>%
                dplyr::group_by(Win.ID) %>%
                dplyr::reframe(temp2 = sum((MAF - tf) ^ 2)) %>%
                as.data.table
        #}

        res4 <- merge(res2, res3) %>%
                dplyr::filter(IS >= minIS) %>%
                as.data.table

        res4<-res4[, ncd1 := sqrt(temp2 / IS), by = Win.ID][,.(Win.ID, SegSites, IS, tf, MidMaf, Mid, ncd1)]
        return(res4)
}
