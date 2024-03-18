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

