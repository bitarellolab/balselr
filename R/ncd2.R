#' Calculate non-central statistic (NCD2)
#'
#' @param x A data.table object
#' @param tf Target frequency
#' @param fold Logical. If TRUE, NCD2 will use minor allele frequencies.
#' @param w Window size in bp. Default is 1000
#' @param by.snp Logical. If TRUE, windows are defined around each SNP
#' in the input data. Else, slidding windows in the range first pos:last pos
#' will be used.
#' @param ncores Number of cores. Increasing this can spead things up for you.
#' Default is 4.
#' @param valormaf Value or maf. "val" returns the window with the lowest test
#' statistic value. "maf" returns the window with the highest MaxMaf. Only useful if tf=0.5. Default is val
#' @param minIS Minimum number of informative sites. Default is 2. Windows with less informative sites than this threshold are discarded.
#' @param label An optional label to include as the last column of the output
#' @param verbose Logical. If TRUE, progress reports will be printed as the
#' function runs.
#' @return A data.table object
#' @export
#'
#' @examples ncd2(x=hc_)
#' @import data.table
#' @importFrom data.table ":="
ncd2 <- function(x,
                 tf = 0.5,
                 fold = T,
                 w = 3000,
                 by.snp = T,
                 ncores = 4,
                 valormaf = "val",
                 minIS = 2,
                 label = NULL, verbose=T) {
        Win.ID <- Mid <- POS <- temp4 <- FDs <- NCD2 <- FD <- NULL
        temp2 <- IS <- SegSites <- CenMaf <- MidMaf <- SNP <- NULL
        CHR <- AF <- tx_1 <- tn_1 <- AF2 <- tx_2 <- tn_2 <- ID <- MAF <- NULL

  assertthat::assert_that(length(unique(x[, CHR])) == 1,
                          msg = "Run one chromosome at a time\n")
        tictoc::tic("Total runtime")
  x[, AF := tx_1 / tn_1]
  x[, AF2 := tx_2 / tn_2]
  x[, ID := seq_along(CHR)]
  w1 <- 3000 / 2
  polpos <- x[AF != 1 & AF != 0]$ID
  fdpos <- sort(c(x[AF == 1 & AF2 == 0]$ID, x[AF == 0 &
    AF2 == 1]$ID))
  x[, SNP := ifelse(ID %in% polpos, T, F)]
  x[, FD := ifelse(ID %in% fdpos, T, F)]
  # if(by.snp==T){
  mylist <-
    parallel::mclapply(x[polpos, ]$POS, function(y) {
      x[POS >= y - w1 &
        POS < y + w1][, .(POS, AF, AF2, ID, SNP, FD, tx_1, tx_2)]
    },
    mc.cores =
      ncores
    )
  mylist <-
    do.call(
      rbind,
      parallel::mclapply(
              1:length(mylist),function(y)
                      mylist[[y]][,Mid:=x[polpos,]$POS[y]][,Win.ID:=y],
                         mc.cores=ncores))
  mylist <- data.table::setDT(mylist)
  mylist[, MAF := ifelse(AF > 0.5, 1 - AF, AF)]

  # }
  mylist <- data.table::setDT(mylist)
  mylist[,tf:=tf]
  res <-
          mylist[, .(SegSites = sum(SNP), FDs = sum(FD), IS = sum(SNP)+sum(FD)),
                 by = Win.ID]
  res[,tf:=tf]

  res1<-mylist %>% dplyr::group_by(Win.ID) %>%
          dplyr::summarise(MidMaf=MAF[which(Mid == POS)], Mid=Mid[1], CenMaf=min(abs(MAF-tf))) %>%
          dplyr::ungroup() %>%
          as.data.table
  res2<-merge(res,res1)
  res3<-mylist %>% dplyr::filter(SNP==T) %>% dplyr::group_by(Win.ID) %>%
          dplyr::summarise(temp2 = sum((MAF-tf)^2)) %>% dplyr::ungroup() %>%
          as.data.table
  res4<-merge(res2, res3) %>% as.data.table
  res4<-res4 %>% dplyr::ungroup() %>% as.data.table
  res4[, temp4:=sum(rep((tf^2),FDs)), by=Win.ID]
  res4[,NCD2:=sqrt((temp2+temp4)/IS), by=Win.ID]
  res4[,temp2:=NULL]
  res4[,temp4:=NULL]
  if (is.null(minIS) == "FALSE") {
    res4 <- res4[IS >= minIS]
  } else {
    return(res4)
  }

  if (valormaf == "val") {
    res5 <- res4[which.min(res4$NCD2), ][, .(Win.ID, SegSites, FDs, IS, CenMaf, Mid, MidMaf, NCD2)][, tf := tf]
  } else if (valormaf == "maf") {
    res5 <- res4[which.min(res4$CenMaf), ][, .(Win.ID, SegSites, FDs, IS, CenMaf, Mid, MidMaf, NCD2)][, tf := tf]
  } else if (valormaf == "all") {
    res5<-res4
  }
  if (is.null(label) == F) {
    res5 <- res5[, label := label]
  }
  if(verbose==T){tictoc::toc()}
  print(res5)
}
