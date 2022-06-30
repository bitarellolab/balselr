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
#' @return A data.table object
#' @export
#'
#' @examples parse_vcf(infile = "inst/example.vcf", nind = c(108, 1), fold=T)
#' @import data.table
#' @importFrom data.table ":="
ncd2 <- function(x = x,
                 tf = 0.5,
                 fold = T,
                 w = 3000,
                 by.snp = T,
                 ncores = 4,
                 valormaf = "val",
                 minIS = 2,
                 label = NULL, verbose=T) {
  assertthat::assert_that(length(unique(x[, CHR])) == 1, msg = "Run one chromosome at a time\n")
        tictoc::tic("Total runtime")
  x[, AF := tx_1 / tn_1]
  x[, AF2 := tx_2 / tn_2]
  x[, ID := seq_along(CHR)]
  w1 <- 3000 / 2
  polpos <- x[MAF != 1 & MAF != 0]$ID
  fdpos <- sort(c(x[AF == 1 & AF2 == 0]$ID, x[AF == 0 &
    AF2 == 1]$ID))
  x[, SNP := ifelse(ID %in% polpos, T, F)]
  x[, FD := ifelse(ID %in% fdpos, T, F)]
  # if(by.snp==T){
  mylist <-
    parallel::mclapply(x[polpos, ]$POS, function(y) {
      x[POS >= y - w1 &
        POS < y + w1][, .(POS, AF, AF2, ID, SNP, FD, tx_1, tx_2)][, Mid := y]
    },
    mc.cores =
      ncores
    )
  mylist <-
    do.call(
      rbind,
      parallel::mclapply(seq_along(mylist), function(y) {
        mylist[[y]][, Win.ID := y]
      }, mc.cores = ncores)
    )
  mylist <- data.table::setDT(mylist)
  mylist[, MAF := ifelse(AF > 0.5, 1 - AF, AF)]
  # }
  res <-
    mylist[, .(SegSites = sum(SNP), FDs = sum(FD)), by = Win.ID][, IS := SegSites +
      FDs]
  mylist[,tf:=tf]
  res1 <- mylist[, (MidMaf <- MAF[which(Mid == POS)]), by = Win.ID]
  res2 <- mylist[, (maf <- MAF[c(which(SNP), which(FD))]), by = Win.ID]
  res3 <- res2[, (CenMaf <- max(abs(V1-tf))), by = Win.ID]
  res4 <- mylist[, (Mid <- Mid[1]), by = Win.ID]

  data.table::setkey(res, Win.ID)
  data.table::setkey(res1, Win.ID)
  data.table::setkey(res3, Win.ID)
  data.table::setkey(res4, Win.ID)
  res5 <- res[res3, ][res4, ][res1, ]
  data.table::setnames(
    res5,
    c(
      "Win.ID",
      "SegSites",
      "FDs",
      "IS",
      "CenMaf",
      "Mid",
      "MidMaf"
    )
  )
  res5[, temp := res2[, ((V1 - tf)^2), by = Win.ID]$V1]
  res5[, NCD2 := sqrt(temp / IS)]

  res5[, temp := NULL]
  if (is.null(minIS) == "FALSE") {
    res5 <- res5[IS >= minIS]
  } else {
    return(res5)
  }

  if (valormaf == "val") {
    res6 <- res5[which.min(res5$NCD2), ][, .(Win.ID, SegSites, FDs, IS, CenMaf, Mid, MidMaf, NCD2)][, tf := tf]
  } else if (valormaf == "maf") {
    res6 <- res5[which.min(res5$CenMaf), ][, .(Win.ID, SegSites, FDs, IS, CenMaf, Mid, MidMaf, NCD2)][, tf := tf]
  } else if (valormaf == "all") {
    res6 <- res5[which.min(res5$NCD1), ][, .(Win.ID, SegSites, FDs, IS,CenMaf, Mid, MidMaf, NCD2)][, tf := tf]
    res6 <- rbind(res6, res5[which.max(res5$MaxMaf), ][, .(Win.ID, SegSites, FDs, IS, CenMaf, Mid, MidMaf, NCD2)][, tf := tf])
  }
  if (is.null(label) == F) {
    res6 <- res5[, label := label]
  }
  if(verbose==T){tictoc::toc()}
  print(res6)
}
