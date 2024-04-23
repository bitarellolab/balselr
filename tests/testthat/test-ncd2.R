# testthat::test_that("ncd1 works with valid input", {
#         parsed_vcf <-  parse_vcf(infile = system.file(package="balselr", "example.vcf"), n0 = 108, type = "ncd1")
#         expect_silent(ncd1(x = parsed_vcf, tf = 0.5, w = 3000, ncores = 2, minIS = 2))
# })


testthat::test_that("ncd2 fails with invalid input", {
        testthat::expect_error(ncd2(x = "not a data.table", tf = 0.5, w = 3000, ncores = 2, minIS = 2))
})

testthat::test_that("parse_vcf correctly parses a VCF file", {
        result <- parse_vcf(infile = system.file(package="balselr", "example.vcf"), n0 = 108, n=2, type = "ncd2")
        testthat::expect_s3_class(result, "data.table")
})



testthat::test_that( "ncd correctly calculates ncd2",{
        x = ncd2_input
        tf = 0.5
        w = 3000
        x[, AF := tx_1 / tn_1]  #allele relative frequency
        x[, AF2 := tx_2 / tn_2]  #allele relative frequency
        x[, ID := seq_along(CHR)]
        w1 <- w / 2
        polpos <- x[AF != 1 & AF != 0]$ID #select positions that are polymorphic
        fdpos <- sort(c(x[AF == 1 & AF2 == 0]$ID, #fixed difference
                        x[AF == 0 & AF2 == 1]$ID)) #also fixed difference
        #to do: check that the there is no intersection between polpos and fdpos.
        x[, SNP := ifelse(ID %in% polpos, T, F)] #logical: True if SNP, False if not.
        x[, fd := ifelse(ID %in% fdpos, T, F)] #logical: True if FD, False if not.
        x[, MAF := ifelse(AF > 0.5, 1 - AF, AF)]
        x[, CHR := stringr::str_replace(stringr::str_replace(CHR,"-",""), ":","")]
        res_0<-x[x$POS >= x$POS[1] & x$POS <= x$POS[1] + w][SNP==T | fd==T][,tf:=0.5][,Win.ID:=paste0(CHR, "_",x$POS[1])]
        res_1<-unique(res_0[,.(POS=POS,
                               S = sum(SNP),
                               FD = sum(fd),
                               IS = sum(SNP) + sum(fd),
                               tf = tf,
                               Win.ID)])
        setkey(res_0, POS, tf, Win.ID)
        setkey(res_1, POS, tf, Win.ID)
        res_2<- setorder(res_0[res_1], cols="POS")
        #to do: check that each Win.ID has the number of rows given by IS
        res_3<-res_2[, .(ncd2 = sqrt(sum((MAF-tf)^2)/IS)), by=Win.ID]
        res_3<-unique(res_3, by = "Win.ID")
        res_3b<-res_2[,.(POS,Win.ID,S, FD, IS, tf)]
        res_3b<-unique(res_3b, by = "Win.ID")
        res_3c <- data.table::merge.data.table(res_3, res_3b)
        expected_ncd2<-res_3c
        calc_ncd2<-ncd2(ncd2_input) %>% dplyr::arrange(desc(Win.ID)) %>% as.data.table()
        testthat::expect_equal(calc_ncd2$ncd2[1], expected_ncd2$ncd2)
})

