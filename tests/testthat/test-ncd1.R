# testthat::test_that("ncd1 works with valid input", {
#         parsed_vcf <-  parse_vcf(infile = system.file(package="balselr", "example.vcf"), n0 = 108, type = "ncd1")
#         expect_silent(ncd1(x = parsed_vcf, tf = 0.5, w = 3000, ncores = 2, minIS = 2))
# })


testthat::test_that("ncd1 fails with invalid input", {
        testthat::expect_error(ncd1(x = "not a data.table", tf = 0.5, w = 3000, ncores = 2, minIS = 2))
})

testthat::test_that("parse_vcf correctly parses a VCF file", {
        result <- parse_vcf(infile = system.file(package="balselr", "example.vcf"), n0 = 108, type = "ncd1")
        testthat::expect_s3_class(result, "data.table")
})



testthat::test_that( "ncd correctly calculates ncd1",{
        x = ncd1_input
        tf = 0.5
        w = 3000
        x[, AF := tx_1 / tn_1] #allele relative frequency
        x[, ID := seq_along(CHR)]
        w1 <- w / 2
        polpos <-
                x[AF != 1 & AF != 0]$ID #select positions that are polymorphic
        x[, SNP := ifelse(ID %in% polpos, T, F)] #logical: True if SNP, False if not.
        x[, MAF := ifelse(AF > 0.5, 1 - AF, AF)]
        res_0<-x[x$POS >= x$POS[1] & x$POS <= x$POS[1] + w][SNP==T][,Win.ID:=paste0(CHR, "_",x$POS[1])]
        res_1<-unique(res_0[,.(S = sum(SNP),
                               IS = sum(SNP),
                               tf = tf,
                               Win.ID)])
        res_2<-res_0[, .(temp2 = sum((MAF - tf) ^ 2)), by=Win.ID]
        res_3 <- data.table::merge.data.table(res_1, res_2)
        expected_ncd1<-res_3[, ncd1:=sqrt((temp2)/IS)][,temp2:=NULL]
        calc_ncd1<-ncd1(ncd1_input) %>% dplyr::arrange(desc(Win.ID)) %>% as.data.table()
        testthat::expect_equal(calc_ncd1$ncd1[1], expected_ncd1$ncd1)
})
