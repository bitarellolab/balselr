# testthat::test_that("ncd1 works with valid input", {
#         parsed_vcf <-  parse_vcf(infile = system.file(package="balselr", "example.vcf"), n0 = 108, type = "ncd1")
#         expect_silent(ncd1(x = parsed_vcf, tf = 0.5, w = 3000, ncores = 2, minIS = 2))
# })


test_that("ncd1 fails with invalid input", {
        expect_error(ncd1(x = "not a data.table", tf = 0.5, w = 3000, ncores = 2, minIS = 2))
})

test_that("parse_vcf correctly parses a VCF file", {
        result <- parse_vcf(infile = system.file(package="balselr", "example.vcf"), n0 = 108, type = "ncd1")
        expect_is(result, "data.table")
})



testthat::test_that( "ncd correctly calculates ncd1"{

})

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
x[SNP==T]




