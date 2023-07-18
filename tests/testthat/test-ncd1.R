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



test_that("ncd1 calculates non-central deviation correctly", {
        parsed_vcf <- parse_vcf(infile = system.file(package="balselr", "example.vcf"), n0 = 108, type = "ncd1")
        result <- ncd1(x = parsed_vcf, tf = 0.5, w = 3000, ncores = 2, minIS = 2)
        window1 <- parsed_vcf[POS >= parsed_vcf$POS[1] - 1500 & POS < parsed_vcf$POS[1] + 1500]
        window1 <- window1[AF != 1 & AF != 0]
        expected_ncd1 <- sqrt(sum((window1$MAF - 0.5) ^ 2) / nrow(window1))

        expect_equal(result$NCD1[1], expected_ncd1)
})
