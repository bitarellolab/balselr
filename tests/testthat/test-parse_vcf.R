test_that("Test error is thrown when file doesn't exist", {
        expect_error(parse_vcf(infile = "filethatdoesnotexist.vcf", n0 = 108, n1 = 1, type = "ncd2"),
                     "VCF file filethatdoesnotexist.vcf does not exist.")
})


testthat::test_that("parse_vcf handles invalid inputs", {
        expect_error(parse_vcf(infile = "non-existent-file", n0 = 108, type = "ncd1"))
        expect_error(parse_vcf(infile = "testvcf", n0 = 108, type = "invalid_type"))
        expect_error(parse_vcf(infile = "testvcf", n0 = 108, type = "ncd2"))
        expect_error(parse_vcf(infile = "testvcf", n0 = a, type = "ncd1"))
})



