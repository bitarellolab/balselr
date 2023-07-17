test_that("Test error is thrown when file doesn't exist", {
        expect_error(parse_vcf(infile = "filethatdoesnotexist.vcf", n0 = 108, n1 = 1, type = "ncd2"),
                     "VCF file filethatdoesnotexist.vcf does not exist.")
})


testthat::test_that("parse_vcf handles invalid inputs", {
        expect_error(parse_vcf(infile = "non-existent-file", n0 = 108, type = "ncd1"))
        expect_error(parse_vcf(infile = "testvcf", n0 = 108, type = "invalid_type"))
        expect_error(parse_vcf(infile = "testvcf", n0 = 108, type = "ncd2"))
        expect_error(parse_vcf(infile = "testvcf", n0 = a, type = "ncd1"))
        expect_error(parse_vcf(infile = "testvcf", n0 = a, n1 = NULL, type = "ncd2"))
})



testthat::test_that("parse_vcf returns a data table", {
        parsed_vcf <-  parse_vcf(infile = "example.vcf", n0 = 108, type = "ncd1")
        expect_s3_class(parsed_vcf, "data.table")
})

# testthat::test_that("parse_vcf correctly generates output file name", {
#         result <- parse_vcf(infile=system.file(package="balselr", "example.vcf"), n0 = 108, type = "ncd1")
#         print(result)
#         outfile <- basename(result$outfile)
#         expect_equal(outfile, "example_ncd1.out")
# })


testthat::test_that("parse_vcf calls .vcf_ncd2 correctly for type = 'ncd2'", {
        result <- parse_vcf(infile = system.file(package="balselr", "example.vcf"), n0 = 108, n1=2, type = "ncd2")
        expect_s3_class(result, "data.table")
        expected_columns <- c("CHR", "POS", "REF", "ALT", "tx_1", "tn_1", "tx_2", "tn_2")
        expect_equal(colnames(result), expected_columns)


})

testthat::test_that("parse_vcf calls .vcf_ncd1 correctly for type = 'ncd1'", {
        result <- parse_vcf(infile = system.file(package="balselr", "example.vcf"), n0 = 108, type = "ncd1")
        expect_s3_class(result, "data.table")
        expected_columns <- c("CHR", "POS", "REF", "ALT", "tx_1", "tn_1")
        expect_equal(colnames(result), expected_columns)


})
