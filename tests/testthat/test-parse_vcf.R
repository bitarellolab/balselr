test_that("Test error is thrown when file doesn't exist", {
        expect_error(parse_vcf(infile = "filethatdoesnotexist.vcf", n0 = 108, n1 = 1, type = "ncd2"),
                     "VCF file filethatdoesnotexist.vcf does not exist.")
})

test_that("Test error is thrown when type is not specified", {
        expect_error(parse_vcf(infile = "example.vcf", n0 = 108, n1 = 1),
                     "You must choose either 'ncd1' or 'ncd2'")
})

test_that("Test error is thrown when n1 is NULL for ncd2 type", {
        expect_error(parse_vcf(infile = "example.vcf", n0 = 108, type = "ncd2"),
                     "n1 cannot be NULL. NCD2 requires an outgroup.")
})

test_that("Test error is thrown when fold is F for ncd1 type", {
        expect_error(parse_vcf(infile = "example.vcf", n0 = 108, n1 = 1, type = "ncd1", fold = F),
                     "Only the folded option is available when Ancestral/Derived states are unknown.")
})

