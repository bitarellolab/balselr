testthat::test_that("output is a data table", {
        readinvcf <-  read_vcf(x="example.vcf")
        expect_s3_class(readinvcf, "data.table")
})
