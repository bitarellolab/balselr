testthat::test_that("output is a data table", {
        readinvcf <-  read_vcf(x=system.file(package="balselr","example.vcf"))
        testthat::expect_s3_class(readinvcf, "data.table")
})

testthat::test_that("read_vcf correctly renames column", {
        result <- read_vcf(x=system.file(package="balselr","example.vcf"))
        testthat:: expect_true("CHR" %in% colnames(result))
})

