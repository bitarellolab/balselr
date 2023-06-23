testthat::test_that("output is a data table", {
        readinvcf <-  read_vcf(x="example.vcf", only.bi = TRUE)
        expect_s3_class(readinvcf, "data.table")
})

testthat::test_that("read_vcf correctly renames column", {
        result <- read_vcf(x="example.vcf", only.bi = TRUE)
        expect_true("CHR" %in% colnames(result))
})

test_that("read_vcf correctly filters data with only.bi parameter", {
        result_bi <- read_vcf(x="example.vcf", only.bi = TRUE)
        result_all <- read_vcf(x="example.vcf", only.bi = FALSE)
        expect_true(nrow(result_bi) <= nrow(result_all))
})

