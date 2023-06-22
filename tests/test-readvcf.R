#the output should be a data.table

testthat::test_that("output is a data table", {
        readinvcf <-  read_vcf(x="inst/example.vcf", only.bi = TRUE)
        testthat::expect_type(readinvcf, "data.frame")
})

testthat::test_that("read_vcf correctly renames column", {
        readinvcf <- read_vcf(vcf_file_path, only.bi = TRUE)
        testthat::expect_true("CHR" %in% colnames(readinvcf))
})
