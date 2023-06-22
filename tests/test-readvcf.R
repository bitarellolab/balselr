#the output should be a data.table

testthat::test_that("output is a data table", {
        readinvcf <-  read_vcf(x="inst/example.vcf")
        testthat::expect_type(readinvcf, "data.frame")
})


