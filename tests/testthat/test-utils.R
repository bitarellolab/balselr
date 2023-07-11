testthat::test_that("outfile_path correctly modifies file names", {
        fixed_time <- as.POSIXct("2023-07-11 12:00:00", tz = "UTC")

        testthat::with_mock(
                `base::Sys.time` = function() fixed_time,
                {
                        expect_equal(outfile_path(), "example_2023-July-11_12-00-00.out")
                        expect_equal(outfile_path("test.txt"), "test_2023-July-11_12-00-00.out")
                        expect_equal(outfile_path("another_test.jpeg"), "another_test_2023-July-11_12-00-00.out")
                }
        )
})
