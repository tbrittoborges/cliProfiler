context("spliceSiteProfile")
library(cliProfiler)

test_that("spliceSiteProfile works as expected",{
    testpath <- system.file("extdata", package = "cliProfiler")
    test <- readRDS(file.path(testpath, "test.rds"))
    test_gff3 <- file.path(testpath, "annotation_test.gff3")

    output <- spliceSiteProfile(
        test, test_gff3, flanking=200, bin=40
    )

    expect_output(str(output), "List of 2")
    expect_is(output$Peaks, "GRanges")
    expect_is(output$Plot, "ggplot")
})
