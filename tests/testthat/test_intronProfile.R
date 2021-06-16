context("intronProfile")
library(cliProfiler)

test_that("intronProfile works as expected",{
    testpath <- system.file("extdata", package = "cliProfiler")
    test <- readRDS(file.path(testpath, "test.rds"))
    test_gff3 <- file.path(testpath, "annotation_test.gff3")

    output <- intronProfile(
        object = test, annotation = test_gff3
    )

    expect_output(str(output), "List of 2")
    expect_is(output$Peaks, "GRanges")
    expect_is(output$Plot, "ggplot")
})
