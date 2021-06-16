context("metaGeneProfile")
library(cliProfiler)

test_that("metaGeneProfile works as expected",{
    testpath <- system.file("extdata", package = "cliProfiler")
    test <- readRDS(file.path(testpath, "test.rds"))
    test_gff3 <- file.path(testpath, "annotation_test.gff3")

    output <- metaGeneProfile(
        object = test, annotation = test_gff3,
        include_intron = FALSE
    )

    expect_output(str(output), "List of 2")
    expect_is(output$Peaks, "GRanges")
    expect_is(output$Plot, "ggplot")
})

test_that("metaGeneProfile gives correct error messages",{
    testpath <- system.file("extdata", package = "cliProfiler")
    test <- readRDS(file.path(testpath, "test.rds"))
    test_gff3 <- file.path(testpath, "annotation_test.gff3")
})
