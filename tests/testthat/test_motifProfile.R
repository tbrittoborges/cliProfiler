context("motifProfile")
library(cliProfiler)

test_that("motifProfile works as expected",{
    testpath <- system.file("extdata", package = "cliProfiler")
    test <- readRDS(file.path(testpath, "test.rds"))
    test_gff3 <- file.path(testpath, "annotation_test.gff3")

    output <- motifProfile(
        object = test, motif = "DRACH",
        genome = "BSgenome.Mmusculus.UCSC.mm10", flanking = 20
    )

    expect_output(str(output), "List of 2")
    expect_is(output$Numbers, "data.frame")
    expect_is(output$Plot, "ggplot")
})
