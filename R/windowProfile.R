## ============================================================================
## The windowProfile function for GRanges objects.
## ----------------------------------------------------------------------------

#' @import dplyr
#' @import methods
#' @import ggplot2
#' @importFrom GenomicRanges GRanges
#' @importFrom Biostrings getSeq
#' @importFrom BSgenome getBSgenome
#' @importFrom GenomicRanges start

## ============================================================================
## Small functions
## ----------------------------------------------------------------------------

.windowPosition <- function(object)
{
    ##-----calculate the position-----##
    object_p <- object[strand(object) == "+"]
    object_p$window_map <- (start(object_p) -
        object_p$window_S)/object_p$window_length
    object_n <- object[strand(object) == "-"]
    object_n$window_map <- (object_n$window_E -
        start(object_n))/object_n$window_length

    object <- c(object_p, object_n)

    ## Give the no mapped peaks a value 3 for their position
    object$window_map[object$window_map == -Inf | object$window_map == Inf] <- 3

    return(object)
}

## Plot
.windowPlot <- function(df, title)
{
    p1 <- ggplot(df, aes(x = window_map)) +
        geom_density(adjust = 0.2, color = "Orange") +
        theme_bw() + xlab("Given Region") +
        scale_x_continuous(breaks = c(0,1),
        labels = c("5'", "3'")) +
        ggtitle(title) + ylab("Density of Peaks")
    return(p1)
}

.windowPlotGroup <- function(df, title)
{
    p1 <- ggplot(df, aes(x = window_map, color=groupIn)) +
        geom_density(adjust = 0.2) +
        theme_bw() + xlab("Given Region") +
        scale_x_continuous(breaks = c(0,1),
        labels = c("5'", "3'")) +
        ggtitle(title) + ylab("Density of Peaks")
    return(p1)
}
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## The "windowProfile" methods for GRanges objects.
##

#' @title windowProfile for the GRanges objects
#'
#' @description An function to check the position of peaks within the given
#'              GRanges windows.
#'
#' @author You Zhou, Kathi Zarnack
#'
#' @param object A GRanges object which contains all the peaks that you
#'                want to check
#' @param annotation A GRanges object that includes the customised genomic
#'                   region.
#' @param title The main title for the output meta gene profile plot.
#' @param group The column name which contains the information of grouping
#'     for making the comparison plot. NA means all the peaks belongs to
#'     the same catagory.
#' @param nomap A logical vector (TRUE or FALSE). It indicates whether you
#'              would like to exclude peaks that cannot assign to annotations
#'              in the plot.
#' @details
#' \itemize{
#'     Here is an explanation of output meta data in the \code{list 1}:
#'     \item \code{center}: The center position of each peaks. This center
#'     position is used for calculating the position of peaks within the
#'     genomic regions.
#'     \item \code{window_S} and \code{window_E}: The boundary of the
#'     annotation that peaks are assigned.
#'     \item \code{window_length}: The length of the annotation feature that
#'     peak assigned.
#'     \item \code{window_map}: The relative position of each peak. This value
#'     close to 0 means this peak located close to the 5' end of the
#'     annotation. The position value close to one means the peak close to
#'     the 3' end. Value 3 means this peaks can not map to any annotation.
#' }
#'
#' @return A list object, the list 1 contains the information of the
#'         assignment of the peaks and their position value within the given
#'         region. The value close to 1 means the peak close to the end of
#'         region in 3' end direction. The list 2 includes the ggplot of
#'         windowProfile.
#' @examples
#' ## Load the test data and get the path to the test gff3 file
#' testpath <- system.file("extdata", package = "cliProfiler")
#' test <- readRDS(file.path(testpath, "test.rds"))
#' test_gff3 <- file.path(testpath, "annotation_test.gff3")
#' test_gff3 <- rtracklayer::import.gff3(test_gff3)
#'
#' output <- windowProfile(test, test_gff3)
#' @export

windowProfile <-
    function(object, annotation, title="Window Profile", group=NA, nomap=FALSE)
    {
        if (missing(object))
            stop("The input GRanges object is missing.")
        if (missing(annotation))
            stop("The GRanges object of annotation file is missing.")
        if (!isS4(object))
            stop("The input object should be a GRanges object.")
        if (!isS4(annotation))
            stop("The annotation should be a GRanges object.")
        if (!is.logical(nomap))
            stop("The nomap should be a logical vector (TRUE or FALSE).")
        if (!is.na(group) &
            (group %in% colnames(elementMetadata(object)) == FALSE))
            stop("When the option group is used, please make sure that
            the group should be a column name of input GRanges object which
            contains the information of which group this peaks belongs to
            and will be used for seperate the for generating the meta gene
            profile curve for different groups")

        ## Center peaks and getting windows
        anno <- annotation

        anno <- anno[order(width(anno), decreasing = TRUE)]

        object <- .centerPeaks(object)

        o <- findOverlaps(object, anno, select = "first")

        object$window_S <-  start(anno)[o]
        object$window_E <-  end(anno)[o]
        object$window_length <-  width(anno)[o]

        object$window_S[is.na(object$window_S)] <- 0
        object$window_E[is.na(object$window_E)] <- 0
        object$window_length[is.na(object$window_length)] <- 0

        object <- .windowPosition(object)

        GenomicRanges::start(object) <- object$oriStart
        GenomicRanges::end(object) <- object$oriEnd
        object$oriStart <- NULL
        object$oriEnd <- NULL

        df <- as.data.frame(object)
        if (is.na(group)) {
            if (nomap == FALSE) {
                df <- df[df$window_map != 3,]
                p1 <- .windowPlot(df, title) +
                    coord_cartesian(xlim = c(0,1))
            }
            if (nomap == TRUE){
                p1 <- .windowPlot(df, title) +
                    coord_cartesian(xlim = c(0,3))
            }
        }
        if (!is.na(group)) {
            df$groupIn <- df[,colnames(df, title) == group]
            if (nomap == FALSE) {
                df <- df[df$window_map != 3,]
                p1 <- .windowPlotGroup(df) +
                    coord_cartesian(xlim = c(0,1))
            }
            if (nomap == TRUE){
                p1 <- .windowPlotGroup(df, title) +
                    coord_cartesian(xlim = c(0,3))
            }
        }
        output <- list(Peaks = object, Plot = p1)
        return(output)
    }
