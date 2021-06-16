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

#' @rdname windowProfile
setMethod("windowProfile", signature(object="GRanges", annotation="GRanges"),
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
)
