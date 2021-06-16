## ============================================================================
## The geneTypeProfile function for GRanges objects.
## ----------------------------------------------------------------------------

#' @import dplyr
#' @import methods
#' @import ggplot2
#' @importFrom GenomicRanges GRanges


## ============================================================================
## Small functions
## ----------------------------------------------------------------------------

.geneTypePlot <- function(df, title)
{
    p1 <- ggplot(df, aes(x=geneType)) +
        geom_bar(fill="orange", color="cornflowerblue") +
        geom_label(stat='count', aes(label=..count..), vjust=0) +
        theme_bw() + ggtitle(title) +
        theme(axis.text.x = element_text(angle = 45, hjust=1))
    return(p1)
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## The "geneTypeProfile" methods for GRanges objects.
##

#' @rdname geneTypeProfile
setMethod("geneTypeProfile", signature(object="GRanges"),
    function(object, annotation, title="Gene Type Profile",
    exlevel=NA, extranscript_support_level=NA)
    {
        if (missing(object))
            stop("The input GRanges object is missing.")
        if (missing(annotation))
            stop("The path to the gff3 annotation file is missing.")
        if (!isS4(object))
            stop("The input object should be a GRanges object.")

        ## Center peaks and getting transcripts
        object <- .centerPeaks(object)

        anno <- rtracklayer::import.gff3(con = annotation)
        trans <- anno[anno$type=="transcript"]
        trans$transcript_support_level[
            is.na(trans$transcript_support_level)] <- 6
        trans$transcript_support_level[
            trans$transcript_support_level == "NA"] <- 6
        trans$transcript_length <- GenomicRanges::width(trans)

        ## annotation filtering
        if (!is.na(exlevel)) {
            trans <- .annoFilterLevel(trans, exlevel)
        }
        if (!is.na(extranscript_support_level)) {
            trans <- .annoFilterTSL(trans, extranscript_support_level)
        }

        ## Peaks assignment
        trans <- trans[order(trans$level,
            trans$transcript_support_level, -trans$transcript_length)]

        o <- findOverlaps(object, trans, select = "first")
        object$geneType <- trans$gene_type[o]
        object$Gene_ID <- trans$gene_id[o]
        df <- as.data.frame(object)

        p1 <- .geneTypePlot(df, title)

        GenomicRanges::start(object) <- object$oriStart
        GenomicRanges::end(object) <- object$oriEnd
        object$oriStart <- NULL
        object$oriEnd <- NULL

        out <- list(Peaks = object, Plot = p1)
        return(out)
    }
)
