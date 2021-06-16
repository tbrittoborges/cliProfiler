## ============================================================================
## The spliceSiteProfile function for GRanges objects.
## ----------------------------------------------------------------------------

#' @import dplyr
#' @import methods
#' @import ggplot2
#' @importFrom GenomicRanges width
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges strand
#' @importFrom GenomicRanges end
#' @importFrom GenomicRanges findOverlaps
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicRanges GRanges
#' @importFrom rtracklayer import.gff3
#' @importFrom S4Vectors elementMetadata

## ============================================================================
## Small functions
## ----------------------------------------------------------------------------

## build the GRange objects of splice site from annotation file
.buildGRend <- function(anno)
{
    anno <- data.frame(seqnames=GenomicRanges::seqnames(anno),
        start=GenomicRanges::end(anno),
        end=GenomicRanges::end(anno),
        strand=GenomicRanges::strand(anno),
        transcript_support_level=anno$transcript_support_level,
        level=anno$level,
        transcript_length=anno$transcript_length,
        ss=GenomicRanges::end(anno))
    anno <- makeGRangesFromDataFrame(anno,keep.extra.columns=TRUE)
    return(anno)
}

.buildGRstart <- function(anno)
{
    anno <- data.frame(seqnames=GenomicRanges::seqnames(anno),
        start=GenomicRanges::start(anno),
        end=GenomicRanges::start(anno),
        strand=GenomicRanges::strand(anno),
        transcript_support_level=anno$transcript_support_level,
        level=anno$level,
        transcript_length=anno$transcript_length,
        ss=GenomicRanges::start(anno))
        anno <- makeGRangesFromDataFrame(anno, keep.extra.columns=TRUE)
    return(anno)
}
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## The "spliceSiteProfile" methods for GRanges objects.
##

#' @rdname spliceSiteProfile
setMethod("spliceSiteProfile", signature(object="GRanges"),
    function(object, annotation, title="Splice Site Profile",
    exlevel=NA, extranscript_support_level=NA, exon_length_filtering=TRUE,
    intron_length_filtering=TRUE, flanking=150, bin=30)
    {
        if (missing(object))
            stop("The input GRanges object is missing.")
        if (missing(annotation))
            stop("The path to the gff3 annotation file is missing.")
        if (!isS4(object))
            stop("The input object should be a GRanges object.")
        if (!is.numeric(flanking))
            stop("The flanking should be a numeric.")
        if (!is.numeric(bin))
            stop("The parameter bin should be a numeric.")
        level <- c(1,2,3,NA)
        tsl <- c(1,2,3,4,5,6,NA)
        if (!exlevel %in% level| !extranscript_support_level %in% tsl)
            warning("The exlevel should be a vector includes the value of 1, 2,
                3 or NA. extranscript_support_level should be a vector includes
                value 1, 2, 3, 4, 5, 6 or NA.")
        if (!is.logical(exon_length_filtering))
            stop("The exon_length_filtering should be a logical value.")
        if (!is.logical(intron_length_filtering))
            stop("The intron_length_filtering should be a logical value.")

        ## Center peaks and getting exons
        anno <- rtracklayer::import.gff3(con = annotation)
        trans <- anno[anno$type=="transcript"]
        introns <- .getIntrons(anno)

        ## Give 6 to the NAs in the transcript support level
        introns$transcript_support_level[
            is.na(introns$transcript_support_level)] <- 6
        introns$transcript_support_level[
            introns$transcript_support_level == "NA"] <- 6
        introns$transcript_length <-
            GenomicRanges::width(trans)[match(introns$transcript_id,
            trans$transcript_id)]
        ## center peaks
        object <- .centerPeaks(object)
        ## annotation filtering
        if (!is.na(exlevel)) {
            introns <- .annoFilterLevel(introns, exlevel)
        }
        if (!is.na(extranscript_support_level)) {
            introns <- .annoFilterTSL(introns, extranscript_support_level)
        }
        if (isTRUE(exon_length_filtering)) {
            introns <- introns[introns$exon1_l >= flanking]
            introns <- introns[introns$exon2_l >= flanking]
        }
        if (isTRUE(intron_length_filtering)) {
            introns <- introns[width(introns) >= flanking]
        }

        ## split positive and negative strand
        introns_p <- introns[strand(introns) == "+"]
        introns_n <- introns[strand(introns) == "-"]

        ## Check the 5' splice site (5'SS)
        ## First generate the coordinates for the 5'SS
        start_p <- .buildGRstart(introns_p)
        start_n <- .buildGRend(introns_n)

        start_p <- start_p[order(start_p$level,
            start_p$transcript_support_level,
            -start_p$transcript_length)]
        start_n <- start_n[order(start_n$level,
            start_n$transcript_support_level,
            -start_n$transcript_length)]

        ## Extend the 5'SS for checking
        start_p <- start_p+flanking
        o <- findOverlaps(object, start_p, select = "first")
        start_n <- start_n+flanking
        o2 <- findOverlaps(object, start_n, select = "first")

        object$start_p <- start_p$ss[o]
        object$start_n <- start_n$ss[o2]

        object_p <- object[strand(object)=="+"]
        object_n <- object[strand(object)=="-"]

        object_p$Position5SS <- start(object_p)-object_p$start_p
        object_n$Position5SS <- object_n$start_n-start(object_n)

        object <- c(object_n, object_p)

        ## Check 3' splice sites (3'SS)
        end_p <- .buildGRend(introns_p)
        end_n <- .buildGRstart(introns_n)

        end_p <- end_p[order(end_p$level,
            end_p$transcript_support_level, -end_p$transcript_length)]
        end_n <- end_n[order(end_n$level,
            end_n$transcript_support_level, -end_n$transcript_length)]

        end_p <- end_p+flanking
        oo <- findOverlaps(object, end_p, select = "first")
        end_n <- end_n+flanking
        oo2 <- findOverlaps(object, end_n, select = "first")

        object$end_p <- end_p$ss[oo]
        object$end_n <- end_n$ss[oo2]

        object_p <- object[strand(object)=="+"]
        object_n <- object[strand(object)=="-"]

        object_p$Position3SS <- start(object_p)-object_p$end_p
        object_n$Position3SS <- object_n$end_n-start(object_n)

        object <- c(object_n, object_p)

        ## Make dataframe without NAs for the ggplot
        df_5SS <- data.frame(Position=
            object$Position5SS[!is.na(object$Position5SS)],
            Type="5' Splice Site")
        df_3SS <- data.frame(Position=
            object$Position3SS[!is.na(object$Position3SS)],
            Type="3' Splice Site")

        df <- rbind(df_5SS, df_3SS)

        df$Type <- factor(df$Type, levels = c("5' Splice Site",
            "3' Splice Site"))
        p <- ggplot(df, aes(x = Position)) +
            geom_histogram(bins = bin, fill = "orange", color = "grey") +
            theme_bw() + facet_wrap(~Type) + ggtitle(title) +
            geom_vline(xintercept = 0, linetype = 2,
            color = "cornflowerblue")

        object$start_p <- NULL
        object$start_n <- NULL
        object$end_p <- NULL
        object$end_n <- NULL
        GenomicRanges::start(object) <- object$oriStart
        GenomicRanges::end(object) <- object$oriEnd
        object$oriStart <- NULL
        object$oriEnd <- NULL

        out <- list(Peaks=object, Plot=p)
        return(out)
    }
)
