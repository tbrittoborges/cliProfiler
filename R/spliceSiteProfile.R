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

#' @title spliceSiteProfile for the GRanges objects
#'
#' @description An function to check the enrichment of peaks around the splice
#'              sites in a absolute distance.
#'
#' @author You Zhou, Kathi Zarnack
#'
#' @param object A GRanges object which should contains all the peaks that you
#'                want to check
#' @param annotation A path way to the annotation file. The format of the
#'                   annotation file should be gff3 and downloaded from
#'                   https://www.gencodegenes.org/
#' @param title The main title for the output meta gene profile plot.
#' @param exlevel A parameter for the annotation filtering. exlevel represents
#'     the level that you would like to exclude. NA means no level filtering
#'     for the annotation file. The level from the annotations refers to
#'     how reliable this annotation is. For more information about level
#'     please check
#'     https://www.gencodegenes.org/pages/data_format.html.
#' @param extranscript_support_level A parameter for the annotation filtering.
#'     extranscript_support_level represents the transcript_support_level
#'     that you would like to exclude (e.g. 4 and 5). NA means no
#'     transcript_support_level filtering for the annotation file.
#'     Transcripts are scored according to how well mRNA and EST alignments
#'     match over its full length. Here the number 6 means the
#'     transcript_support_level NA. For more information about level please
#'     check
#'     https://www.gencodegenes.org/pages/data_format.html.
#' @param exon_length_filtering The exon_length_filtering should be a logical
#'     value which indicated whether user would like to exclude the exons
#'     that have a length less than flanking value. Set this parameter to TRUE
#'     to turn on this filtering step.
#' @param intron_length_filtering The intron_length_filtering should be a
#'     logical value which indicated whether user would like to exclude the
#'     introns that have a length less than flanking value. Set this parameter
#'     to TRUE to turn on this filtering step.
#' @param flanking The size of the flanking windows that you would like to
#'                 check. Flanking=5 will give you the result of the 10+1nt
#'                 windows around the center of peaks.
#' @param bin A number that indicates how many bins would you like to use in
#'            the histogram.
#'
#' @return A list object, the list 1 contains the information of the
#'         position of peaks around 5' or 3' splice sites. The list 2
#'         includes the plot of spliceSiteProfile
#' @examples
#' ## Load the test data and get the path to the test gff3 file
#' testpath <- system.file("extdata", package = "cliProfiler")
#' test <- readRDS(file.path(testpath, "test.rds"))
#' test_gff3 <- file.path(testpath, "annotation_test.gff3")
#'
#' output <- spliceSiteProfile(test, test_gff3,
#'   flanking = 200, bin = 40
#' )
#' @export

spliceSiteProfile <- function(object, annotation, title="Splice Site Profile",
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
        if (sum(!exlevel %in% level) > 0 |
            sum(!extranscript_support_level %in% tsl) > 0)
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
        if (sum(is.na(exlevel)) == 0) {
            introns <- .annoFilterLevel(introns, exlevel)
        }
        if (sum(is.na(extranscript_support_level)) == 0) {
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
        ## Filter out the duplicates

        start <- c(start_p, start_n)
        start <- unique(start)

        start <- start[order(start$level,
            start$transcript_support_level,
            -start$transcript_length)]

        ## Extend the 5'SS for checking
        start <- start+flanking
        o <- findOverlaps(object, start, select = "first")

        object$start_2 <- start$ss[o]

        ## Check 3' splice sites (3'SS)
        end_p <- .buildGRend(introns_p)
        end_n <- .buildGRstart(introns_n)
        end <- c(end_p, end_n)
        end <- unique(end)

        end <- end[order(end$level,
            end$transcript_support_level,
            -end$transcript_length)]

        end <- end+flanking
        oo <- findOverlaps(object, end, select = "first")

        object$end_2 <- end$ss[oo]

        object_p <- object[strand(object)=="+"]
        object_n <- object[strand(object)=="-"]

        ## calculation
        object_p$Position5SS <- start(object_p)-object_p$start_2
        object_n$Position5SS <- object_n$start_2-start(object_n)

        object_p$Position3SS <- start(object_p)-object_p$end_2
        object_n$Position3SS <- object_n$end_2-start(object_n)

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

        object$start_2 <- NULL
        object$end_2 <- NULL
        GenomicRanges::start(object) <- object$oriStart
        GenomicRanges::end(object) <- object$oriEnd
        object$oriStart <- NULL
        object$oriEnd <- NULL

        out <- list(Peaks=object, Plot=p)
        return(out)
    }
