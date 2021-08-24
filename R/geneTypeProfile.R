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
        theme(axis.text.x = element_text(angle = 45, hjust=0))
    return(p1)
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## The "geneTypeProfile" methods for GRanges objects.
##

#' @title geneTypeProfile for the GRanges objects
#'
#' @description An function to check the gene type belonging for the peaks.
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
#' @details
#' \itemize{
#'     Here is an explanation of output meta data in the \code{list 1}:
#'     \item \code{center}: The center position of each peaks. This center
#'     position is used for calculating the position of peaks within the
#'     genomic regions.
#'     \item \code{geneType}: The gene type of the gene that input peak belongs
#'     to.
#'     \item \code{Gene_ID}: The gene ID of the gene that input peak
#'     belongs to.
#' }
#'
#' @return A list object, the list 1 contains the information of the
#'         assignment of the peaks and the gene type of their located genes.
#'         The list 2 includes the plot of geneTypeProfile
#' @examples
#' ## Load the test data and get the path to the test gff3 file
#' testpath <- system.file("extdata", package = "cliProfiler")
#' test <- readRDS(file.path(testpath, "test.rds"))
#' test_gff3 <- file.path(testpath, "annotation_test.gff3")
#'
#' output <- geneTypeProfile(test, test_gff3)
#' @export

geneTypeProfile <- function(object, annotation, title="Gene Type Profile",
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
        anno_exon <- anno[anno$type == "exon"]
        anno_exon <- as.data.frame(anno_exon) %>%
            group_by(transcript_id) %>%
            mutate(trans_len = sum(width))  %>%
            filter(exon_number == "1") %>%
            makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

        trans$transcript_length <- anno_exon$trans_len[
            match(trans$transcript_id, anno_exon$transcript_id)]

        ## annotation filtering
        if (sum(is.na(exlevel)) == 0) {
            trans <- .annoFilterLevel(trans, exlevel)
        }
        if (sum(is.na(extranscript_support_level)) == 0) {
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
