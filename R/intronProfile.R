## ============================================================================
## The intronProfile function for GRanges objects.
## ----------------------------------------------------------------------------
#' @import dplyr
#' @import methods
#' @import ggplot2
#' @importFrom S4Vectors split
#' @importFrom GenomicRanges psetdiff
#' @importFrom GenomicRanges width
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges strand
#' @importFrom GenomicRanges end

## ============================================================================
## Small functions
## ----------------------------------------------------------------------------

.getIntrons <- function(anno)
{
    trans <- anno[anno$type=="transcript"]
    exon.gr <- anno[anno$type=="exon"]
    exon.gr$exonID <- paste0(exon.gr$transcript_id, ",", exon.gr$exon_number)
    exonsByTranscript <- split(anno[anno$type=="exon"],
        anno[anno$type=="exon"]$transcript_id)
    trans <- trans[match(names(exonsByTranscript),
        trans$transcript_id)]

    introns <- psetdiff(trans, exonsByTranscript)
    introns <- as(introns, "GRangesList") %>% unlist

    introns$transcript_id <- names(introns)
    names(introns) <- NULL

    ## make intron number
    introns.n <- introns[strand(introns) == "-"] %>%
        sort(., decreasing = TRUE) %>% as.data.frame() %>%
        group_by(transcript_id) %>%
        mutate(intron_number = seq_len(n())) %>%
        makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
    introns.p <- introns[strand(introns) == "+"] %>%
        sort(.) %>% as.data.frame() %>% group_by(transcript_id) %>%
        mutate(intron_number = seq_len(n())) %>%
        makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
    introns <- c(introns.n, introns.p)

    introns$intronID1 <- paste0(introns$transcript_id,
        ",",introns$intron_number)
    introns$intronID2 <- paste0(introns$transcript_id,
        ",",introns$intron_number+1)
    ## Assign the length of flanking exons this value will be used in ssprofile
    introns$exon1_l <- width(exon.gr)[match(introns$intronID1,
        exon.gr$exonID)]
    introns$exon2_l <- width(exon.gr)[match(introns$intronID2,
        exon.gr$exonID)]

    introns$transcript_length <-
        width(trans)[match(introns$transcript_id,
        trans$transcript_id)]
    introns$level <-
        trans$level[match(introns$transcript_id,
        trans$transcript_id)]

    ## tsl means transcript_support_level
    introns$tsl <-
        trans$transcript_support_level[match(introns$transcript_id,
        trans$transcript_id)]

    ## Remove NA avoid warning message
    introns$level[is.na(introns$level)] <- 6
    introns$tsl[is.na(introns$tsl)] <- 6
    introns$tsl[introns$tsl == "NA"] <-6

    ## Change the type of data for the ranking step
    introns$level <- as.integer(introns$level)
    introns$transcript_support_level <- as.integer(introns$tsl)
    introns$tsl <- NULL

    return(introns)
}

.intronPosition <- function(object)
{
    ##-----calculate the position-----##
    object_p <- object[strand(object) == "+"]
    object_p$Intron_map <- (start(object_p) -
        object_p$Intron_S)/object_p$Intron_length
    object_n <- object[strand(object) == "-"]
    object_n$Intron_map <- (object_n$Intron_E -
        start(object_n))/object_n$Intron_length

    object <- c(object_p, object_n)

    ## Give the no mapped peaks a value 3 for their position
    object$Intron_map[object$Intron_map == -Inf|object$Intron_map == Inf] <- 3

    return(object)
}

## Plot
.intronPlot <- function(df, title)
{
    p1 <- ggplot(df, aes(x = Intron_map)) +
        geom_density(adjust = 0.2, color = "Orange") +
        theme_bw() + xlab("Intron") +
        scale_x_continuous(breaks = c(0,1),
        labels = c("5'SS", "3'SS")) +
        ggtitle(title) + ylab("Density of Peaks")
    return(p1)
}

.intronPlotGroup <- function(df, title)
{
    p1 <- ggplot(df, aes(x = Intron_map, color = groupIn)) +
        geom_density(adjust = 0.2) +
        theme_bw() + xlab("Intron") +
        scale_x_continuous(breaks = c(0,1),
        labels = c("5'SS", "3'SS")) +
        ggtitle(title) + ylab("Density of Peaks")
    return(p1)
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## The "intronProfile" methods for GRanges objects.
##

#' @title intronProfile for the GRanges objects
#'
#' @description An function to check the position of peaks in the intronic
#'              region.
#'
#' @author You Zhou, Kathi Zarnack
#'
#' @param object A GRanges object which should contains all the peaks that you
#'                want to check
#' @param annotation A path way to the annotation file. The format of the
#'                   annotation file should be gff3 and downloaded from
#'                   https://www.gencodegenes.org/
#' @param title The main title for the output meta gene profile plot.
#' @param group The column name which contains the information of grouping
#'     for making the comparison plot. NA means all the peaks belongs to
#'     the same catagory.
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
#' @param maxLength A numeric value which indicate the maximum value of exon
#'                  length for the annotation filtering. Or a NA which will
#'                  turn off the max length annotation filtering.
#' @param minLength A numeric value which indicate the minimum value of exon
#'                  length for the annotation filtering. Or a NA which will
#'                  turn off the min length annotation filtering.
#' @param nomap A logical vector (TRUE or FALSE). It indicates whether you
#'              would like to exclude peaks that cannot assign to annotations
#'              in the plot.
#' @details
#' \itemize{
#'     Here is an explanation of output meta data in the \code{list 1}:
#'     \item \code{center}: The center position of each peaks. This center
#'     position is used for calculating the position of peaks within the
#'     genomic regions.
#'     \item \code{Intron_S} and \code{Intron_E}: The location of 5' and 3'
#'     splice sites (SS) of the intron.
#'     \item \code{Intron_length}: The length of the intron that peak assigned.
#'     \item \code{Intron_transcript_id}: The transcript ID for the intron.
#'     \item \code{Intron_map}: The relative position of each peak. This value
#'     close to 0 means this peak located close to the 5' SS. The position
#'     value close to one means the peak close to the 3' SS. Value 3 means
#'     this peaks can not map to any annotation.
#' }
#'
#' @return A list object, the list 1 contains the information of the
#'     assignment of the peaks and their position value within the intron.
#'     The value close to 1 means the peak close to the 3' splice site.
#'     The list 2 includes the plot of intronProfile.
#' @examples
#' ## Load the test data and get the path to the test gff3 file
#' testpath <- system.file("extdata", package = "cliProfiler")
#' test <- readRDS(file.path(testpath, "test.rds"))
#' test_gff3 <- file.path(testpath, "annotation_test.gff3")
#'
#' output <- intronProfile(test, test_gff3)
#' @export
#'

intronProfile <- function(object, annotation, title="Intron Profile", group=NA,
    exlevel=NA, extranscript_support_level=NA, maxLength=NA, minLength=NA,
    nomap=FALSE)
    {
        if (missing(object))
            stop("The input GRanges object is missing.")
        if (missing(annotation))
            stop("The path to the gff3 annotation file is missing.")
        if (!isS4(object))
            stop("The input object should be a GRanges object.")
        if (!is.logical(nomap))
            stop("The nomap should be a logical vector (TRUE or FALSE).")
        level <- c(1,2,3,NA)
        tsl <- c(1,2,3,4,5,6,NA)
        if (sum(!exlevel %in% level) > 0 |
            sum(!extranscript_support_level %in% tsl) > 0)
            warning("The exlevel should be a vector includes the value of 1, 2,
                3 or NA. extranscript_support_level should be a vector includes
                value 1, 2, 3, 4, 5, 6 or NA.")
        if (!is.numeric(maxLength)&!is.na(maxLength))
            stop("The maxLength should be a numeric value which indicate the
                maximum length of the exon or NA.")
        if (!is.numeric(minLength)&!is.na(minLength))
            stop("The minLength should be a numeric value which indicate the
                minimum length of the exon or NA.")
        if (!is.na(group) &
            (group %in% colnames(elementMetadata(object)) == FALSE))
            stop("When the option group is used, please make sure that
            the group should be a column name of input GRanges object which
            contains the information of which group this peaks belongs to
            and will be used for seperate the for generating the meta gene
            profile curve for different groups")

        ## Center peaks and getting introns
        anno <- rtracklayer::import.gff3(con = annotation)
        object <- .centerPeaks(object)
        introns <- .getIntrons(anno)
        introns$intronID <- seq_len(length(introns))

        ## Annotation filtering
        if (sum(is.na(exlevel)) == 0) {
            introns <- .annoFilterLevel(introns, exlevel)
        }
        if (sum(is.na(extranscript_support_level)) == 0) {
            introns <- .annoFilterTSL(introns, extranscript_support_level)
        }
        if (!is.na(maxLength)) {
            introns <- .lenFilteringMax(introns, maxLength)
        }
        if (!is.na(minLength)) {
            introns <- .lenFilteringMin(introns, minLength)
        }

        ## Ranking intron for the peaks assignment
        introns <- introns[order(introns$level,
            introns$transcript_support_level, -introns$transcript_length)]

        o2 <- findOverlaps(object, introns, select = "first")

        ## Assign the annotation details to the input GRanges
        object$Intron_S <-  start(introns)[o2]
        object$Intron_E <-  end(introns)[o2]
        object$Intron_length <-  width(introns)[o2]
        object$intronID <- introns$intronID[o2]
        object$Intron_transcript_id <- introns$transcript_id[o2]

        object$Intron_S[is.na(object$Intron_S)] <- 0
        object$Intron_E[is.na(object$Intron_E)] <- 0
        object$Intron_length[is.na(object$Intron_length)] <- 0
        object$intronID[is.na(object$intronID)] <- "NO"
        object$Intron_transcript_id[is.na(object$Intron_transcript_id)]<-
            "NO"

        object <- .intronPosition(object)
        object$intronID <- NULL

        GenomicRanges::start(object) <- object$oriStart
        GenomicRanges::end(object) <- object$oriEnd
        object$oriStart <- NULL
        object$oriEnd <- NULL

        df <- as.data.frame(object)
        if (is.na(group)) {
            if (nomap == FALSE) {
                df <- df[df$Intron_map != 3,]
                p1 <- .intronPlot(df, title) +
                    coord_cartesian(xlim = c(0,1))
            }
            if (nomap == TRUE){
                p1 <- .intronPlot(df, title) +
                    coord_cartesian(xlim = c(0,3))
            }
        }
        if (!is.na(group)) {
            df$groupIn <- df[,colnames(df) == group]
            if (nomap == FALSE) {
                df <- df[df$Intron_map != 3,]
                p1 <- .intronPlotGroup(df, title) +
                    coord_cartesian(xlim = c(0,1))
            }
            if (nomap == TRUE){
                p1 <- .intronPlotGroup(df, title) +
                    coord_cartesian(xlim = c(0,3))
            }
        }

        output <- list(Peaks = object, Plot = p1)
        return(output)
    }





