## ============================================================================
## The exonProfile function for GRanges objects.
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

## Length filtering
.lenFilteringMax <- function(anno, maxLength)
{
    anno <- anno[width(anno) <= maxLength]
    return(anno)
}
.lenFilteringMin <- function(anno, minLength)
{
    anno <- anno[width(anno) >= minLength]
    return(anno)
}

## Calculation of exon position
.exonPosition <- function(object)
{
    ##-----calculate the position-----##
    object_p <- object[strand(object) == "+"]
    object_p$exon_map <- (start(object_p) -
                                object_p$exon_S)/object_p$exon_length
    object_n <- object[strand(object) == "-"]
    object_n$exon_map <- (object_n$exon_E -
                                start(object_n))/object_n$exon_length

    object <- c(object_p, object_n)

    ## Give the no mapped peaks a value 3 for their position
    object$exon_map[object$exon_map == -Inf | object$exon_map == Inf] <- 3

    return(object)
}

## Exclude the first and last exon for each transcript
.exonExtract <- function(anno)
{
    anno <- as.data.frame(anno) %>% group_by(transcript_id) %>%
        mutate(maxExon = max(exon_number))
    anno <- makeGRangesFromDataFrame(anno, keep.extra.columns = TRUE)
    anno <- anno[anno$exon_number != 1 & anno$exon_number != anno$maxExon]
    return(anno)
}

## plot
.exonPlot <- function(df, title)
{
    p1 <- ggplot(df, aes(x = exon_map)) +
        geom_density(adjust = 0.2, color = "Orange") +
        theme_bw() + xlab("Exon") +
        scale_x_continuous(breaks = c(0,1),
        labels = c("3'SS", "5'SS")) +
        ggtitle(title) + ylab("Density of Peaks")
    return(p1)
}

.exonPlotGroup <- function(df, title)
{
    p1 <- ggplot(df, aes(x = exon_map, color = groupIn)) +
        geom_density(adjust = 0.2) +
        theme_bw() + xlab("Exon") +
        scale_x_continuous(breaks = c(0,1),
        labels = c("3'SS", "5'SS")) +
        ggtitle(title) + ylab("Density of Peaks")
    return(p1)
}
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## The "exonProfile" methods for GRanges objects.
##

#' @rdname exonProfile
setMethod("exonProfile", signature(object="GRanges"),
    function(object, annotation, title="Exon Profile", group=NA,
    exlevel=NA, extranscript_support_level=NA, maxLength=NA,minLength=NA,
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
        if (!exlevel %in% level| !extranscript_support_level %in% tsl)
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

        ## Center peaks and getting exons
        anno <- rtracklayer::import.gff3(con = annotation)
        anno_exon <- anno[anno$type == "exon"]
        trans <- anno[anno$type=="transcript"]
        anno_exon$transcript_length <-
            GenomicRanges::width(trans)[match(anno_exon$transcript_id,
            trans$transcript_id)]

        ## Give 6 to the NAs in the transcript support level
        anno_exon$transcript_support_level[
            is.na(anno_exon$transcript_support_level)] <- 6
        anno_exon$transcript_support_level[
            anno_exon$transcript_support_level == "NA"] <- 6

        object <- .centerPeaks(object)

        ## Annotation filtering
        if (!is.na(exlevel)) {
            anno_exon <- .annoFilterLevel(anno_exon, exlevel)
        }
        if (!is.na(extranscript_support_level)) {
            anno_exon <- .annoFilterTSL(anno_exon, extranscript_support_level)
        }
        if (!is.na(maxLength)) {
            anno_exon <- .lenFilteringMax(anno_exon, maxLength)
        }
        if (!is.na(minLength)) {
            anno_exon <- .lenFilteringMin(anno_exon, minLength)
        }

        anno_exon <- .exonExtract(anno_exon)

        anno_exon <- anno_exon[order(anno_exon$level,
            anno_exon$transcript_support_level,
            -anno_exon$transcript_length)]
        o <- findOverlaps(object, anno_exon, select = "first")

        object$exon_S <-  start(anno_exon)[o]
        object$exon_E <-  end(anno_exon)[o]
        object$exon_length <-  width(anno_exon)[o]
        object$exon_transcript_id <- anno_exon$transcript_id[o]

        object$exon_S[is.na(object$exon_S)] <- 0
        object$exon_E[is.na(object$exon_E)] <- 0
        object$exon_length[is.na(object$exon_length)] <- 0
        object$exon_transcript_id[is.na(object$exon_transcript_id)]<-"NO"

        object <- .exonPosition(object)

        GenomicRanges::start(object) <- object$oriStart
        GenomicRanges::end(object) <- object$oriEnd
        object$oriStart <- NULL
        object$oriEnd <- NULL
        df <- as.data.frame(object)

        if (is.na(group)) {
            if (nomap == FALSE) {
                df <- df[df$exon_map != 3,]
                p1 <- .exonPlot(df, title) + coord_cartesian(xlim = c(0,1))
            }
            if (nomap == TRUE){
                p1 <- .exonPlot(df, title) + coord_cartesian(xlim = c(0,3))
            }
        }
        if (!is.na(group)) {
            df$groupIn <- df[,colnames(df) == group]
            if (nomap == FALSE) {
                df <- df[df$exon_map != 3,]
                p1 <- .exonPlotGroup(df, title) + coord_cartesian(xlim = c(0,1))
            }
            if (nomap == TRUE){
                p1 <- .exonPlotGroup(df, title) + coord_cartesian(xlim = c(0,3))
            }
        }

        output <- list(Peaks = object, Plot = p1)
        return(output)
    }
)
