## ============================================================================
## The metaGeneProfile function for GRanges objects.
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
#' @importFrom utils globalVariables


## To fix the global variable note
utils::globalVariables(c("transcript_id", "geneType", "Fraction", "Number",
    "groupIn", "exon_map", "Intron_map", "..count..", "exlevel", "maxLength",
    "extranscript_support_level", "level","extranscript_support_level",
    "minLength", "window_map","exon_number", "location", "title","."))

## Extract the center position of all input peaks
.centerPeaks <- function(granges)
{
    ## Store the original information for the peaks
    granges$oriStart <- start(granges)
    granges$oriEnd <- end(granges)

    ## Extract the center of all peaks
    granges$half_length <- width(granges)/2
    ## For the peaks which width = 1, keep them as what they are
    granges$half_length[granges$half_length == 0.5] <- 0
    granges$half_length <- round(granges$half_length)
    GenomicRanges::start(granges) <- GenomicRanges::start(granges) +
        granges$half_length
    GenomicRanges::end(granges) <- start(granges)
    granges$half_length <- NULL

    granges$center <- start(granges)
    return(granges)
}

## Function for filtering annotations
.annoFilterLevel <- function(anno, exlevel = exlevel)
{
    anno <- anno[!anno$level %in% exlevel]
    return(anno)
}

.annoFilterTSL <- function(anno, extsl = extranscript_support_level)
{
    anno <- anno[!anno$transcript_support_level %in% extsl]
    return(anno)
}

## Extract the position of the each annotation fragment (positive strand)
## without the introns.
.annoCalculateP <- function(annoP){
    annoP <- as.data.frame(annoP)
    annoP <- annoP %>% group_by(transcript_id) %>%
        mutate(full_length = sum(width)) %>%
        mutate(Rstart = min(start)) %>% mutate(Rend = max(end))
    annoP <- makeGRangesFromDataFrame(annoP,keep.extra.columns = TRUE)
    annoP <- sort(annoP) %>% as.data.frame()
    annoP <- annoP %>% group_by(transcript_id) %>%
        mutate(Add = c(0, cumsum(width)[-length(width)]))
    annoP <- makeGRangesFromDataFrame(annoP,keep.extra.columns = TRUE)
    return(annoP)
}

## Extract the position of the each annotation fragment (negative strand)
## without the introns.
.annoCalculateN <- function(annoN)
{
    annoN <- as.data.frame(annoN)
    annoN <- annoN %>% group_by(transcript_id) %>%
        mutate(full_length = sum(width)) %>%
        mutate(Rstart = min(start)) %>% mutate(Rend = max(end))
    annoN <- makeGRangesFromDataFrame(annoN,keep.extra.columns = TRUE)
    annoN <- sort(annoN, decreasing=TRUE) %>% as.data.frame()
    annoN <- annoN %>% group_by(transcript_id) %>%
        mutate(Add = c(0, cumsum(width)[-length(width)]))
    annoN <- makeGRangesFromDataFrame(annoN,keep.extra.columns = TRUE)
    return(annoN)
}

## The functions for assigning peaks
.peakAssignment <- function(granges, overlap, fullAnno)
{
    granges$location <- fullAnno$type2[overlap]
    granges$length <- fullAnno$full_length[overlap]
    granges$Rstart <- fullAnno$Rstart[overlap]
    granges$Rend <- fullAnno$Rend[overlap]
    granges$S1 <- start(fullAnno)[overlap]
    granges$E1 <- end(fullAnno)[overlap]
    granges$Add <- fullAnno$Add[overlap]
    granges$Gene_ID <- fullAnno$gene_id[overlap]
    return(granges)
}

## Calculate the position of peaks without intron
.intronOutPosition <- function(granges)
{
    granges$Position <- 5
    granges_n <- granges[strand(granges) == "-"]
    granges_p <- granges[strand(granges) == "+"]

    granges_n$Position <- (granges_n$E1-start(granges_n) +
        granges_n$Add)/granges_n$length
    granges_p$Position <- (start(granges_p) - granges_p$S1 +
        granges_p$Add) /granges_p$length
    granges <- c(granges_n,granges_p)

    granges$Position[granges$Position == -Inf] <- 5
    granges$Position[granges$Position ==  Inf] <- 5
    granges$length <- NULL
    granges$S1 <- NULL
    granges$E1 <- NULL
    granges$Add <- NULL
    granges$Rstart <- NULL
    granges$Rend <- NULL

    return(granges)
}

## Making new annotations for including the intron region
.makeTranscriptAnno <- function(anno)
{
    anno <- as.data.frame(anno)
    anno <- anno %>% group_by(transcript_id) %>%
        mutate(Rstart = min(start)) %>% mutate(Rend = max(end))
    anno$order <- seq_len(nrow(anno))
    anno <- anno %>% group_by(transcript_id) %>% filter(order == min(order))
    anno <- makeGRangesFromDataFrame(anno,keep.extra.columns = TRUE)
    GenomicRanges::start(anno) <- anno$Rstart
    GenomicRanges::end(anno) <- anno$Rend
    return(anno)
}

## Calculate position of peaks for the full transcript region (include intron)
.intronInPosition <- function(granges)
{
    granges$Position <- 5
    granges$all_length <- abs(granges$Rstart-granges$Rend)
    granges_n <- granges[strand(granges) == "-"]
    granges_p <- granges[strand(granges) == "+"]

    granges_n$Position <- (granges_n$Rend-end(granges_n))/granges_n$all_length
    granges_p$Position <- (start(granges_p) -
        granges_p$Rstart)/granges_p$all_length
    granges <- c(granges_n,granges_p)
    granges$Position[granges$Position == -Inf] <- 5
    granges$Position[granges$Position ==  Inf] <- 5
    granges$length <- NULL
    granges$all_length <- NULL
    granges$Rstart <- NULL
    granges$Rend <- NULL
    granges$S1 <- NULL
    granges$E1 <- NULL

    return(granges)
}

.plotMeta <- function(df, title, adjust)
{
    df$Position[df$location == "CDS"] <-
        df$Position[df$location == "CDS"] + 1
    df$Position[df$location == "UTR3"] <-
        df$Position[df$location == "UTR3"] + 2
    p1 <- ggplot(df, aes(x = Position)) +
        geom_density(fill = "white", alpha= 0.05,
        adjust = adjust, color = "Orange")  +
        ggtitle(title) +
        scale_x_continuous(breaks = c(0.5,1.5,2.5, 5),
        labels = c("5'UTR","CDS","3'UTR", "No Map")) +
        geom_vline(xintercept=c(1,2), linetype=2, color="cornflowerblue") +
        theme_bw()
    return(p1)
}

.plotMetaSplit <- function(df, title, adjust)
{
    df$Position[df$location == "CDS"] <-
        df$Position[df$location == "CDS"] + 1
    df$Position[df$location == "UTR3"] <-
        df$Position[df$location == "UTR3"] + 2
    p1 <- ggplot(df, aes(x = Position, color = location)) +
        geom_density(fill = "white", alpha= 0.05,
        adjust = adjust)  +
        ggtitle(title) +
        scale_x_continuous(breaks = c(0.5,1.5,2.5, 5),
        labels = c("5'UTR","CDS","3'UTR", "No Map")) +
        geom_vline(xintercept=c(1,2), linetype=2, color="cornflowerblue") +
        theme_bw() +
        geom_density(data=df, aes(x = Position), color="grey",
        linetype = 2, fill="grey", alpha=0.1)
    return(p1)
}

.plotMetaGroup <- function(df, title, adjust)
{
    df$Position[df$location == "CDS"] <-
        df$Position[df$location == "CDS"] + 1
    df$Position[df$location == "UTR3"] <-
        df$Position[df$location == "UTR3"] + 2
    p1 <- ggplot(df, aes(x=Position, color=groupIn)) +
        geom_density(fill="white",alpha= 0.05,adjust=adjust)  +
        ggtitle(title) +
        scale_x_continuous(breaks=c(0.5,1.5,2.5,5),
        labels = c("5'UTR","CDS","3'UTR","No Map"))  +
        geom_vline(xintercept = c(1,2), linetype = 2,
        color = "cornflowerblue") +
        theme_bw() + theme(legend.position="bottom") +
        guides(color=guide_legend(title="Group"))
    return(p1)
}

.plotTranscript <- function(df, title, adjust)
{
    df$Position[df$location == "CDS"] <-
        df$Position[df$location == "CDS"] + 1
    df$Position[df$location == "UTR3"] <-
        df$Position[df$location == "UTR3"] + 2
    p1 <- ggplot(df, aes(x=Position)) +
        geom_density(fill="white", alpha=0.05,
        adjust = adjust, color="Orange")  +
        ggtitle(title) +
        scale_x_continuous(breaks = c(0,1,2,3, 5),
        labels = c("5' End","Start Codon", "Stop Codon", "3' End", "No Map")) +
        geom_vline(xintercept=c(1,2), linetype=2, color="cornflowerblue") +
        theme_bw()
    return(p1)
}

.plotTranscriptSplit <- function(df, title, adjust)
{
    df$Position[df$location == "CDS"] <-
        df$Position[df$location == "CDS"] + 1
    df$Position[df$location == "UTR3"] <-
        df$Position[df$location == "UTR3"] + 2
    p1 <- ggplot(df, aes(x = Position, color = location)) +
        geom_density(fill = "white", alpha= 0.05,
        adjust = adjust)  +
        ggtitle(title) +
        scale_x_continuous(breaks = c(0,1,2,3, 5),
        labels = c("5' End","Start Codon",
        "Stop Codon", "3' End", "No Map")) +
        geom_vline(xintercept=c(1,2), linetype=2, color="cornflowerblue") +
        theme_bw() +
        geom_density(data=df, aes(x = Position), color="grey",
        linetype = 2, fill="grey", alpha=0.1)
    return(p1)
}

.plotTranscriptGroup <- function(df, title, adjust)
{
    df$Position[df$location == "CDS"] <-
        df$Position[df$location == "CDS"] + 1
    df$Position[df$location == "UTR3"] <-
        df$Position[df$location == "UTR3"] + 2
    p1 <- ggplot(df, aes(x = Position, color = groupIn)) +
        geom_density(fill = "white", alpha= 0.05, adjust = adjust) +
        ggtitle(title) +
        scale_x_continuous(breaks = c(0, 1 ,2 ,3 ,5),
        labels = c("5' End","Start Codon", "Stop Codon", "3' End", "No Map")) +
        geom_vline(xintercept = c(1,2), linetype = 2,
        color = "cornflowerblue") +
        theme_bw() + theme(legend.position="bottom") +
        guides(color=guide_legend(title="Group"))
    return(p1)
}



## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## The "metaGeneProfile" methods for GRanges objects.
##

#' @rdname metaGeneProfile
setMethod("metaGeneProfile", signature(object="GRanges"),
    function(object, annotation, include_intron=FALSE,
    title="Meta Gene Profile", group=NA, split=FALSE,
    exlevel=NA, extranscript_support_level=NA,
    adjust=1, nomap=FALSE)
    {
        if(missing(object))
            stop("The input GRanges object is missing.")
        if(!isS4(object))
            stop("The input should be a S4 GRanges object.")
        if(missing(annotation))
            stop("The pathway to the annotation file is needed.")
        if(!is.logical(include_intron))
            stop("The include_intron should be a logical vector,
            i.e. TRUE or FALSE.")
        if(!is.logical(nomap))
            stop("The nomap should be a logical vector, i.e. TRUE or FALSE.")
        if(!is.numeric(adjust))
            stop("The adjust should be a numberic, e.g. 1 or 0.5.")
        if(!is.character(title))
            stop("The title should be a character string.")
        if(length(GenomicRanges::seqnames(object)) == 0)
            stop("The input object should not be empty.")
        if (!is.na(group) &
            (group %in% colnames(elementMetadata(object)) == FALSE))
            stop("When the option group is used, please make sure that
            the group should be a column name of input GRanges object which
            contains the information of which group this peaks belongs to
            and will be used for seperate the for generating the meta gene
            profile curve for different groups")
        level <- c(1,2,3,NA)
        tsl <- c(1,2,3,4,5,6,NA)
        if (!exlevel %in% level| !extranscript_support_level %in% tsl)
            warning("The exlevel should be a vector includes the value of 1, 2,
                3 or NA. extranscript_support_level should be a vector includes
                value 1, 2, 3, 4, 5, 6 or NA.")
        if (isTRUE(split)& !is.na(group))
            stop("If you would like to use the group option please set split to
            FALSE. Conversely, when you set split to TRUE please set group
            to NA to make sure the plot is readable.")
        if(isS4(object)&length(GenomicRanges::seqnames(object)) != 0){
            ## import annotation file
            anno <- rtracklayer::import.gff3(con = annotation)

            ## Annotation filtering
            anno$transcript_support_level[
                is.na(anno$transcript_support_level)] <- 6
            anno$transcript_support_level[
                anno$transcript_support_level == "NA"] <- 6

            if (!is.na(exlevel)) {
                anno <- .annoFilterLevel(anno, exlevel)
            }
            if (!is.na(extranscript_support_level)) {
                anno <- .annoFilterTSL(anno, extranscript_support_level)
            }

            object <- .centerPeaks(object)

            CDS <- anno[anno$type == "CDS"]
            UTR3 <- anno[anno$type == "three_prime_UTR"]
            UTR5 <- anno[anno$type == "five_prime_UTR"]

            if(include_intron == FALSE){
                ## Separate different annotation regions for the
                ## calculation
                CDS_p  <- CDS[strand(CDS) == "+"] %>% sort()
                CDS_n  <- CDS[strand(CDS) == "-"] %>% sort()

                UTR3_p <- UTR3[strand(UTR3) == "+"] %>% sort()
                UTR3_n <- UTR3[strand(UTR3) == "-"] %>% sort()

                UTR5_p <- UTR5[strand(UTR5) == "+"] %>% sort()
                UTR5_n <- UTR5[strand(UTR5) == "-"] %>% sort()

                ## Get the location of annotation fragments for intronic
                ## region exclusive version.
                ## 5'UTR
                UTR5_p <- .annoCalculateP(UTR5_p)
                UTR5_n <- .annoCalculateN(UTR5_n)
                UTR5 <- sort(c(UTR5_n,UTR5_p))
                UTR5$type2 <- "UTR5"

                ## CDS
                CDS_p <- .annoCalculateP(CDS_p)
                CDS_n <- .annoCalculateN(CDS_n)
                CDS <- sort(c(CDS_n,CDS_p))
                CDS$type2 <- "CDS"

                ## 3' UTR
                UTR3_p <- .annoCalculateP(UTR3_p)
                UTR3_n <- .annoCalculateN(UTR3_n)
                UTR3 <- sort(c(UTR3_n,UTR3_p))
                UTR3$type2 <- "UTR3"
                ## Assign the transcript length for the following peak
                ## assignment
                anno.transcript <- anno[anno$type == "transcript"]
                UTR5$transcript_length <- width(
                    anno.transcript[match(UTR5$transcript_id,
                    anno.transcript$transcript_id)])
                UTR3$transcript_length <- width(
                    anno.transcript[match(UTR3$transcript_id,
                    anno.transcript$transcript_id)])
                CDS$transcript_length <- width(
                    anno.transcript[match(CDS$transcript_id,
                    anno.transcript$transcript_id)])

                ## Rank the transcript fragments base on the level,
                ## transcript support level and transcript length ##
                fullAnno <- c(UTR5,UTR3,CDS)

                fullAnno <-
                    fullAnno[order(fullAnno$level,
                    fullAnno$transcript_support_level,
                    -fullAnno$transcript_length)]

                ## Assign peaks to annotation according to
                ## level->transcript->transcript length
                ## Use select = "first" to assign peaks to the first
                ## ranked annotation
                o <- findOverlaps(object, fullAnno, select = "first")
                object <- .peakAssignment(object, o,fullAnno)

                ## Remove all the NA from the object
                object$location[is.na(object$location)] <- "NO"
                object$length[is.na(object$length)] <- 0
                object$Rstart[is.na(object$Rstart)] <- 0
                object$Rend[is.na(object$Rend)] <- 0
                object$S1[is.na(object$S1)] <- 0
                object$E1[is.na(object$E1)] <- 0
                object$Add[is.na(object$Add)] <- 0
                object$Gene_ID[is.na(object$Gene_ID)] <- "Nan"

                GenomicRanges::start(object) <- object$oriStart
                GenomicRanges::end(object) <- object$oriEnd
                object$oriStart <- NULL
                object$oriEnd <- NULL

                ## Get result for the intron excluded version
                object <- .intronOutPosition(object)
                df <- as.data.frame(object)
                if (is.na(group)) {
                    if (split == TRUE) {
                        if(nomap == FALSE){
                            df <- df[df$Position != 5, ]
                            df$location <- factor(df$location,
                                levels = c("UTR5", "CDS", "UTR3"))
                        }
                        ourpic <- .plotMetaSplit(df, title=title, adjust=adjust)
                    }
                    if (split == FALSE) {
                        ourpic <- .plotMeta(df, title=title, adjust=adjust)
                    }
                }
                if (!is.na(group)) {
                    df$groupIn <- df[,colnames(df) == group]
                    if (nomap == FALSE) {
                        ourpic <- .plotMetaGroup(df, title=title,
                            adjust=adjust)
                    }
                    if (nomap == TRUE) {
                        ourpic <- .plotMetaGroup(df, title=title, adjust=adjust)
                    }
                }
                if (nomap == FALSE) {
                    ourpic <- ourpic + coord_cartesian(xlim = c(0,3))
                }
                if (nomap == TRUE) {
                    ourpic <- ourpic + coord_cartesian(xlim = c(0,5))
                }
                ## get output
                output <- list(Peaks = object, Plot = ourpic)
            }
            if (include_intron == TRUE) {

                ## Making new annotation includes intronic region for the
                ## peaks assignment
                UTR5 <- .makeTranscriptAnno(UTR5)
                UTR5$type2 <- "UTR5"
                UTR3 <- .makeTranscriptAnno(UTR3)
                UTR3$type2 <- "UTR3"
                CDS  <- .makeTranscriptAnno(CDS)
                CDS$type2 <- "CDS"

                anno.transcript <- anno[anno$type == "transcript"]
                UTR5$transcript_length <-
                    width(anno.transcript[match(UTR5$transcript_id,
                    anno.transcript$transcript_id)])
                UTR3$transcript_length <-
                    width(anno.transcript[match(UTR3$transcript_id,
                    anno.transcript$transcript_id)])
                CDS$transcript_length <-
                    width(anno.transcript[match(CDS$transcript_id,
                    anno.transcript$transcript_id)])

                ## Making all the annotation together and rank them
                ## according to level->transcript->transcript length
                fullAnno <- c(UTR5,UTR3,CDS)

                fullAnno <-
                    fullAnno[order(fullAnno$level,
                    fullAnno$transcript_support_level,
                    -fullAnno$transcript_length)]

                ## Assign peaks
                o <- findOverlaps(object, fullAnno, select = "first")
                object <- .peakAssignment(object, o, fullAnno)

                ## Remove all the NA from the object
                object$location[is.na(object$location)] <- "NO"
                object$Rstart[is.na(object$Rstart)] <- 0
                object$Rend[is.na(object$Rend)] <- 0
                object$S1[is.na(object$S1)] <- 0
                object$E1[is.na(object$E1)] <- 0
                object$Gene_ID[is.na(object$Gene_ID)] <- "Nan"

                object <- .intronInPosition(object)

                GenomicRanges::start(object) <- object$oriStart
                GenomicRanges::end(object) <- object$oriEnd
                object$oriStart <- NULL
                object$oriEnd <- NULL

                df <- as.data.frame(object)
                if (is.na(group)) {
                    if (split == TRUE) {
                        if(nomap == FALSE){
                            df <- df[df$Position != 5, ]
                            df$location <- factor(df$location,
                                levels = c("UTR5", "CDS", "UTR3"))
                        }
                        ourpic <- .plotTranscriptSplit(df, title=title,
                            adjust=adjust)
                    }
                    if (split == FALSE) {
                        ourpic <- .plotTranscript(df, title=title,
                            adjust=adjust)
                    }
                }
                if (!is.na(group)) {
                    df$groupIn <- df[,colnames(df) == group]
                    ourpic <- .plotTranscriptGroup(df, title=title,
                            adjust=adjust)
                }
                if (nomap == FALSE) {
                    ourpic <- ourpic + coord_cartesian(xlim = c(0,3))
                }
                if (nomap == TRUE) {
                    ourpic <- ourpic + coord_cartesian(xlim = c(0,5))
                }
                output <- list(Peaks = object, Plot = ourpic)
            }
            return(output)
        }
    }
)
