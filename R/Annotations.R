#' @include AllClasses.R CAGEexp.R

#' @name plotAnnot
#'
#' @title Plot annotation statistics
#' 
#' @description Plot maping statistics of an object containing mapping statistics in
#' counts as percentages in stacked barplots.
#' 
#' @param x An object from which can be extracted a table with columns named
#'        `promoter`, `exon`, `intron`, `mapped`, `extracted`, `rdna`, and
#'        `tagdust`, that will be passed to the `mapStats` function.
#' @param scope The name of a \dQuote{scope}, that defines which data is plotted
#'        and how it is normalised, or a function implementing that scope.
#'        See [mapStatsScopes()] for details on each scope.
#' @param title The title of the plot.
#' @param group A factor to group the samples, or the name of a `colData`
#'        column of a `CAGEexp` object, or a formula giving the names of columns
#'        to be pasted together.
#' @param facet A factor or the name of a `colData` column of a
#'        `CAGEexp` object, to facet the samples in the sense of
#'        `ggplot2`'s [facet_wrap][ggplot2::facet_wrap()] function.
#' @param normalise Whether to normalise or not. Default: TRUE.
#' 
#' @return Returns invisibly a `ggplot2` object of class `c("gg", "ggplot")`.
#' 
#' @details Stacked barplots with error bars inspired from
#' <http://stackoverflow.com/questions/10417003/stacked-barplot-with-errorbars-using-ggplot2>.
#' See <http://www.biomedcentral.com/1471-2164/14/665/figure/F1> for example.
#' 
#' @seealso [mapStats()] for a list of _scopes_.
#' @family CAGEr annotation functions
#' @family CAGEr plot functions
#' 
#' @author Charles Plessy
#' 
#' @examples
#' p <- plotAnnot(exampleCAGEexp, 'counts', 'Here is the title')
#' print(p)
#' p + ggplot2::theme_bw()
#' ggplot2::theme_set(ggplot2::theme_bw()) ; p
#' plotAnnot(exampleCAGEexp, 'counts', 'Same, non-normalised', normalise = FALSE)
#' exampleCAGEexp$myGroups <- factor(c("A", "A", "B", "B", "C"))
#' plotAnnot(exampleCAGEexp, 'counts', group = "myGroups")
#' plotAnnot(exampleCAGEexp, 'counts', group = ~myGroups)
#' plotAnnot(exampleCAGEexp, 'counts', group = ~sampleLabels + myGroups)
#' plotAnnot(exampleCAGEexp, CAGEr:::msScope_counts , group = "myGroups")
#' 
#' @docType methods
#' @importFrom ggplot2 aes_string coord_flip geom_bar geom_segment geom_point
#' @importFrom ggplot2 ggplot ggtitle position_stack facet_wrap
#' @export

setGeneric("plotAnnot", function( x, scope, title
                                , group = "default", facet = NULL
                                , normalise = TRUE)
  standardGeneric("plotAnnot"))

#' @rdname plotAnnot

setMethod("plotAnnot", "data.frame",
  function( x, scope, title, group, facet, normalise) {
  p <- ggplot( mapStats( x, scope = scope
                       , group = group, facet = facet, normalise = normalise)
        , aes_string( x    = "group"
                    , y    = "value"
                    , fill = "variable")
        , main = all) +
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    geom_segment(aes_string( xend = "group"
                           , y    = "ystart"
                           , yend = "yend")) +
    geom_point( aes_string( x = "group"
                          , y = "yend")
              , shape = "|") +
    coord_flip() +
    ggtitle(title)
  if(! is.null(facet)) {
    # mapStats always returns the facet in a column named "facet"
    p + facet_wrap("facet")
  } else {
    p
  }
})

#' @rdname plotAnnot

setMethod("plotAnnot", "DataFrame",
  function( x, scope, title, group, facet, normalise) {
  plotAnnot( data.frame(x, check.names = FALSE)
           , scope = scope, title = title
           , group = group, facet = facet, normalise = normalise)
})

#' @rdname plotAnnot
#' @importFrom formula.tools lhs rhs.vars

setMethod("plotAnnot", "CAGEexp",
  function( x, scope, title, group, facet, normalise) {
  if (missing(group)) {
    group <- sampleLabels(x)
  } else if (length(group) == 1) {
    if (! exists(group, data.frame(colData(x))))
      stop("Could not find group ", dQuote(group), ".")
    group <- colData(x)[[group]]
  } else if (is(group, "formula")) {
    if(!is.null(lhs(group)))
      stop("Formula must start with a tilde.")
    vars <- rhs.vars(group)
    for (var in vars)
      if(! var %in% colnames(colData(x)))
        stop("Column ", dQuote(var), " not found in sample metadata table.")
    group <- do.call(paste, colData(x)[vars])
  } else {
    stop("Could not find group ", dQuote(group), ".")
  }
  if (missing(title)) title <- paste("CAGEr object", dQuote(deparse(substitute(x))))
  plotAnnot( colData(x)
           , scope = scope, title = title
           , group = group, facet = facet, normalise = normalise)
})

#' @name mapStats
#' 
#' @title Process mapping statistics
#' 
#' @description Using a data frame containing mapping statistics in counts,
#' transform the data in percentages that can be used for stacked barplots.
#' 
#' @param libs A data frame with containing columns required by the `scope`
#'        chosen.
#' @param scope The name of a \dQuote{scope}, that defines which data is plotted
#'        and how it is normalised, or a function that implements a custom scope.
#'        See [mapStatsScopes()] for details on each scope.
#' @param group A vector of factors defining groups in the data.  By default,
#'        the \dQuote{group} column of the \dQuote{libs} table.
#' @param facet A vector of factors defining facets in the data (in the sense
#'        of `ggplot2`'s [facet_wrap][ggplot2::facet_wrap()] function).
#' @param normalise Whether to normalise or not. Default: `TRUE`.
#'
#' @return Returns a data frame with mean and standard deviation of normalised
#' mapping statistics, plus absolute positions for the error bars.  The first
#' column, `group`, is a vector of factors sorted with the [gtools::mixedorder()]
#' function.  The facet column, if any, is always called `facet`.
#' 
#' @details See the plotAnnot vignette and the [mapStatsScopes()]
#' help page for details on what the scopes are.
#' 
#' See <http://stackoverflow.com/questions/10417003/stacked-barplot-with-errorbars-using-ggplot2> about stacked barplot.
#' 
#' @author Charles Plessy
#' 
#' @seealso [plotAnnot], [mapStatsScopes]
#' 
#' @examples
#' library(SummarizedExperiment)
#' CAGEr:::mapStats(as.data.frame(colData(exampleCAGEexp)), "counts", sampleLabels(exampleCAGEexp))
#' CAGEr:::mapStats(as.data.frame(colData(exampleCAGEexp)), "counts", c("A", "A", "B", "B", "C"))
#' 
#' @importFrom gtools mixedorder
#' @importFrom plyr ddply
#' @importFrom reshape melt

mapStats <- function( libs
                    , scope
                    , group="default"
                    , facet = NULL
                    , normalise = TRUE) {
  
  if (identical(group, "default")) {
      if        ("group" %in% colnames(libs)) {
        group <- libs$group
      } else if ("Group" %in% colnames(libs)) {
        group <- libs$Group
    } else
      stop(paste("Missing", dQuote("group"), "column in the data frame."))
  }
  
  # Backup levels for later.  Coerce to factor if it was not.  This way,
  # numerical ordering is preserved despite the conversion to characters
  # when facetting
  group.levels <- levels(factor(group))
  if (!is.null(facet))
    facet.levels <- levels(factor(libs[,facet]))

  if (! ("tagdust" %in% colnames(libs))) libs[, "tagdust"] <- 0
  
  if (! is.function(scope))
    scope <- switch( scope
                   , all        = msScope_all
                   , annotation = msScope_annotation
                   , counts     = msScope_counts
                   , mapped     = msScope_mapped
                   , qc         = msScope_qc
                   , steps      = msScope_steps
                   , function(libs) stop("Unknown scope", call. = FALSE))
  custom.list <- scope(libs)
  libs    <- custom.list$libs
  columns <- custom.list$columns
  total   <- custom.list$total
  
  if (normalise == FALSE) total <- 1

  doMean <- function (X) tapply(libs[,X] / total, group, mean)
  doSd   <- function (X) tapply(libs[,X] / total, group, sd  )
  
  if (! is.null(facet)) {
    if(! facet %in% colnames(libs)) stop("Missing ", dQuote(facet), " column.")
    group <- paste(group, libs[,facet], sep = "__FACET__")
  }
  
  # "simplify" needs to be FALSE so that conversion to data frame works even
  # when the group contains only a single level.
  mapstats          <- data.frame(sapply(columns, doMean, simplify = FALSE))
  mapstats$group    <- rownames(mapstats)
  mapstats[gtools::mixedorder(mapstats$group), ]
  mapstats$group    <- factor(mapstats$group, unique(mapstats$group))
  
  mapstats.sd       <- data.frame(sapply(columns, doSd, simplify = FALSE))
  mapstats.sd$group <- rownames(mapstats.sd)
  
  mapstats          <- reshape::melt(mapstats,    id.vars="group")
  mapstats$sd       <- reshape::melt(mapstats.sd, id.vars="group")$value
  
  value <- NULL # To silence "no visible binding for global variable" error in R CMD check.
  mapstats          <- plyr::ddply( mapstats
                                  , plyr::.(group)
                                  , transform
                                  , ystart = cumsum(value)
                                  , yend   = cumsum(value) + sd)
  if (! is.null(facet)) {
    mapstats$facet <- sub(".*__FACET__", "", mapstats$group)
    mapstats$facet <- factor(mapstats$facet, levels = facet.levels)
    
    mapstats$group <- sub("__FACET__.*", "", mapstats$group)
    mapstats$group <- factor(mapstats$group, levels = group.levels)
  }
  mapstats
}

.checkLibsDataFrame <- function(libs, columns) {
  if (! all(columns %in% colnames(libs)))
    stop( "Input data frame needs the following columns:\n"
        , paste(columns, collapse = " "))
}

#' @name mapStatsScopes
#' @aliases mapStatsScopes
#' 
#' @title mapStats scopes
#' 
#' @description Functions implementing the `scope` parameter of the
#' `\link{mapStats}` function.
#' 
#' @param libs A data frame containing metadata describing samples in sequence
#'        libraries.
#' 
#' @return Returns a list with three elements: `libs` contains a modified
#' version of the input data frame where columns have been reorganised as needed,
#' `colums` contains the names of the columns to use for plotting and
#' provides the order of the stacked bars of the `plotAnnot` function,
#' `total` indicates the total counts used for normalising the data.

#' @rdname mapStatsScopes
#' @details The **`counts`** scope reports the number of molecules aligning in
#' _promoter_, _exon_, _intron_ and otherwise _intergenic_.
#' regions.

msScope_counts <- function(libs) {
  .checkLibsDataFrame(libs, c("promoter", "exon", "intron", "librarySizes"))
  libs$Promoter   <- libs$promoter
  libs$Exon       <- libs$exon
  libs$Intron     <- libs$intron
  libs$Intergenic <- libs$librarySizes - libs$promoter - libs$intron - libs$exon
  list( libs    = libs
      , columns = c("Promoter", "Exon", "Intron", "Intergenic")
      , total   = libs$librarySizes)
}

#' @rdname mapStatsScopes
#' @details The **`mapped`** scope reports the number of molecules aligning in
#' _promoter_, _exon_, _intron_ and otherwise _intergenic_,
#' plus the number of PCR duplicates (mapped tags minus molecule counts), plus
#' the number of non-properly paired mapped tags.

msScope_mapped <- function(libs) {
  .checkLibsDataFrame(libs, c( "promoter", "exon", "intron"
                             , "mapped", "properpairs", "librarySizes"))
  libs$Non_proper <- libs$mapped       - libs$properpairs
  libs$Duplicates <- libs$properpairs  - libs$librarySizes
  libs$Intergenic <- libs$librarySizes - libs$promoter - libs$intron - libs$exon
  libs$Intron     <- libs$intron
  libs$Exon       <- libs$exon
  libs$Promoter   <- libs$promoter
  list( libs    = libs
      , columns = c( "Promoter", "Exon", "Intron", "Intergenic"
                   , "Duplicates", "Non_proper")
      , total   = libs$mapped)
}

#' @rdname mapStatsScopes
#' @details The **`qc`** scope reports the number of tags removed as 
#' _tag dust_, _rRNA_, _spikes_, plus the _unmapped_ tags,
#' plus the number of non-properly paired mapped tags, plus the number of PCR
#' duplicates (mapped tags minus molecule counts), plus the number of unique
#' molecule counts.

msScope_qc <- function(libs) {
  .checkLibsDataFrame(libs, c("extracted", "rdna", "spikes", "cleaned", "mapped", "properpairs", "librarySizes"))
  libs$Tag_dust     <- libs$extracted   - libs$rdna - libs$spikes - libs$cleaned
  libs$rDNA         <- libs$rdna
  libs$Spikes       <- libs$spikes
  libs$Unmapped     <- libs$cleaned     - libs$mapped
  libs$Non_proper   <- libs$mapped      - libs$properpairs
  libs$Duplicates   <- libs$properpairs - libs$librarySizes
  libs$Counts       <- libs$librarySizes
  list( libs    = libs
      , columns = c( "Tag_dust", "rDNA", "Spikes", "Unmapped"
                   , "Non_proper", "Duplicates", "Counts")
      , total   = libs$extracted)
}

#' @rdname mapStatsScopes
#' @details The **`steps`** scope reports the number of tags removed by
#' _cleaning_, _mapping_, and _deduplication_, plus the number
#' of _unique molecule counts_.

msScope_steps <- function(libs) {
  .checkLibsDataFrame(libs, c( "extracted", "cleaned", "properpairs"
                             , "librarySizes", "extracted"))
  libs$Cleaning      <- libs$extracted   - libs$cleaned
  libs$Mapping       <- libs$cleaned     - libs$properpairs
  libs$Deduplication <- libs$properpairs - libs$librarySizes
  libs$Counts        <- libs$librarySizes
  total   <- libs$extracted
  columns <- c("Cleaning", "Mapping", "Deduplication", "Counts")
  if ("total" %in% colnames(libs)) {
    total <- libs$total
    libs$Extraction <- with(libs, total - extracted)
    columns <- c("Extraction", columns)
  }
  list( libs    = libs
      , columns = columns
      , total   = total)
}

#' @rdname mapStatsScopes
#' @details The legacy **`all`** scope reports the number of tags in
#' _promoters_, _exons_, _introns_, or _mapped_ elswhere, or removed because
#' they match rRNA or are likely primer artefacts, normalised by the total
#' nubmer of extracted tags.

msScope_all <- function(libs) {
  .checkLibsDataFrame(libs, c( "mapped", "promoter", "intron", "exon", "rdna"
                             , "tagdust", "extracted"))
  libs$mapped <- libs$mapped - libs$promoter - libs$intron - libs$exon
  list( libs    = libs
      , columns = c("promoter", "exon", "intron", "mapped", "rdna", "tagdust")
      , total   = libs$extracted)
}

#' @rdname mapStatsScopes
#' @details The legacy **`annotation`** scope reports the number of tags in
#' _promoters_, _exons_, _introns_, or _mapped_ elswhere, or removed because
#' they match rRNA or are likely primer artefacts, normalised by the total
#' nubmer of mapped tags.

msScope_annotation <- function(libs) {
  .checkLibsDataFrame(libs, c( "mapped", "promoter", "intron", "exon", "rdna"
                             , "tagdust", "extracted"))
  libs$mapped <- libs$mapped - libs$promoter - libs$intron - libs$exon
  list( libs    = libs
      , columns = c("promoter", "exon", "intron", "mapped", "rdna", "tagdust")
      , total   = libs$mapped)
}


#' @name annotateCTSS
#' 
#' @title Annotate and compute summary statistics
#' 
#' @description `annotateCTSS` annotates the _CTSS_ of a [`CAGEexp`] object and
#' computes annotation statistics.
#' 
#' @param object `CAGEexp` object (`CAGEset`s are not supported).
#'   
#' @param ranges A [`GRanges`] object, optionally containing `gene_name`,
#'   `type` and `transcript_type` metadata.
#'   
#' @return `annotateCTSS` returns the input object with the following
#'   modifications:
#' 
#'  * The Genomic Ranges of the `tagCountMatrix` experiment gains an
#'    `annotation` metadata column, with levels such as `promoter`,
#'     `exon`, `intron` and `unknown`.  If the annotation has a `gene_name`
#'     metadata, then a `genes` column is also added, with gene symbols from
#'     the annotation.
#'  * The sample metadata gets new columns, indicating total counts in each of
#'    the annotation levels.  If the annotation has a `gene_name` metadata, then
#'     a `genes` column is added to indicate the number of different gene symbols
#'     detected.
#' 
#' @seealso [`CTSStoGenes`], and the [`exampleZv9_annot`] example data.
#' @family CAGEr object modifiers
#' @family CAGEr annotation functions
#' 
#' @author Charles Plessy
#' 
#' @examples 
#' library(SummarizedExperiment)
#' annotateCTSS(exampleCAGEexp, exampleZv9_annot)
#' colData(exampleCAGEexp)
#' 
#' @export

setGeneric("annotateCTSS", function(object, ranges) standardGeneric("annotateCTSS"))

#' @rdname annotateCTSS

setMethod("annotateCTSS", "CAGEset", function (object, ranges){
  stop("CAGEset objects not supported.")})

#' @rdname annotateCTSS

setMethod("annotateCTSS", c("CAGEexp", "GRanges"), function (object, ranges){
  objName <- deparse(substitute(object))
  
  CTSScoordinatesGR(object)$genes      <- ranges2genes(CTSScoordinatesGR(object), ranges)
  CTSScoordinatesGR(object)$annotation <- ranges2annot(CTSScoordinatesGR(object), ranges)
  
  annot <- sapply( CTSStagCountDF(object)
                 , function(X) tapply(X, CTSScoordinatesGR(object)$annotation, sum))
  colData(object)[levels(CTSScoordinatesGR(object)$annotation)] <- DataFrame(t(annot))
  
  if (validObject(object)) {
    assign(objName, object, envir = parent.frame())
    invisible(1)
  }
})

#' @name annotateConsensusClusters
#' 
#' @rdname annotateCTSS
#' 
#' @description `annotateConsensusClusters` annotates the _consensus clusters_
#' of a CAGEr object.
#'   
#' @return `annotateConsensusClusters` returns the input object with the same
#' modifications as above.
#' 
#' @examples 
#' annotateConsensusClusters(exampleCAGEexp, exampleZv9_annot)
#' consensusClustersGR(exampleCAGEexp)
#' 
#' @export

setGeneric("annotateConsensusClusters", function(object, ranges) standardGeneric("annotateConsensusClusters"))

#' @rdname annotateCTSS

setMethod("annotateConsensusClusters", "CAGEset", function (object, ranges){
  stop("CAGEset objects not supported.")})

#' @rdname annotateCTSS

setMethod("annotateConsensusClusters", c("CAGEexp", "GRanges"), function (object, ranges){
  objName <- deparse(substitute(object))
  if(is.null(experiments(object)$tagCountMatrix))
    stop(objName, " does not contain CTSS expressiond data, see ", dQuote("getCTSS()"), ".")
  consensusClustersGR(object)$annotation <- ranges2annot(consensusClustersGR(object), ranges)
  if(!is.null(ranges$gene_name))
    consensusClustersGR(object)$genes    <- ranges2genes(consensusClustersGR(object), ranges)
  if (validObject(object)) {
    assign(objName, object, envir = parent.frame())
    invisible(1)
  }
})


#' @name ranges2annot
#' 
#' @title Hierarchical annotation of CTSSes
#' 
#' @description Assigns region types such as `promoter`, `exon` or `unknown`
#'              to CTSSes.
#' 
#' @param ranges A [`CTSS`] object, for example extracted from a
#'        `RangedSummarizedExperiment` object with the [`rowRanges`]
#'        command.
#' 
#' @param annot A [`GRanges`] from which promoter positions will be inferred.
#'        Typically GENCODE.  If the `type` metadata is present, it should
#'        contain `gene`, `exon` and `transcript` among its values.  Otherwise,
#'        all entries are considered transcripts. If the `transcript_type`
#'        metadata is available, the entries that may not be primary products
#'        (for instance \sQuote{snoRNA}) are discarded.
#' 
#' @return A Run-length-encoded ([`Rle`]) factor of same length as the `CTSS`
#'         object, indicating if the interval is `promoter`, `exon`, `intron` or
#'         `unknown`, or just `promoter`, `gene`, `unknown` if the `type`
#'         metadata is absent.
#' 
#' @details Only the biotypes that are likely to have a pol II promoter will be
#' filtered in.  This is currently hardcoded in the function; see its source
#' code.  Example of biotypes without a pol II promoter: VDJ segments, miRNA,
#' but also snoRNA, etc.  Thus, the _Intergenic_ category displayed in output of
#' the [`plotAnnot`] may include counts overlaping with real exons of discarded
#' transcribed regions: be careful that large percentages do not necessarly
#' suggest abundance of novel promoters.
#' 
#'         
#' @family CAGEr annotation functions
#' @seealso [`CTSScoordinatesGR`], [`exampleZv9_annot`]
#' 
#' @author Charles Plessy
#' 
#' @examples
#' CAGEr:::ranges2annot(CTSScoordinatesGR(exampleCAGEexp), exampleZv9_annot)
#' 
#' ctss <- GenomicRanges::GRanges("chr1", IRanges::IPos(c(1,100,200,1500)), "+")
#' ctss <- GenomicRanges::GPos(ctss, stitch = FALSE)
#' ctss <- CAGEr:::.CTSS(ctss)
#' gr1   <- GenomicRanges::GRanges( "chr1"
#'                                , IRanges::IRanges(c(650, 650, 1400), 2000), "+")
#' CAGEr:::ranges2annot(ctss, gr1)
#' gr2 <- gr1
#' gr2$type            <- c("transcript",     "exon",           "transcript")
#' gr2$transcript_type <- c("protein_coding", "protein_coding", "miRNA")
#' CAGEr:::ranges2annot(ctss, gr2)
#' 
#' @importFrom GenomicRanges findOverlaps promoters
#' @importFrom S4Vectors Rle

ranges2annot <- function(ranges, annot) {
  typesWithPromoter <- c( "protein_coding", "processed_transcript", "lincRNA"
                        , "antisense", "processed_pseudogene"
                        , "unprocessed_pseudogene")
  if(!is.null(annot$transcript_type))
    annot <- annot[annot$transcript_type %in% typesWithPromoter]
  
  findOverlapsBool <- function(A, B) {
    overlap <- findOverlaps(A, B)
    overlap <- as(overlap, "List")
    any(overlap)
  }
  
  if(!is.null(annot$type)) {
    classes <- c("promoter", "exon", "intron", "unknown")
    p <- findOverlapsBool(ranges, promoters(annot[annot$type == "transcript"], 500, 500))
    e <- findOverlapsBool(ranges, annot[annot$type == "exon"])
    t <- findOverlapsBool(ranges, annot[annot$type == "transcript"])
    annot <- sapply( 1:length(ranges), function(i) {
      if      (p[i]) {classes[1]}
      else if (e[i]) {classes[2]}
      else if (t[i]) {classes[3]}
      else           {classes[4]}
    })
  } else {
    classes <- c("promoter", "gene", "unknown")
    p <- findOverlapsBool(ranges, promoters(annot, 500, 500))
    g <- findOverlapsBool(ranges, annot)
    annot <- sapply( 1:length(ranges), function(i) {
      if      (p[i]) {classes[1]}
      else if (g[i]) {classes[2]}
      else           {classes[3]}
    })
  }

  annot <- factor(annot, levels = classes)
  Rle(annot)
}

#' ranges2genes
#' 
#' Assign gene symbol(s) to Genomic Ranges.
#' 
#' This private (non-exported) function is used to assign gene symbols
#' to genomic ranges.  It is run by \code{\link{annotateCTSS}}, which has to
#' be run before \code{\link{CTSStoGenes}}.
#' 
#' @param ranges Genomics Ranges object, for example extracted from a
#'               RangedSummarizedExperiment object with the \code{rowRanges}
#'               command.
#' 
#' @param genes A \code{\link{GRanges}} object containing \code{gene_name} metadata.
#' 
#' @return A \code{\link{Rle}} character vector of same length as the GRanges object,
#' indicating one gene symbol or a semicolon-separated list of gene symbols for each
#' range.
#'         
#' @family CAGEr annotation functions
#' @family CAGEr gene expression analysis functions
#' @seealso \code{\link{CTSScoordinatesGR}}, \code{\link{exampleZv9_annot}}
#' 
#' @author Charles Plessy
#' 
#' @examples
#' CAGEr:::ranges2genes(CTSScoordinatesGR(exampleCAGEexp), exampleZv9_annot)
#' 
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors List Rle unstrsplit
#' @importFrom IRanges extractList

ranges2genes <- function(ranges, genes) {
  if (is.null(genes$gene_name))
    stop("Annotation must contain ", dQuote("gene_name"), " metdata.")
  names(genes) <- genes$gene_name
  ranges2names(ranges, genes)
}


#' ranges2names
#' 
#' Intersection of genomic ranges
#' 
#' This private (non-exported) function intersects two genomic ranges and
#' for each element of the first object returns the name of the elements of
#' the second object that it intersects with.
#' 
#' @param rangesA A \code{\link{GRanges}} object.
#' @param rangesB A second GRanges object.
#' 
#' @return A \code{\link{Rle}} character vector of same length as the \code{rangesA}
#' GRanges object, indicating one name or a semicolon-separated list of names from
#' the each \code{rangesB} object.
#'         
#' @family CAGEr annotation functions
#' 
#' @author Charles Plessy
#' 
#' @examples
#' names(exampleZv9_annot) <- exampleZv9_annot$gene_name
#' CAGEr:::ranges2names(CTSScoordinatesGR(exampleCAGEexp), exampleZv9_annot)
#' 
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors List Rle unstrsplit
#' @importFrom IRanges extractList

ranges2names <- function(rangesA, rangesB) {
  if (is.null(names(rangesB)))
    stop(sQuote("rangesB"), " must contain have names.")
  names <- findOverlaps(rangesA, rangesB)
  names <- as(names, "List")
  names <- extractList(names(rangesB), names)
  names <- unique(names)
  names <- unstrsplit(names, ";")
  Rle(names)
}


#' Example zebrafish annotation data
#' 
#' Annotation data for zebrafish's chromosome 17's interval  26000000-54000000
#' (Zv9/danRer7 genome), to be used in documentation examples.
#'
#' @author Prepared by Charles Plessy \email{plessy@riken.jp} using archive ENSEMBL data.
#' @references \url{http://mar2015.archive.ensembl.org/biomart/}
#' 
#' @details Data was retreived from ENSEMBL's Biomart server using a query to extract
#' gene, transcripts and exon coordinates.  For the record, here it is as URL
#' (long, possibly overflowing).
#' 
#' http://mar2015.archive.ensembl.org/biomart/martview/78d86c1d6b4ef51568ba6d46f7d8b254?VIRTUALSCHEMANAME=default&ATTRIBUTES=drerio_gene_ensembl.default.structure.ensembl_gene_id|drerio_gene_ensembl.default.structure.ensembl_transcript_id|drerio_gene_ensembl.default.structure.start_position|drerio_gene_ensembl.default.structure.end_position|drerio_gene_ensembl.default.structure.transcript_start|drerio_gene_ensembl.default.structure.transcript_end|drerio_gene_ensembl.default.structure.strand|drerio_gene_ensembl.default.structure.chromosome_name|drerio_gene_ensembl.default.structure.external_gene_name|drerio_gene_ensembl.default.structure.gene_biotype|drerio_gene_ensembl.default.structure.exon_chrom_start|drerio_gene_ensembl.default.structure.exon_chrom_end|drerio_gene_ensembl.default.structure.is_constitutive|drerio_gene_ensembl.default.structure.rank&FILTERS=&VISIBLEPANEL=resultspanel
#' 
#' And here it is as XML.
#' 
#' \preformatted{<?xml version="1.0" encoding="UTF-8"?>
#' <!DOCTYPE Query>
#' <Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
#'   <Dataset name = "drerio_gene_ensembl" interface = "default" >
#'     <Attribute name = "ensembl_gene_id" />
#'     <Attribute name = "ensembl_transcript_id" />
#'     <Attribute name = "start_position" />
#'     <Attribute name = "end_position" />
#'     <Attribute name = "transcript_start" />
#'     <Attribute name = "transcript_end" />
#'     <Attribute name = "strand" />
#'     <Attribute name = "chromosome_name" />
#'     <Attribute name = "external_gene_name" />
#'     <Attribute name = "gene_biotype" />
#'     <Attribute name = "exon_chrom_start" />
#'     <Attribute name = "exon_chrom_end" />
#'     <Attribute name = "is_constitutive" />
#'     <Attribute name = "rank" />
#'   </Dataset>
#' </Query>}
#' 
#' The downloaded file was then transformed as follows.
#' 
#' \preformatted{x <- read.delim("~/Downloads/mart_export.txt", stringsAsFactors = FALSE)
#' e <- GRanges(paste0("chr", x$Chromosome.Name), IRanges(x$Exon.Chr.Start..bp., x$Exon.Chr.End..bp.), ifelse(x$Strand + 1, "+", "-"))
#' e$gene_name <- Rle(x$Associated.Gene.Name)
#' e$transcript_type <- Rle(x$Gene.type)
#' e$type <- "exon"
#' e$type <- Rle(e$type)
#'
#' e <- GRanges(paste0("chr", x$Chromosome.Name), IRanges(x$Exon.Chr.Start..bp., x$Exon.Chr.End..bp.), ifelse(x$Strand + 1, "+", "-"))
#' e$gene_name <- Rle(x$Associated.Gene.Name)
#' e$transcript_type <- Rle(x$Gene.type)
#' e$type <- "exon"
#' e$type <- Rle(e$type)
#' e <- sort(unique(e))
#' 
#' g <- GRanges( paste0("chr", x$Chromosome.Name)
#'             , IRanges(x$Gene.Start..bp., x$Gene.End..bp.)
#'             , ifelse( x$Strand + 1, "+", "-"))
#'             
#' g$gene_name <- Rle(x$Associated.Gene.Name)
#' g$transcript_type <- Rle(x$Gene.type)
#' g$type <- "gene"
#' g$type <- Rle(g$type)
#' g <- sort(unique(g))
#' 
#' t <- GRanges( paste0("chr", x$Chromosome.Name)
#'             , IRanges(x$Transcript.Start..bp., x$Transcript.End..bp.)
#'             , ifelse( x$Strand + 1, "+", "-"))
#'             
#' t$gene_name <- Rle(x$Associated.Gene.Name)
#' t$transcript_type <- Rle(x$Gene.type)
#' t$type <- "transcript"
#' t$type <- Rle(t$type)
#' t <- sort(unique(t))
#' 
#' gff <- sort(c(g, t, e))
#' gff <- gff[seqnames(gff) == "chr17"]
#' gff <- gff[start(gff) > 26000000 & end(gff) < 54000000]
#' seqlevels(gff) <- seqlevelsInUse(gff)
#' 
#' save(gff, "data/exampleZv9_annot.RData", compress = "xz")}

"exampleZv9_annot"


#' Make a gene expression table.
#' 
#' Add a gene expression table in the `GeneExpSE` experiment slot of an
#' annotated [`CAGEexp`] object.  [`CAGEset`] objects are not supported.
#' 
#' @param object A `CAGEexp` object that was annotated with the [annotateCTSS()]
#'        function.
#' 
#' @return The input object with the following modifications:
#' 
#'  * A new `geneExpMatrix` experiment containing gene expression levels as
#'    a [`SummarizedExperiment`] object with one assay called `counts`, which
#'    is plain `matrix` of integers.  (This plays better than `Rle DataFrames`
#'    when interfacing with downstream packages like DESeq2, and since the number of
#'    genes is limited, a `matrix` will not cause problems of performance.)
#'  * New `genes` column data added, indicating total number of gene symbols
#'    detected per library.
#'  * New `unannotated` column data added, indicating for each sample the
#'     number of counts that did not overlap with a known gene.
#' 
#' @author Charles Plessy
#' 
#' @seealso [annotateCTSS()].
#' 
#' @family CAGEr object modifiers
#' @family CAGEr gene expression analysis functions
#' 
#' @examples 
#' CTSStoGenes(exampleCAGEexp)
#' all( librarySizes(exampleCAGEexp) -
#'      colSums(SummarizedExperiment::assay(GeneExpSE(exampleCAGEexp))) ==
#'      exampleCAGEexp$unannotated)
#' 
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export

setGeneric("CTSStoGenes", function(object) standardGeneric("CTSStoGenes"))

#' @rdname CTSStoGenes

setMethod("CTSStoGenes", "CAGEset", function (object)
  stop("Not supported for ", dQuote("CAGEset"), " objects."))

#' @rdname CTSStoGenes

setMethod("CTSStoGenes", "CAGEexp", function (object) {
  objName <- deparse(substitute(object))
  if (is.null(CTSScoordinatesGR(object)$genes))
    stop(objName, " is not annotated, see ", dQuote("annotateCTSS()"), ".")
  genes <- rowsum(CTSStagCountDf(object), as.factor(CTSScoordinatesGR(object)$genes))
  object$unannotated <- unname(unlist(genes[1,]))
  genes <- genes[-1,]
  GeneExpSE(object) <- SummarizedExperiment( assays  = SimpleList(counts = as.matrix(genes))
                                           , rowData = DataFrame(symbol = rownames(genes)))
  object$genes      <- colSums(assay(GeneExpSE(object)) > 0)
  # object$geneSymbols <- countSymbols(assay(GeneExpSE(object)) %>% as.data.frame)
  if (validObject(object)) {
    assign(objName, object, envir = parent.frame())
    invisible(1)
  }
})
