#' Generate a file list for ZENBU
#' 
#' A batch of files various formats such as _BED_, _BED12_, _BAM_ etc. can
#' be uploaded to a [ZENBU](https://fantom.gsc.riken.jp/zenbu/) genome browser
#' using the `zenbu_upload` command-line tool.  This function creates a
#' tab-separated data frame and optionally writes to prepare the upload, using
#' the metatdata from the _column data_ of the `CAGEexp` object.
#' 
#' The _filelist_ table has 4 columns: 1) file path, 2) short name, 3) long
#' description, and 4) metadata.  The metadata is in semicolon-separated
#' key/value format with optional quotes like in GFF files.
#' 
#' @param x A `[CAGEexp]` object or a `DataFrame` produced by subsetting the
#'        `colData` of a `CAGEexp` object.
#' @param file The patch of the file in which to write the file list.
#' @param suffix File name suffix.  Typically `BED`, but can be longer for
#'        instance if you exported multiple kinds of BED files.
#' @param prefix A prefix to be added to the file path.
#' 
#' @family CAGEr export functions
#' 
#' @author Charles Plessy
#' 
#' @export

.ZENBU_filelist <- function(DF, file = NULL, suffix, prefix = NULL) {
  DF_sub <- DF
  DF_sub$sampleLabels <- DF_sub$inputFiles <- DF_sub$inputFilesType <- DF_sub$Description <- NULL
  if(is.null(DF$Description)) DF$Description <- ""
  out <- data.frame(path = paste(DF$sampleLabels|>unname(), suffix, sep = "."),
                    name = DF$sampleLabels|>unname(),
                    desc = DF$Description)
  out$meta <- sapply(1:nrow(DF_sub), \(n) {
    paste(colnames(DF_sub), sapply(DF_sub[n,,drop =TRUE], unname), sep = "=", collapse = ";")
  })
  if(is.null(file)) {
    return(out)
  } else {
    write.table(out, file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  }
  invisible(out)
  # TODO: add the pwd as prefix by default ?
}

