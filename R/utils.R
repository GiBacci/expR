#' Format output function
#'
#' A prefix will be added to the name fo the
#' file provided and its directory will be changed
#' as specified in the [outdir] param
#'
#' @param prefix character. The prefix to add
#' @param outdir character. The output directory
#'
#' @return a file path
#' @export
#'
#' @examples
formatOutput <- function(prefix, outdir = NULL){
  function(files){
    if(is.null(outdir)){
      outdir <- dirname(files)
    }
    filename <- basename(files)
    output <- file.path(outdir, paste0(prefix, filename))
    output
  }
}

#' Format function for sample definition
#'
#' Sample names will be converted from file
#' names using this function. File names
#' will be splitted depending on the [sep]
#' param and the element corresponding to the [index]
#' will be returned
#'
#' @param index numeric. The element to return
#' @param sep ragular expression. Input files will be splitted
#'  accordingly
#'
#' @return formatted sample names
#' @export
#'
#' @examples
formatSample <- function(index = 1, sep = "[[:punct:]]"){
  function(files){
    f <- basename(files)
    splitted <- strsplit(f, split = sep)
    sapply(splitted, "[", index)
  }
}
