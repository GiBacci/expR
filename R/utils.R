#' Format output function
#'
#' A prefix will be added to the name fo the
#' file provided and its directory will be changed
#' as specified in the [outdir] param
#'
#' @param prefix character. The prefix to add
#' @param outdir character. The output directory
#' @param suffix character. The suffix that will be
#'  added at the end of the file name
#'
#' @return a file path
#' @export
#'
#' @examples
formatOutput <- function(prefix, outdir = NULL, suffix = ""){
  function(files){
    if(is.null(outdir)){
      outdir <- dirname(files)
    }
    filename <- basename(files)
    output <- file.path(outdir, paste0(prefix, filename, suffix))
    output
  }
}

#' Format function for sample definition
#'
#' Sample names will be converted from file
#' names using this function. File names
#' will be splitted depending on the \code{sep}
#' param and the element corresponding to the \code{index}
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

#' Get run id from sequence header
#'
#' This function is used to detect the run id of each file.
#' The id is parsed from the sequence header.
#'
#' @param path Character vector, the paths to the sequence file/s
#' @param nrecords Numeric, the number of records that must be read
#'  before returning the run id
#'
#' @export
#'
#' @return the run id as character
runFromHeader <- function(nrecords = 1){
  function(files){
    sapply(files, function(p){
      con <- file(p)
      if(summary(con)$class == "gzfile"){
        close(con)
        con <- gzfile(p)
      }
      on.exit(close(con))

      res <- lapply(1:nrecords, readLines,
                    con = con, n = 4)

      ids <- sapply(res, "[", 1)
      run <- unique(sapply(strsplit(ids, ":"),
                           "[", 2))

      if(length(run) > 1){
        msg <- sprintf("More than one run found in file: %s",
                       basename(path))
        stop(msg)
      }

      return(run)
    }, USE.NAMES = F)
  }
}

#' Get run id from directory
#'
#' This function returns the run id
#' by parsign the name of the directory
#' in which the file/s is/are
#'
#' @return the run id
#' @export
#'
#' @examples
runFromBasedir <- function(){
  function(x){
    basename(dirname(x))
  }
}
