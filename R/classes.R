#' @import methods
NULL

#' Main targeted experiment class
#'
#' This class represents a sequencing experiment usually
#' involving sequencing a PRC products (such as 16S rRNA
#' gene). Each sample may have one or two files (single-end
#' or paired-end sequencing) and a run identifier.
#'
#' @slot samples Character vector, sample names
#' @slot run Character vector, the run id for each sample
#' @slot step Character, the step of analysis
#' @slot n Numeric vector, the number of samples
#'
#' @return a Targeted experiment
#'
#' @example
.TargetedExperimet <- setClass("TargetedExperimet",
                               slots = c(
                                 samples = "character",
                                 run = "character",
                                 step = "character",
                                 n = "numeric"
                               ))

#' Single end experiment class
#'
#' This class represent a standard single-end experiment.
#' In addition to the [TargetedExperimet][.TargetedExperimet]
#' class it provides a vector of file paths and a list of
#' output paths. They can both be used in the context
#' of other methods such as [by_sample()], for
#' sample wise computation and [by_run()] for run wise
#' computation.
#'
#' Instances of this class should be created using the
#' [getSingleEndFromFilename()] function.
#'
#' @slot files the paths to sequence files
#' @slot output the paths to output file
#' @seealso [getSingleEndFromFilename()], for class
#'  creation and [by_sample()] and [by_run()] for
#'  computation
#'
#' @return a single end experiment
#'
#' @example
.SingleEndSamples <- setClass("SingleEndSamples",
                              contains = "TargetedExperimet",
                              slots = c(
                                files = "character",
                                output = "character"
                              ))

#' Paired sample class
#'
#' This class represent a standard paired-end experiment.
#' In addition to the [TargetedExperimet][.TargetedExperimet]
#' class it provides a vector of forward and reverse file paths
#' togheter with a list of forward and reverse output
#' paths. All slots can be used in the context
#' of other methods such as [by_sample()], for
#' sample wise computation and [by_run()] for run wise
#' computation.
#'
#' Instances of this class should be created using the
#' [getPairedFromFilename()] function.
#'
#' @slot forward character vector. Forward file paths
#' @slot reverse character vector. Reverse file paths
#' @slot forward.out character vector. Forward output files
#' @slot reverse.out character vector. Reverse output files
#'
#' @return a paired experiment representation
#'
#' @examples
.PairedSamples <- setClass("PairedSamples",
                           contains = "TargetedExperimet",
                           slots = c(
                             forward = "character",
                             reverse = "character",
                             forward.out = "character",
                             reverse.out = "character"
                           ))


#' Validation function for [TergetedExperiment]
#'
#' @param object a [TergetedExperiment]
#'
#' @return logical. TRUE if the object is valid
validTargetedExperimetObject <- function(object){
  if(length(object@run) != length(object@samples)){
    "Runs and samples differ in length"
  }else if(length(object@step) > 1){
    "Step must be a single string"
  }else if(length(object@n) != 1){
    "n must be a single integer number"
  }else if(length(object@samples) != object@n){
    "Number of samples is not equal to n"
  }else{
    TRUE
  }
}
# Set validity function
setValidity("TargetedExperimet", validTargetedExperimetObject)

#' Validation function for [PairedSamples]
#'
#' @param object a [PairedSamples]
#'
#' @return logical. TRUE if the object is valid
validPairedSamplesObject <- function(object){
  nforward <- length(object@forward)
  nreverse <- length(object@reverse)

  nforward.out <- length(object@forward.out)
  nreverse.out <- length(object@reverse.out)

  n <- object@n

  len <- c(nforward, nreverse,
           nforward.out, nreverse.out)

  if(all(len == n)){
    TRUE
  }else{
    "Forward, reverse and samples differ in length!"
  }
}
# Set validity function
setValidity("PairedSamples", validPairedSamplesObject)

#' Validation function for [SingleEndSamples]
#'
#' @param object a [SingleEndSamples]
#'
#' @return logical. TRUE if the object is valid
validSingleEndSamplesObject <- function(object){
  n <- object@n

  nfiles <- length(object@files)
  nout <- length(object@output)

  if(all(c(nfiles, nout) == n)){
    TRUE
  }else{
    "Files and samples differ in length"
  }
}
# Set validity
setValidity("SingleEndSamples", validSingleEndSamplesObject)


setMethod("show", signature = "TargetedExperimet", function(object){
  cat("An object of class", class(object), "\n", sep = " ")
  cat(" ", object@n, " samples\n", sep = "")
  cat(" ", length(unique(object@run)), " run/s", sep = "")
})

setMethod("show", signature = "SingleEndSamples", function(object){
  callNextMethod(object)
  outSet <- all(!is.na(object@files))
  cat("\n")
  cat(" ", object@n, " single-end files\n", sep = "")
  cat(" output set ", outSet, sep = "")
})

setMethod("show", signature = "PairedSamples", function(object){
  callNextMethod(object)
  outSet <- all(!is.na(object@forward.out)) | all(!is.na(object@reverse.out))
  cat("\n")
  cat(" ", object@n, " paired-end files ", "(", object@n*2, " files)\n", sep = "")
  cat(" output set ", outSet, sep = "")
})

