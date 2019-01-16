#' Constructor for [TargetedExperiment].
#'
#' This functiona is used implicitly by other constructors
#' and should not be used directly.
#'
#' @param samples character vector. Sample names
#' @param run character vector. Run ids
#' @param step character. The name of the analysis step
#'
#' @return an instance of [TargetedExperiment]
.getTargetedExperiment <- function(samples, run, step){
  .TargetedExperimet(samples = samples, run = run, step = step,
                     n = length(samples))
}


#' Get run id from sequence header
#'
#' This function is used internally by other
#' constructors to detect the run id of each file.
#' The id is parsed from the sequence header.
#'
#' @param path character vector. The paths to the sequence file/s
#' @param nrecords numeric. The number of record that must be read
#'  before returning the run id
#'
#' @return the run id as character
.getRunFromHeader <- function(path, nrecords = 1){
  sapply(path, function(p){
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


#' Constructor for single-end experiments
#'
#' This function will return a single end experiment based
#' on the name of the files selected.
#'
#' @param maindir character. The directory in which to search for
#'  sequence files.
#' @param pattern regular expression. Only file names which match
#'  the regular expression will be returned.
#' @param samples it can be NULL, a character vector, or a funciton.
#'  If it's NULL then file names will be used as sample names. If
#'  it's a character vector sample will be called using the names
#'  provided. If it's a function file names will be processed
#'  using the function provided and samples will have the resulting
#'  names.
#' @param recursive logical. Should the listing recurse into directories?
#' @param run one of "none", "folder", or "header". If "none" no
#'  run will be specified, if "folder" the run id will be inferred from
#'  the folder structure, and if "header" the run id will be parsed from
#'  the header of the sequences.
#' @param output similar to samples parameter
#'
#' @return a Single-end experiment
#' @export
#'
#' @examples
#' \dontrun{
#' single <- getSingleEndFromFilename(samples = formatSample(sep = "_"),
#'                                    recursive = F,
#'                                    run = "header",
#'                                    output = formatOutptut(prefix = "out_"))
#' }
#' @seealso [formatSample()] and [formatOutptut()] for sample and output
#'  formatting.
getSingleEndFromFilename <- function(maindir = ".",
                                     pattern = ".fastq",
                                     samples = NULL,
                                     recursive = F,
                                     run = c("none", "folder", "header"),
                                     output = NULL){

  run <- match.arg(run)

  files <- sort(list.files(path = maindir,
                           pattern = pattern,
                           recursive = recursive,
                           full.names = T))

  if(length(files) == 0){
    stop("No file/s were found")
  }

  if(run == "folder"){
    run <- basename(dirname(files))
  }else if(run == "header"){
    run <- .getRunFromHeader(files)
  }else{
    run <- rep("none", length(files))
  }

  if(is.null(output)){
    o <- as.character(rep(NA, length(files)))
  }else if(is.character(output)){
    if(length(output) != length(files)){
      stop("Length of output files differs from sample length")
    }
    o <- output
  }else if(is.function(output)){
    o <- output(files)
  }

  if(is.null(samples)){
    s <- basename(files)
  }else if(is.function(samples)){
    s <- samples(files)
  }else if(is.character(samples)){
    s <- samples
  }

  experiment <- .getTargetedExperiment(samples = s, step = "raw",
                                       run = run)
  .SingleEndSamples(experiment, files = files, output = o)
}


#' Constructor for paired-end experiments
#'
#' This funciton is similar to [getSingleEndFromFilename] but
#' it can parse forward and reverse files separately.
#'
#' @param maindir character. The directory in which to search for
#'  sequence files.
#' @param forward regular expression. Only file names which match
#'  the regular expression will be returned as forward files.
#' @param reverse regular expression. Only file names which match
#'  the regular expression will be returned as reverse files.
#' @param samples it can be NULL, a character vector, or a funciton.
#'  If it's NULL then file names will be used as sample names. If
#'  it's a character vector sample will be called using the names
#'  provided. If it's a function file names will be processed
#'  using the function provided and samples will have the resulting
#'  names.
#' @param recursive logical. Should the listing recurse into directories?
#' @param run one of "none", "folder", or "header". If "none" no
#'  run will be specified, if "folder" the run id will be inferred from
#'  the folder structure, and if "header" the run id will be parsed from
#'  the header of the sequences.
#' @param output similar to samples parameter
#'
#' @return a paired-end experiment
#' @export
#'
#' @examples
getPairedFromFilename <- function(maindir = ".",
                                  forward = "_R1_",
                                  reverse = "_R2_",
                                  samples = NULL,
                                  recursive = F,
                                  run = c("none", "folder", "header"),
                                  output = NULL){

  run <- match.arg(run)

  fwr <- getSingleEndFromFilename(maindir = maindir,
                                  pattern = forward,
                                  samples = samples,
                                  recursive = recursive,
                                  run = run,
                                  output = output)

  rev <- getSingleEndFromFilename(maindir = maindir,
                                  pattern = reverse,
                                  samples = samples,
                                  recursive = recursive,
                                  run = run,
                                  output = output)

  if(all(fwr@samples == rev@samples)){
    if(all(fwr@run == rev@run)){
      experiment <- .getTargetedExperiment(samples = fwr@samples, step = "raw",
                                           run = fwr@run)
      .PairedSamples(experiment, forward = fwr@files, reverse = rev@files,
                     forward.out = fwr@output, reverse.out = rev@output)
    }else{
      stop("Forward and reverse runs are different")
    }
  }else{
    stop("Forward and reverse samples are different")
  }
}

#' Set of methods to access or to set
#' slots of a TargetedExperiment
#'
#' @param obj a TargetedExperiment object
#' @param value the value to set
#' @name accessors
NULL

#' Step of Analysis
#'
#' \code{step} returns/set the step of analysis
#' from a [TargetedExperiment] object.
#'
#' @export
#'
#' @rdname accessors
setGeneric("step", function(obj) standardGeneric("step"))
setMethod("step", "TargetedExperimet", function(obj){
  return(obj@step)
})

#' @usage step(obj) <- value
#'
#' @export
#' @rdname accessors
setGeneric("step<-", function(obj, value) standardGeneric("step<-"))
setMethod("step<-", "TargetedExperimet", function(obj, value){
  obj@step <- value
  obj
})


#' Number of samples
#'
#' \code{N} returns the number of samples.
#'
#' @export
#'
#' @rdname accessors
setGeneric("N", function(obj) standardGeneric("N"))
setMethod("N", "TargetedExperimet", function(obj){
  return(obj@n)
})


#' Get input files
#'
#' \code{files} and \code{output} return the input/output
#' files from a [TargetedExperiment] object. If object
#' is a [SingleEndSamples] it will return a
#' vector of files otherwise it will return
#' a data frame with forward and reverse files.
#'
#' @export
#'
#' @rdname accessors
setGeneric("files", function(obj) standardGeneric("files"))
setMethod("files", "SingleEndSamples", function(obj){
  obj@files
})
setMethod("files", "PairedSamples", function(obj){
  data.frame(forward = obj@forward,
             reverse = obj@reverse,
             stringsAsFactors = F)
})

#' This method return the output files from
#' a TargetedExperiment object. If object
#' is a [SingleEndSamples] it will return a
#' vector of files otherwise it will return
#' a data frame with forward and reverse files.
#'
#' @export
#'
#' @rdname accessors
setGeneric("output", function(obj) standardGeneric("output"))
setMethod("output", "SingleEndSamples", function(obj){
  obj@output
})
setMethod("output", "PairedSamples", function(obj){
  data.frame(forward.out = obj@forward.out,
             reverse.out = obj@reverse.out,
             stringsAsFactors = F)
})

#' Get samples from experiment
#'
#' \code{samples} returns the samples from
#' a [TargetedExperiment] object.
#'
#' @export
#'
#' @rdname accessors
setGeneric("samples", function(obj) standardGeneric("samples"))
setMethod("samples", "TargetedExperimet", function(obj){
  obj@samples
})

#' Get run from experiment
#'
#' \code{run} returns the samples from
#' a [TargetedExperiment] object.
#'
#' @export
#'
#' @rdname accessors
setGeneric("run", function(obj) standardGeneric("run"))
setMethod("run", "TargetedExperimet", function(obj){
  obj@run
})


#' Mothod used internally to eval
#' function
#'
#' @param obj a [TargetedExperiment] object
#' @param fct a factor
#' @param f a name object. Function parameters
#' @param thread numeric. The number of thread
#'
#' @return a list of bject returned by the function
#' @export
.evalFunBy <- function(obj, fct, f, thread){
  if(!is.factor(fct)){
    fct <- as.factor(fct)
  }

  data <- as.data.frame(obj)
  list <- split(data, fct)

  if(thread > 1 & requireNamespace("foreach", quietly = TRUE) &
     requireNamespace("doParallel", quietly = TRUE)){

    `%dopar%` <- foreach::`%dopar%`

    cl <- parallel::makeCluster(thread)
    parallel::clusterExport(cl = cl, varlist = ls(.GlobalEnv))
    doParallel::registerDoParallel(cl, cores = thread)
    on.exit(parallel::stopCluster(cl))

    foreach::foreach(i = seq_along(list)) %dopar% {
      d <- list[[i]]
      eval(f, envir = d, enclos = .GlobalEnv)
    }
  } else {

    lapply(list, function(d) {
      eval(f, envir = d, enclos = .GlobalEnv)
    })

  }
}

#' Mothod used internally to eval
#' a system call
#'
#' @param list a list of environments (data frame)
#' @param command the command to be executed
#' @param args a name object. Thsi will be evalueted
#'  and passed to the args param of [system2] function
#' @param stdout passed to stdout of [system2]
#' @param stderr passed to stderr of [system2]
#' @param thread numeric. The number of thread
#'
#' @return a list of bject returned by [system2] function
#' @export
.evalCommandList <- function(list, command, args, stdout, stderr, thread){
  print(ls(.GlobalEnv))
  if(thread > 1 & requireNamespace("foreach", quietly = TRUE) &
     requireNamespace("doParallel", quietly = TRUE)){

    # foreach <- foreach::foreach
    `%dopar%` <- foreach::`%dopar%`

    cl <- parallel::makeCluster(thread)
    parallel::clusterExport(cl = cl, varlist = ls(.GlobalEnv))
    doParallel::registerDoParallel(cl, cores = thread)
    on.exit(parallel::stopCluster(cl))

    foreach::foreach(i = seq_along(list)) %dopar% {
      d <- list[[i]]
      a <- eval(args, envir = d, enclos = .GlobalEnv)
      system2(command = command, args = a,
              stdout = stdout, stderr = stderr)
    }
  } else {
    lapply(list, function(d) {
      a <- eval(args, envir = d, enclos = .GlobalEnv)
      system2(command = command, args = a,
              stdout = stdout, stderr = stderr)
    })
  }
}

#' Run a function for each sample
#'
#' This function will execute a function call
#' for each sample of the experiment provided.
#' If threads > 1 and the doParallel package
#' has been installed, multiple processes will
#' be executed at the same time.
#'
#' @param obj the experiment
#' @param fun the function
#' @param ... additional parameters to fun
#' @param thread the number of parallel process to execute
#'
#' @return a list of object
#' @export
#'
#' @examples
setGeneric("bySample", function(obj, fun, ..., thread = 1) standardGeneric("bySample"))
setMethod("bySample", "TargetedExperimet", function(obj, fun, ..., thread = 1){
  f <- substitute(fun(...))
  .evalFunBy(obj, obj@samples, f, thread)
})

#' Run a function for each run
#'
#' This function will execute a function call
#' for each run of the experiment provided.
#' If threads > 1 and the doParallel package
#' has been installed, multiple processes will
#' be executed at the same time.
#'
#' @param obj the experiment
#' @param fun the function
#' @param ... additional parameters to fun
#' @param thread the number of parallel process to execute
#'
#' @return a list of object
#' @export
#'
#' @examples
setGeneric("byRun", function(obj, fun, ..., thread = 1) standardGeneric("byRun"))
setMethod("byRun", "TargetedExperimet", function(obj, fun, ..., thread = 1){
  f <- substitute(fun(...))
  .evalFunBy(obj, obj@run, f, thread)
})

#' Run a function splitting the experiment
#' based on a factor
#'
#' This function will execute a function call
#' for each set of the experiment.
#' If threads > 1 and the doParallel package
#' has been installed, multiple processes will
#' be executed at the same time.
#'
#' @param obj the experiment
#' @param fun the function
#' @param by a factor. The experiment will be spitted according
#'  to this factor
#' @param ... additional parameters to fun
#' @param thread the number of parallel process to execute
#'
#' @return a list of object
#' @export
#'
#' @examples
setGeneric("byFct", function(obj, fun, by, ..., thread = 1) standardGeneric("byFct"))
setMethod("byFct", "TargetedExperimet", function(obj, fun, by, ..., thread = 1){
  f <- substitute(fun(...))
  .evalFunBy(obj, by, f, thread)
})


#' Launches a system call for each sample
#'
#' This function will execute a system call
#' for each sample of the experiment provided.
#' If threads > 1 and the doParallel package
#' has been installed, multiple processes will
#' be executed at the same time.
#'
#' @param obj the experiment
#' @param cmd the system command
#' @param ... argument/s to the command
#' @param stdout passed to stdout of [system2]
#' @param stderr passed to stderr of [system2]
#' @param thread the number of parallel process to execute
#'
#' @return a list of object
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # Run only in unix systems
#' sysBySample(exp, "echo", samples)
#'
#' }
#' @export
setGeneric("sysBySample", function(obj, cmd, ..., stdout = "", stderr = "", thread = 1)
  standardGeneric("sysBySample"))
setMethod("sysBySample", "TargetedExperimet", function(obj, cmd, ..., stdout = "",
                                                       stderr = "", thread = 1){
  f <- substitute(system2(command = cmd, args = c(...),
                          stdout = stdout, stderr = stderr))
  .evalFunBy(obj, obj@samples, f, thread)
})


#' Launches a system call for each run
#'
#' This function will execute a system call
#' for each run of the experiment provided.
#' If threads > 1 and the doParallel package
#' has been installed, multiple processes will
#' be executed at the same time.
#'
#' @param obj the experiment
#' @param cmd the system command
#' @param ... argument/s to the command
#' @param stdout passed to stdout of [system2]
#' @param stderr passed to stderr of [system2]
#' @param thread the number of parallel process to execute
#'
#' @return a list of object
#' @export
#'
#' @examples
setGeneric("sysByRun", function(obj, cmd, ..., stdout = "", stderr = "", thread = 1)
  standardGeneric("sysByRun"))
setMethod("sysByRun", "TargetedExperimet", function(obj, cmd, ..., stdout = "",
                                                      stderr = "", thread = 1){
  f <- substitute(system2(command = cmd, args = c(...),
                          stdout = stdout, stderr = stderr))
  .evalFunBy(obj, obj@run, f, thread)
})

#' Launches a system call by splitting the experiment
#' based on a factor
#'
#' This function will execute a system call
#' for each set of the experiment.
#' If threads > 1 and the doParallel package
#' has been installed, multiple processes will
#' be executed at the same time.
#'
#' @param obj the experiment
#' @param cmd the system command
#' @param by a factor. The experiment will be spitted according
#'  to this factor
#' @param ... argument/s to the command
#' @param stdout passed to stdout of [system2]
#' @param stderr passed to stderr of [system2]
#' @param thread the number of parallel process to execute
#'
#' @return a list of object
#' @export
#'
#' @examples
setGeneric("sysByFct", function(obj, cmd, by, ..., stdout = "", stderr = "", thread = 1)
  standardGeneric("sysByFct"))
setMethod("sysByFct", "TargetedExperimet", function(obj, cmd, ..., stdout = "",
                                                  stderr = "", thread = 1){
  f <- substitute(system2(command = cmd, args = c(...),
                          stdout = stdout, stderr = stderr))
  .evalFunBy(obj, by, f, thread)
})


#' Get an experiment from the output of another one.
#'
#' This function will convert an experiment to another one
#' by inverting the output files as input. The new output files
#' will be generate depending on the [new.output] param
#'
#' @param obj a TargetedExperiment. The experiment to convert
#' @param new.output either a character vector (or a matrix with two
#'  columns for paired experiments), NULL, or a function.
#'  If a character vector is provided, output file will be named
#'  accordingly otherwise the specified function will be used
#'  to convert input files into output ones. If NULL, the
#'  resulting experiment will hove no output.
#'
#' @return a TargetedExperiment
#' @export
#'
#' @examples
setGeneric("getExperimentFromOutput", function(obj, new.output = NULL)
  standardGeneric("getExperimentFromOutput"))
setMethod("getExperimentFromOutput", "SingleEndSamples", function(obj, new.output = NULL){
  if(is.null(new.output)){
    new.o <- as.character(rep(NA, obj@n))
  }else if(is.character(new.output)){
    if(length(new.output) != obj@n){
      stop("Length of output files differs from sample length")
    }
    new.o <- new.output
  }else if(is.function(new.output)){
    new.o <- new.output(obj@files)
  }
  obj@files <- obj@output
  obj@output <- new.o
  return(obj)
})

setMethod("getExperimentFromOutput", "PairedSamples", function(obj, new.output = NULL){
  if(is.null(new.output) | is.null(dim(new.output))){
    new.o.fwr <- as.character(rep(NA, obj@n))
    new.o.rev <- as.character(rep(NA, obj@n))
  }else if(is.function(new.output)){
    new.o.fwr <- new.output(obj@forward)
    new.o.rev <- new.output(obj@reverse)
  }else if(ncol(new.output) == 2){
    m <- as.matrix(new.output)
    if(is.character(m)){
      if(nrow(m) != obj@n){
        stop("Length of output files differs from sample length")
      }
      new.o.fwr <- m[,1]
      new.o.rev <- m[,2]
    }
  }

  obj@forward <- obj@forward.out
  obj@reverse <- obj@reverse.out

  obj@forward.out <- new.o.fwr
  obj@reverse.out <- new.o.rev

  return(obj)
})

## SUBSETTING
setMethod("[", "TargetedExperimet", function(x, i, drop="missing"){
  .samples = x@samples[i]
  .run = x@run[i]
  .n = length(.samples)
  .TargetedExperimet(x, samples = .samples, run = .run, n = .n)
})

setMethod("[", "SingleEndSamples", function(x, i, drop="missing"){
  .samples = x@samples[i]
  .run = x@run[i]
  .n = length(.samples)
  .files <- x@files[i]
  .output <- x@output[i]

  .SingleEndSamples(x, samples = .samples, run = .run, n = .n,
                    files = .files, output = .output)
})

setMethod("[", "PairedSamples", function(x, i, drop="missing"){
  .samples = x@samples[i]
  .run = x@run[i]
  .n = length(.samples)

  .forward <- x@forward[i]
  .reverse <- x@reverse[i]

  .forward.out <- x@forward.out[i]
  .reverse.out <- x@reverse.out[i]

  .PairedSamples(x, samples = .samples, run = .run, n = .n,
                 forward = .forward, reverse = .reverse,
                 forward.out = .forward.out, reverse.out = .reverse.out)
})

#' Convert a Targeted experiment into a data frame
#'
#' @param x the experiment
#' @param row.names ignored
#' @param optional ignored
#' @param ... ignored
#'
#' @return a data frame
#' @export
#'
#' @examples
as.data.frame.TargetedExperimet <- function(x, row.names=NULL, optional=FALSE, ...){
  data.frame(samples = x@samples, run = x@run, step = rep(x@step, x@n),
             stringsAsFactors = F)
}

#' @rdname as.data.frame.TargetedExperimet
#' @export
as.data.frame.SingleEndSamples <- function(x, row.names=NULL, optional=FALSE, ...){
  data.frame(as.data.frame.TargetedExperimet(x), files = x@files,
             output = x@output, stringsAsFactors = F)
}

#' @rdname as.data.frame.TargetedExperimet
#' @export
as.data.frame.PairedSamples <- function(x, row.names=NULL, optional=FALSE, ...){
  data.frame(as.data.frame.TargetedExperimet(x), forward = x@forward,
             reverse = x@reverse, forward.out = x@forward.out,
             reverse.out = x@reverse.out, stringsAsFactors = F)
}
