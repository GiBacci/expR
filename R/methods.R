#### Constructors ########################

#' Constructor for [Experiment].
#'
#' This functiona is used implicitly by other constructors
#' and should not be used directly.
#'
#' @param samples character vector. Sample names
#' @param run character vector. Run ids
#' @param step character. The name of the analysis step
#'
#' @return an instance of [Experiment]
.getExperiment <- function(samples, run, step){
  Experiment(samples = samples, run = run, step = step,
             n = length(samples))
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
#' @param run it can be NULL, a character vector, or a funciton.
#'  If it's NULL then no run will be assigned. If it's a character
#'  vector run will be called using the names
#'  provided. If it's a function file names will be processed
#'  using the function provided and samples will have the resulting
#'  names.
#' @param output similar to samples parameter
#'
#' @return a Single-end experiment
#'
#' @examples
#' @seealso [formatSample()] and [formatOutptut()] for sample and output
#'  formatting.
.getSingleEndFromFilename <- function(maindir = ".",
                                     pattern = ".fastq",
                                     samples = NULL,
                                     run = NULL,
                                     output = NULL,
                                     recursive = F){

  files <- sort(list.files(path = maindir,
                           pattern = pattern,
                           recursive = recursive,
                           full.names = T))

  if(length(files) == 0){
    stop("No file/s were found")
  }

  if(anyDuplicated(files)){
    stop("Duplicated input file/s")
  }

  if(is.null(run)){
    r <- as.character(rep(NA, length(files)))
  }else if(is.character(run)){
    if(length(run) != length(files)){
      stop("Run and files differ in length")
    }
    r <- run
  }else if(is.function(run)){
    r <- run(files)
  }else{
    stop("Cannot recognize run option")
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

  if(anyDuplicated(s)){
    warning("Duplicated sample names", call. = F, noBreaks. = T)
  }

  experiment <- .getExperiment(samples = s, step = "raw", run = r)
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
#'
#' @examples
#' @seealso [formatSample()], [formatOutptut()], and [runFromHeader()]
#'  for sample, output, and run formatting.
.getPairedFromFilename <- function(maindir = ".",
                                   forward = "_R1_",
                                   reverse = "_R2_",
                                   samples = NULL,
                                   run = NULL,
                                   output = NULL,
                                   recursive = F){

  fwr <- .getSingleEndFromFilename(maindir = maindir,
                                   pattern = forward,
                                   samples = samples,
                                   recursive = recursive,
                                   run = run,
                                   output = output)

  rev <- suppressWarnings(
    .getSingleEndFromFilename(maindir = maindir,
                              pattern = reverse,
                              samples = samples,
                              recursive = recursive,
                              run = run,
                              output = output))

  if(all(fwr@samples == rev@samples)){
    if(all(fwr@run == rev@run)){
      experiment <- .getExperiment(samples = fwr@samples, step = "raw",
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

#' Construct an Experiment from file names
#'
#' This methods constructs an [Experiment] object
#' from file names in a given directory.
#'
#' @param maindir Character, the path to the main directory
#' @param pattern Character or Character vector of length two.
#'  If a single character, a single-end experiment will be created
#'  otherwise the first element will be used to detected forward
#'  sequence files and the second element will be used to detect
#'  reverse sequence files.
#' @param samples either \code{NULL}, a Character vector, or a
#'  function. If \code{NULL}, samples will have the same name
#'  as the input files. If a vector, it will be used as vector
#'  name for samples. If a function, input files will be converted
#'  to sample names using the function provided.
#' @param run either \code{NULL}, a Character vector, or a
#'  function. In the first case samples will hove no run assigned,
#'  whereas in the other two the behaviour will be the same of the
#'  \code{sample} option.
#' @param output either \code{NULL}, a Character vector, or a
#'  function. Same behaviour as the \code{run} parameter.
#' @param recursive logical, should the listing recurse into directories?
#'
#' @return an [Experiment] object
#' @export
#'
#' @examples
#' @seealso [formatSample()], [formatOutptut()], and [runFromHeader()]
#' or [runFromBasedir()] for sample, output, and run formatting.
expFromFiles <-  function(maindir = ".",
                          pattern = ".fastq",
                          samples = NULL,
                          run = NULL,
                          output = NULL,
                          recursive = F){
  if(length(pattern) == 1){
    .getSingleEndFromFilename(maindir = maindir,
                              pattern = pattern,
                              samples = samples,
                              run = run,
                              output = output,
                              recursive = recursive)
  }else if(length(pattern) == 2){
    .getPairedFromFilename(maindir = maindir,
                           forward = pattern[1],
                           reverse = pattern[2],
                           samples = samples,
                           run = run,
                           output = output,
                           recursive = recursive)
  }else{
    stop("Pattern length must be 1 or 2")
  }
}

#' Task constructor
#'
#' It creates a [Task] object given
#' the experiment that will be
#' used as input.
#'
#' @param exp [Experiment], the experiment on which the
#'  task will be executed
#'
#' @return a [Task] object
#' @export
#'
#' @examples
newTask <- function(exp){
  Task(exp = exp, out = list())
}


## ACCESSORS #################################
#' Set of methods to access or to set
#' slots of a Experiment
#'
#' @param obj [Experiment]
#' @param value the value to set. For the
#'  \code{setOutput} function, value must be a
#'  character vector (or a matrix for paired samples)
#'  or a function
#' @name accessors
NULL

#' Step of Analysis
#'
#' \code{step} returns/set the step of analysis
#' from a [Experiment] object.
#'
#' @export
#'
#' @rdname accessors
setGeneric("step", function(obj) standardGeneric("step"))
setMethod("step", "Experiment", function(obj){
  return(obj@step)
})

#' @usage step(obj) <- value
#'
#' @export
#' @rdname accessors
setGeneric("step<-", function(obj, value) standardGeneric("step<-"))
setMethod("step<-", "Experiment", function(obj, value){
  obj@step <- value
  obj
})

#' Step of Analysis
#'
#' \code{setStep} sets the step of analysis
#' of an [Experiment] object.
#'
#' @export
#' @rdname accessors
setGeneric("setStep", function(obj, value) standardGeneric("setStep"))
setMethod("setStep", "Experiment", function(obj, value){
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
setMethod("N", "Experiment", function(obj){
  return(obj@n)
})


#' Get input files
#'
#' \code{files} and \code{output} return the input/output
#' files from a [Experiment] object. If object
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
#' a Experiment object. If object
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

#' Set the output of a given experiment.
#'
#' \code{setOutput} is used ot change tho output of
#' a given experiment. The \code{value} parameter must be a
#' character vector (for single end experiments)
#' or a character matrix with two columns
#' (for paired end samples). In alternative
#' it can be a function that will be used for
#' building output names.
#'
#' @export
#'
#' @rdname accessors
setGeneric("setOutput", function(obj, value) standardGeneric("setOutput"))
setMethod("setOutput", "SingleEndSamples", function(obj, value){
  if(is.character(value)){
    if(length(value) == N(obj)){
      obj@output <- value
    }else{
      stop("Output length and sample number differ")
    }
  }else if(is.function(value)){
    obj@output <- value(obj@files)
  }else{
    stop("Output must be a character vector or a function")
  }
  return(obj)
})

setMethod("setOutput", "PairedSamples", function(obj, value){
  if(is.character(value)){
    value <- as.matrix(value)
    if((ncol(value) == 2) &
       (nrow(value) == N(obj))){
      obj@forward.out <- value[,1]
      obj@reverse.out <- value[,2]
    }else{
      stop("Output dimensions are different than sample number")
    }
  }else if(is.function(value)){
    obj@forward.out <- value(obj@forward)
    obj@reverse.out <- value(obj@reverse)
  }else{
    stop("Output must be a character matrix or a function")
  }
  return(obj)
})

#' Get samples from experiment
#'
#' \code{samples} returns the samples from
#' a [Experiment] object.
#'
#' @export
#'
#' @rdname accessors
setGeneric("samples", function(obj) standardGeneric("samples"))
setMethod("samples", "Experiment", function(obj){
  obj@samples
})

#' Get run from experiment
#'
#' \code{run} returns the samples from
#' a [Experiment] object.
#'
#' @export
#'
#' @rdname accessors
setGeneric("run", function(obj) standardGeneric("run"))
setMethod("run", "Experiment", function(obj){
  obj@run
})

#' Get the number of run from experiment
#'
#' \code{nrun} returns the number of run from
#' a [Experiment] object.
#'
#' @export
#'
#' @rdname accessors
setGeneric("nrun", function(obj) standardGeneric("nrun"))
setMethod("nrun", "Experiment", function(obj){
  length(unique(obj@run))
})

#' Check if output exists
#'
#' \code{output.exists} returns \code{TRUE}
#' for each output file that is found
#' on the hard drive. If the experiment
#' is a paired-end experiment,
#' the method will return \code{TRUE} only
#' if both forward and reverse
#' output files exist.
#'
#' @param obj [Experiment]
#'
#' @export
#'
#' @rdname accessors
setGeneric("output.exists", function(obj) standardGeneric("output.exists"))
setMethod("output.exists", "SingleEndSamples", function(obj){
  file.exists(output(obj))
})

setMethod("output.exists", "PairedSamples", function(obj){
  M <- as.matrix(output(obj))
  apply(M, 1, function(x) all(file.exists(x)))
})

#' Check if input exists
#'
#' \code{input.exists} returns \code{TRUE}
#' for each input file that is found
#' on the hard drive. If the experiment
#' is a paired-end experiment,
#' the method will return \code{TRUE} only
#' if both forward and reverse
#' input files exist.
#'
#' @param obj [Experiment]
#'
#' @export
#'
#' @rdname accessors
setGeneric("input.exists", function(obj) standardGeneric("input.exists"))
setMethod("input.exists", "SingleEndSamples", function(obj){
  file.exists(files(obj))
})

setMethod("input.exists", "PairedSamples", function(obj){
  M <- as.matrix(files(obj))
  apply(M, 1, function(x) all(file.exists(x)))
})

#' Check if input and output files are
#' up to date
#'
#' \code{up2date} returns \code{TRUE}
#' for each input and output file up to date.
#' Namely, the timestamp of each input file
#' must be smaller than its output file.
#'
#' @param obj [Experiment]
#'
#' @export
#'
#' @rdname accessors
setGeneric("up2date", function(obj) standardGeneric("up2date"))
setMethod("up2date", "SingleEndSamples", function(obj){
  res <- file.mtime(files(obj)) < file.mtime(output(obj))
  res[is.na(res)] <- FALSE
  res
})

setMethod("up2date", "PairedSamples", function(obj){
  I <- as.matrix(files(obj))
  O <- as.matrix(output(obj))

  I.f <- file.mtime(I[,1])
  I.r <- file.mtime(I[,2])

  O.f <- file.mtime(O[,1])
  O.r <- file.mtime(O[,2])

  res <- I.f < O.f & I.r < O.r
  res[is.na(res)] <- FALSE
  res
})

### FUNCTION METHODS #################
#' Set of methods to apply a function
#' over an Experiment object
#'
#' If threads > 1 and the [doMC] package
#' has been installed, multiple processes will
#' be executed at the same time.
#'
#' @param obj [Experiment]
#' @param fun function, the function that will be applied
#'  to each subset of the [Experiment]
#' @param by a factor. The experiment will be spitted according
#'  to this factor
#' @param ... additional parameters to fun
#' @param thread integer, the number of parallel process
#'  to execute
#'
#' @return a list. The length of the list depends on the
#'  factor used to split the experiment. \code{byAll} returns
#'  a single object.
#' @name byMethods
NULL

#' Mothod used internally to eval
#' function
#'
#' @param obj a [Experiment] object
#' @param fct a factor
#' @param f a name object. Function parameters
#' @param thread numeric. The number of thread
#'
#' @return a list of bject returned by the function
.evalFunBy <- function(obj, fct, f, thread){
  if(!is.factor(fct)){
    fct <- as.factor(fct)
  }

  data <- as.data.frame(obj)
  list <- split(data, fct)

  if(thread > 1 & requireNamespace("foreach", quietly = TRUE) &
     requireNamespace("doParallel", quietly = TRUE) &
     requireNamespace("doMC", quietly = TRUE)){

    `%dopar%` <- foreach::`%dopar%`

    old <- foreach::getDoParWorkers()
    doMC::registerDoMC(thread)
    on.exit(doMC::registerDoMC(old))

    foreach::foreach(i = seq_along(list)) %dopar% {
      d <- list[[i]]
      eval(f, envir = d, enclos = .GlobalEnv)
    }
  } else {

    if(thread > 1){
      warning("Cannot execute parallel processes due to the
              absence of one (or more) packages needed")
    }

    lapply(list, function(d) {
      eval(f, envir = d, enclos = .GlobalEnv)
    })

  }
}

#' Run a function for each sample
#'
#' \code{bySample} the function is called splitting
#'  the experiment into sample/s.
#'
#' @export
#' @rdname byMethods
setGeneric("bySample", function(obj, fun, ..., thread = 1) standardGeneric("bySample"))
setMethod("bySample", "Experiment", function(obj, fun, ..., thread = 1){
  f <- substitute(fun(...))
  .evalFunBy(obj, obj@samples, f, thread)
})

#' Run a function for each run
#'
#' \code{byRun} the function is called splitting
#'  the experiment into run/s.
#'
#' @export
#'
#' @rdname byMethods
setGeneric("byRun", function(obj, fun, ..., thread = 1) standardGeneric("byRun"))
setMethod("byRun", "Experiment", function(obj, fun, ..., thread = 1){
  f <- substitute(fun(...))
  .evalFunBy(obj, obj@run, f, thread)
})

#' Run a function splitting the experiment
#' based on a factor
#'
#' \code{byFct} the function is called splitting
#'  the experiment according to the factor provided.
#'
#' @export
#'
#' @rdname byMethods
setGeneric("byFct", function(obj, fun, by, ..., thread = 1) standardGeneric("byFct"))
setMethod("byFct", "Experiment", function(obj, fun, by, ..., thread = 1){
  f <- substitute(fun(...))
  .evalFunBy(obj, by, f, thread)
})


#' Run a function on the whole experiment
#'
#' \code{byAll} the function is called using all files
#'  as input/output.
#'
#' @export
#'
#' @rdname byMethods
setGeneric("byAll", function(obj, fun, ..., thread = 1) standardGeneric("byAll"))
setMethod("byAll", "Experiment", function(obj, fun, ..., thread = 1){
  f <- substitute(fun(...))
  by <- rep("all", N(obj))
  res <- .evalFunBy(obj, by, f, thread)
  res[[1]]
})

#' Run a function on each file/s
#'
#' \code{byFile} the function is called on each
#'  file separately.
#'
#' @export
#'
#' @rdname byMethods
setGeneric("byFile", function(obj, fun, ..., thread = 1) standardGeneric("byFile"))
setMethod("byFile", "Experiment", function(obj, fun, ..., thread = 1){
  f <- substitute(fun(...))
  by <- paste0("file_", seq(1, N(obj)))
  by <- factor(by, levels = by)
  .evalFunBy(obj, by, f, thread)
})

### SYSTEM FUNCTION METHODS #################
#' Set of methods to apply a system call
#' over an Experiment object
#'
#' If threads > 1 and the [doMC] package
#' has been installed, multiple processes will
#' be executed at the same time.
#'
#' @param obj [Experiment]
#' @param cmd character, a system command
#' @param by factor, the experiment will be spitted according
#'  to this factor
#' @param ... argument/s to the command
#' @param stdout passed to stdout of [system2]
#' @param stderr passed to stderr of [system2]
#' @param thread integer, the number of parallel process to execute
#'
#' @return a list of bject returned by the function
#'
#' @name sysByMethods
NULL

#' Launches a system call for each sample
#'
#' \code{sysBySample} the function is called splitting
#'  the experiment into sample/s.
#'
#' @export
#' @rdname sysByMethods
setGeneric("sysBySample", function(obj, cmd, ..., stdout = "", stderr = "", thread = 1)
  standardGeneric("sysBySample"))
setMethod("sysBySample", "Experiment", function(obj, cmd, ..., stdout = "",
                                                       stderr = "", thread = 1){
  f <- substitute(system2(command = cmd, args = c(...),
                          stdout = stdout, stderr = stderr))
  .evalFunBy(obj, obj@samples, f, thread)
})


#' Launches a system call for each run
#'
#' \code{sysByRun} the function is called splitting
#'  the experiment into run/s.
#' @export
#'
#' @rdname sysByMethods
setGeneric("sysByRun", function(obj, cmd, ..., stdout = "", stderr = "", thread = 1)
  standardGeneric("sysByRun"))
setMethod("sysByRun", "Experiment", function(obj, cmd, ..., stdout = "",
                                                      stderr = "", thread = 1){
  f <- substitute(system2(command = cmd, args = c(...),
                          stdout = stdout, stderr = stderr))
  .evalFunBy(obj, obj@run, f, thread)
})

#' Launches a system call by splitting the experiment
#' based on a factor
#'
#' \code{sysByFct} the function is called splitting
#'  the experiment according to the factor provided.
#'
#' @export
#'
#' @rdname sysByMethods
setGeneric("sysByFct", function(obj, cmd, by, ..., stdout = "", stderr = "", thread = 1)
  standardGeneric("sysByFct"))
setMethod("sysByFct", "Experiment", function(obj, cmd, by, ..., stdout = "",
                                                  stderr = "", thread = 1){
  f <- substitute(system2(command = cmd, args = c(...),
                          stdout = stdout, stderr = stderr))
  .evalFunBy(obj, by, f, thread)
})

#' Run a system call on the whole experiment
#'
#' \code{sysByAll} the function is called using all files
#'  as input/output.
#'
#' @export
#'
#' @rdname sysByMethods
setGeneric("sysByAll", function(obj, cmd, ..., stdout = "", stderr = "", thread = 1)
  standardGeneric("sysByAll"))
setMethod("sysByAll", "Experiment", function(obj, cmd, ..., stdout = "",
                                                    stderr = "", thread = 1){
  f <- substitute(system2(command = cmd, args = c(...),
                          stdout = stdout, stderr = stderr))

  by <- rep("all", N(obj))
  res <- .evalFunBy(obj, by, f, thread)
  res[[1]]
})

#' Run a system call on each file of the experiment
#'
#' \code{sysByFile} the function is called on each
#'  file separately.
#'
#' @export
#'
#' @rdname sysByMethods
setGeneric("sysByFile", function(obj, cmd, ..., stdout = "", stderr = "", thread = 1)
  standardGeneric("sysByFile"))
setMethod("sysByFile", "Experiment", function(obj, cmd, ..., stdout = "",
                                             stderr = "", thread = 1){
  f <- substitute(system2(command = cmd, args = c(...),
                          stdout = stdout, stderr = stderr))

  by <- paste0("file_", seq(1, N(obj)))
  by <- factor(by, levels = by)

  .evalFunBy(obj, by, f, thread)
})

#' Get an experiment from the output of another one.
#'
#' This function will convert an experiment to another one
#' by inverting the output files as input. The new output files
#' will be generate depending on the [new.output] param
#'
#' @param obj a Experiment. The experiment to convert
#' @param new.output either a character vector (or a matrix with two
#'  columns for paired experiments), NULL, or a function.
#'  If a character vector is provided, output file will be named
#'  accordingly otherwise the specified function will be used
#'  to convert input files into output ones. If NULL, the
#'  resulting experiment will hove no output.
#'
#' @return a Experiment
#' @export
#'
#' @examples
setGeneric("out2exp", function(obj, new.output = NULL)
  standardGeneric("out2exp"))
setMethod("out2exp", "SingleEndSamples", function(obj, new.output = NULL){
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

setMethod("out2exp", "PairedSamples", function(obj, new.output = NULL){
  if(is.null(new.output)){
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


#' Removal of output file/s from a [Experiment]
#' object
#'
#' This method will remove any output file/s
#' stored in the given in the [Experiment]
#' object
#'
#' @param obj a [Experiment]
#'
#' @return \code{TRUE} for each file correctly
#'  removed
#'
#' @export
#'
#' @examples
setGeneric("clearOutput", function(obj) standardGeneric("clearOutput"))
setMethod("clearOutput", "SingleEndSamples", function(obj){
  v <- unique(as.vector(output(obj)))
  v <- v[file.exists(v)]
  r <- file.remove(v)

  names(r) <- v
  r
})
setMethod("clearOutput", "PairedSamples", function(obj){
  v <- unique(as.vector(as.matrix(output(obj))))
  v <- v[file.exists(v)]
  r <- file.remove(v)

  names(r) <- v
  r
})

## SUBSETTING
setMethod("[", "Experiment", function(x, i, drop="missing"){
  .samples = x@samples[i]
  .run = x@run[i]
  .n = length(.samples)
  Experiment(x, samples = .samples, run = .run, n = .n)
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

setMethod("length", "Experiment", function(x){
  N(x)
})

#' Convert an experiment into a data frame
#'
#' @param x [Experiment]
#' @param row.names ignored
#' @param optional ignored
#' @param ... ignored
#'
#' @return a data frame
#' @export
#'
#' @examples
as.data.frame.Experiment <- function(x, row.names=NULL, optional=FALSE, ...){
  data.frame(samples = x@samples, run = x@run, step = rep(x@step, x@n),
             stringsAsFactors = F)
}

#' @rdname as.data.frame.Experiment
#' @export
as.data.frame.SingleEndSamples <- function(x, row.names=NULL, optional=FALSE, ...){
  data.frame(as.data.frame.Experiment(x), files = x@files,
             output = x@output, stringsAsFactors = F)
}

#' @rdname as.data.frame.Experiment
#' @export
as.data.frame.PairedSamples <- function(x, row.names=NULL, optional=FALSE, ...){
  data.frame(as.data.frame.Experiment(x), forward = x@forward,
             reverse = x@reverse, forward.out = x@forward.out,
             reverse.out = x@reverse.out, stringsAsFactors = F)
}

##### TASK METHODS

#' Run a Task on the experiment
#'
#' a Task (namely a function) will be
#' called using the experiment as first
#' input.
#'
#' @param obj a [Task] object
#' @param task a function. The first argument of
#'  this function will be the [Experiment]
#'  object contained in the [Task]
#' @param ... additianl arguments to be passed
#'  to the funciton call
#' @param check.out if \code{TRUE} the [Task] won't be
#'  executedon existing output file/s
#' @param quiet if \code{TRUE} no message will be
#'  displayed
#'
#' @return a [Task] object
#' @export
#'
#' @examples
setGeneric("runTask", function(obj, task, ..., check.out = F, quiet = T) standardGeneric("runTask"))
setMethod("runTask", signature = "Task", function(obj, task, ..., check.out = F, quiet = T){

  exp <- obj@exp
  out <- obj@out
  index <- length(out) + 1

  done <- output.exists(exp) & up2date(exp)

  if(check.out & sum(done) != 0){
    if(!quiet){
      message("Some samples have been already processed:")
      s <- paste(samples(exp)[done], collapse = "\n")
      message(s)

      v <- as.vector(as.matrix(output(exp)[done,]))
      v <- paste(v[file.exists(v)], collapse = "\n")
      message("Remove the following files if you want to repeat the task:")
      message(v)
    }
    e <- exp[!done]
  }else{
    e <- exp
  }

  s <- substitute(task(e, ...))
  out[[index]] <- eval(s, enclos = .GlobalEnv)

  obj@out <- out
  return(obj)
})


#' Change the experimen withing a [Task] object
#'
#' The experiment contained in a [Task] object
#' will be changed according to the \code{new}
#' parameter.
#'
#' @param obj a [Task] object
#' @param new either a new [Experiment]
#'  or a function that will be used to process
#'  the [Experiment] contained in the Task
#' @param ... additional argument that will be
#'  passed to \code{new} if it is function
#'  (ignored otherwise).
#'
#' @return a [Task] object
#' @export
#'
#' @examples
setGeneric("changeExp", function(obj, new, ...) standardGeneric("changeExp"))
setMethod("changeExp", "Task", function(obj, new, ...){
  if(is.function(new)){
    exp <- obj@exp
    s <- substitute(new(exp, ...))
    obj@exp <- eval(s)
  }else{
    obj@exp <- new
  }
  return(obj)
})
## TASK OUTPUT ACCESSORS
#' Set of methods to access the
#' slots of a [TaskOutput] object
#'
#' @param obj a TaskOutput object
#'
#' @name task.accessors
NULL


#' Get experiment from a [Task] object
#'
#' \code{getExp} returns the experiment from
#' a [Task] object
#'
#' @param obj a [Task] object.
#'
#' @return the Experiment
#' @export
#'
#' @rdname task.accessors
#' @examples
setGeneric("getExp", function(obj) standardGeneric("getExp"))
setMethod("getExp", "Task", function(obj){
  obj@exp
})

#' Gte the output from a [Task] object
#'
#' \code{getOut} returns the output
#'  from a [Task] object
#'
#' @param obj
#'
#' @return the output from a [Task] object
#'
#' @export
#'
#' @rdname task.accessors
#' @examples
setGeneric("getOut", function(obj) standardGeneric("getOut"))
setMethod("getOut", "Task", function(obj){
  obj@out
})
