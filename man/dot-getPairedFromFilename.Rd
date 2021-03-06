% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{.getPairedFromFilename}
\alias{.getPairedFromFilename}
\title{Constructor for paired-end experiments}
\usage{
.getPairedFromFilename(maindir = ".", forward = "_R1_",
  reverse = "_R2_", samples = NULL, run = NULL, output = NULL,
  recursive = F)
}
\arguments{
\item{maindir}{character. The directory in which to search for
sequence files.}

\item{forward}{regular expression. Only file names which match
the regular expression will be returned as forward files.}

\item{reverse}{regular expression. Only file names which match
the regular expression will be returned as reverse files.}

\item{samples}{it can be NULL, a character vector, or a funciton.
If it's NULL then file names will be used as sample names. If
it's a character vector sample will be called using the names
provided. If it's a function file names will be processed
using the function provided and samples will have the resulting
names.}

\item{run}{one of "none", "folder", or "header". If "none" no
run will be specified, if "folder" the run id will be inferred from
the folder structure, and if "header" the run id will be parsed from
the header of the sequences.}

\item{output}{similar to samples parameter}

\item{recursive}{logical. Should the listing recurse into directories?}
}
\value{
a paired-end experiment
}
\description{
This funciton is similar to \link{getSingleEndFromFilename} but
it can parse forward and reverse files separately.
}
\seealso{
\code{\link[=formatSample]{formatSample()}}, \code{\link[=formatOutptut]{formatOutptut()}}, and \code{\link[=runFromHeader]{runFromHeader()}}
for sample, output, and run formatting.
}
