% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{task.accessors}
\alias{task.accessors}
\alias{getExp}
\alias{getOut}
\title{Set of methods to access the
slots of a \link{TaskOutput} object}
\usage{
getExp(obj)

getOut(obj, index = NULL)
}
\arguments{
\item{obj}{a TaskOutput object}

\item{index}{integer, the index of the output to get.
If \code{NULL} the entire list is retrived.}

\item{obj}{a \link{Task} object.}
}
\value{
the Experiment

the output from a \link{Task} object
}
\description{
\code{getExp} returns the experiment from
a \link{Task} object

\code{getOut} returns the output
from a \link{Task} object
}
