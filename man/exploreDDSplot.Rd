% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exploreDDSplot.R
\name{exploreDDSplot}
\alias{exploreDDSplot}
\title{exploreDDSplot}
\usage{
exploreDDSplot(
  countMatrix,
  targets,
  cmp = cmp[[1]],
  preFilter = NULL,
  samples,
  blind = TRUE,
  scattermatrix = FALSE,
  plotly = FALSE,
  savePlot = FALSE,
  filePlot = NULL
)
}
\arguments{
\item{countMatrix}{\code{date.frame} or \code{matrix} containing raw read counts}

\item{targets}{targets \code{data.frame}}

\item{cmp}{\verb{character matrix} where comparisons are defined in two columns.
This matrix should be generated with the \code{\link[systemPipeR:readComp]{systemPipeR::readComp()}} function
from the targets file. Values used for comparisons need to match those in the
\code{Factor} column of the targets file.}

\item{preFilter}{allows removing rows in which there are very few reads.
Accepts a numeric value with the minimum of total reads to keep.
Default is \code{NULL}.}

\item{samples}{a \verb{character vector} of two samples or \code{ALL} samples in the
dataset. Could be specified the \code{SampleName} column name of the targets
file or the respective numeric values. Also, if set as \code{ALL},
a correlation matrix it will be plot.}

\item{blind}{logical, whether to blind the transformation to the experimental
design (see varianceStabilizingTransformation), from \code{\link[DESeq2:vst]{DESeq2::vst()}}
or \code{\link[DESeq2:rlog]{DESeq2::rlog()}}.}

\item{scattermatrix}{if \code{samples} set as \code{ALL}, requires to assign \code{TRUE}
to build    a correlation matrix and plot the correlogram of all the samples.}

\item{plotly}{logical: when \code{FALSE} (default), the \code{ggplot2} plot will be
returned. \code{TRUE} returns the \code{plotly} version of the plot.}

\item{savePlot}{logical: when \code{FALSE} (default), the plot will not be saved.
If \code{TRUE} the plot will be saved, and requires the \code{filePlot} argument.}

\item{filePlot}{file name where the plot will be saved. For more information,
please consult the \code{\link[ggplot2:ggsave]{ggplot2::ggsave()}} function.}
}
\value{
returns an object of \verb{ggplot2 plot}.
}
\description{
Scatterplot of transformed counts from
two samples or grid of all samples
}
\examples{
## Targets file
targetspath <- system.file("extdata", "targets.txt", package = "systemPipeR")
targets <- read.delim(targetspath, comment = "#")
cmp <- systemPipeR::readComp(
    file = targetspath,
    format = "matrix", delim = "-"
)
## Count table file
countMatrixPath <- system.file("extdata", "countDFeByg.xls",
    package = "systemPipeR"
)
countMatrix <- read.delim(countMatrixPath, row.names = 1)
## Plot
exploreDDSplot(countMatrix, targets,
    cmp = cmp[[1]], preFilter = NULL,
    samples = c(3, 4)
)
}
