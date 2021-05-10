################
## exploreDDS ##
################
#' @title exploreDDS
#' @description Convenience wrapper function to transform raw read counts using 
#' the [DESeq2::DESeq2-package()] package transformations methods. The input 
#' file has to contain all the genes, not just differentially expressed ones.
#'
#' @param countMatrix `date.frame` or `matrix` containing raw read counts.
#' @param targets targets `data.frame`.
#' @param cmp `character matrix` where comparisons are defined in two columns.
#' This matrix should be generated with the [systemPipeR::readComp()] function
#' from the targets file. Values used for comparisons need to match those in the
#' `Factor` column of the targets file.
#' @param preFilter allows removing rows in which there are very few reads.
#' Accepts a numeric value with the minimum of total reads to keep. Default is 
#' `NULL`.
#' @param transformationMethod a `character string` indicating which 
#' transformation method it will be used on the raw read counts. Supported 
#' methods include `rlog` and `vst` using the `DESeq2` package or default `raw`
#' for no data transformation.
#' @param blind logical, whether to blind the transformation to the experimental
#' design (see varianceStabilizingTransformation),
#' from [DESeq2::vst()] or [DESeq2::rlog()].
#'
#' @details Note that the recommendation is to use the resulting transformed
#' values in the `transformationMethod` argument only for visualization and
#' clustering, not for differential expression analysis which needs raw counts.
#'  Users are strongly encouraged to consult the [DESeq2::DESeq2-package()]
#' vignette for
#' more detailed information on this topic and how to properly run `DESeq2` on 
#' data sets with more complex experimental designs.
#'
#' @references For more details on `DESeq2`, please consult the following
#' page: 
#'\href{http://bioconductor.org/packages/release/bioc/html/DESeq2.html}{DESeq2}.
#' For more details on `targets` file definition, please consult the following
#' page: 
#' \href{http://www.bioconductor.org/packages/release/bioc/vignettes/systemPipeR/inst/doc/systemPipeR.html#25_structure_of_targets_file}{systemPipeR}.
#'
#' @author Daniela Cassol
#'
#' @return returns an object of class [DESeq2::DESeqTransform()].
#'
#' @examples
#' ## Targets file
#' targetspath <- system.file("extdata", "targets.txt", package = "systemPipeR")
#' targets <- read.delim(targetspath, comment = "#")
#' cmp <- systemPipeR::readComp(file = targetspath, format = "matrix", 
#' delim = "-")
#' ## Count table file
#' countMatrixPath <- system.file("extdata", "countDFeByg.xls", 
#' package = "systemPipeR")
#' countMatrix <- read.delim(countMatrixPath, row.names = 1)
#' ## Run
#' exploredds <- exploreDDS(countMatrix, targets,
#'     cmp = cmp[[1]],
#'     preFilter = NULL, transformationMethod = "raw"
#' )
#' exploredds
#' @export
#' @importFrom DESeq2 DESeqDataSetFromMatrix 
#' @importFrom DESeq2 counts DESeq rlog varianceStabilizingTransformation
exploreDDS <- function(countMatrix, targets, cmp = cmp[[1]],
                       preFilter = NULL, transformationMethod = "raw",
                       blind = TRUE) {
    if (!transformationMethod %in% c("raw", "rlog", "vst")) {
        stop("Supported methods include 'raw', 'rlog' and 'vst'")
    }
    if (is.data.frame(countMatrix)) {
        countMatrix <- Biostrings::as.matrix(countMatrix)
    } else if (is.matrix(countMatrix)) {
        countMatrix <- countMatrix
    } else {
        stop("countMatrix needs to be assigned an object of class 'data.frame'
             OR 'matrix'")
    }
    if (!is.data.frame(targets)) stop("targets needs to be assignes an object
                                      of class 'data.frame'")
    if (all(!is.matrix(cmp) & length(cmp) == 2)) cmp <- t(as.matrix(cmp))
    samples <- as.character(targets$Factor)
    names(samples) <- paste(as.character(targets$SampleName), "", sep = "")
    suppressWarnings({
        dds <- DESeq2::DESeqDataSetFromMatrix(
            countData = countMatrix,
            colData = data.frame(condition = samples), design = ~condition
        )
    })
    if (!is.null(preFilter)) {
        if (!is.numeric(preFilter)) stop("'preFilter' needs to be numeric 
                                         value.")
        keep <- BiocGenerics::rowSums(DESeq2::counts(dds)) >= preFilter
        dds <- dds[keep, ]
    }
    dds_deseq2 <- DESeq2::DESeq(dds, quiet = TRUE)
    ## Count data transformations
    if (transformationMethod == "rlog") {
        normdata <- DESeq2::rlog(dds_deseq2, blind = TRUE)
    } else if (transformationMethod == "vst") {
        normdata <- DESeq2::varianceStabilizingTransformation(dds_deseq2,
            blind = TRUE)
    } else if (transformationMethod == "raw") {
        normdata <- dds
    }
    return(normdata)
}
