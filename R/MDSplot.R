#############
## MDSplot ##
#############
#' @title Multidimensional scaling with MDSplot
#' @description This function computes and plots multidimensional scaling
#' analysis for dimension reduction of count expression matrix. Internally, it is
#' applied the [stats::dist()] function to the transformed count matrix to get sample-to-sample distances.
#'
#' @param exploredds object of class [DESeq2::DESeqDataSet()], generated from `exploreDDS` function.
#' @param method a `character string` indicating which correlation coefficient is to be computed,
#' based on the [stats::cor()] function. Options are: c("pearson" "kendall", "spearman").
#' @param plotly logical: when `FALSE` (default), the `ggplot2` plot will be returned.
#' `TRUE` option returns the `plotly` version of the plot.
#' @param savePlot logical: when `FALSE` (default), the plot will not be saved.
#' If `TRUE` the plot will be saved, and requires the `filePlot` argument.
#' @param filePlot file name where the plot will be saved. For more information, please consult the
#' [ggplot2::ggsave()] function.
#'
#' @return returns an object of `ggplot` or `plotly` class.
#'
#' @examples
#' ## Targets file
#' targetspath <- system.file("extdata", "targets.txt", package = "systemPipeR")
#' targets <- read.delim(targetspath, comment = "#")
#' cmp <- systemPipeR::readComp(file = targetspath, format = "matrix", delim = "-")
#' ## Count table file
#' countMatrixPath <- system.file("extdata", "countDFeByg.xls", package = "systemPipeR")
#' countMatrix <- read.delim(countMatrixPath, row.names = 1)
#' ## Plot
#' exploredds <- exploreDDS(countMatrix, targets, cmp = cmp[[1]], preFilter = NULL, transformationMethod = "rlog")
#' MDSplot(exploredds, plotly = FALSE)
#' @export
#' @importFrom DESeq2 DESeqTransform
#' @importFrom ggplot2 ggplot aes_string geom_point scale_y_reverse ggtitle ggsave
#' @importFrom plotly ggplotly
#' @importFrom stats cor dist cmdscale
#' @importFrom SummarizedExperiment assay colData
MDSplot <- function(exploredds, method = "spearman", plotly = FALSE, savePlot = FALSE, filePlot = NULL) {
    ## Add validation
    if (all(!methods::is(exploredds) == "DESeqTransform")) {
        warning("'exploredds' needs to be assignes an object of class 'DESeqTransform'.
        Here we are converting the object into a 'DESeqTransform'class for
        downstream analysis. For more information check 'help(exploreDDS)'")
        exploredds <- DESeq2::DESeqTransform(exploredds)
    }
    ## transformation to a distance matrix
    d <- stats::cor(SummarizedExperiment::assay(exploredds), method = method)
    distmat <- stats::dist(1 - d)
    ## perform MDS
    mdsData <- data.frame(stats::cmdscale(distmat))
    mds <- cbind(mdsData, as.data.frame(SummarizedExperiment::colData(exploredds)))
    Sample <- exploredds$condition
    ## plot
    plot <- ggplot2::ggplot(mds, ggplot2::aes_string("X1", "X2", color = Sample)) +
        ggplot2::geom_point(size = 3) +
        ggplot2::scale_y_reverse() +
        ggplot2::ggtitle("Multidimensional Scaling (MDS)")
    ## Save plot
    if (savePlot == TRUE) {
        ggplot2::ggsave(plot = plot, filename = filePlot)
    }
    ## Return
    if (plotly == TRUE) {
        return(plotly::ggplotly(plot))
    }
    return(plot)
}
