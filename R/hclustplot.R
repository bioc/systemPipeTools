################
## hclustplot ##
################
#' @title Hierarchical Clustering Dendrogram (hclustplot)
#' @description This function computes the sample-wise correlation coefficients 
#' using the [stats::cor()] function from the transformed expression values. 
#' After transformation to a distance matrix, hierarchical clustering is 
#' performed with the [stats::hclust()] function, and the result is plotted as 
#' a dendrogram.
#'
#' @param exploredds object of class [DESeq2::DESeqTransform()].
#' @param method a `character string` indicating which correlation coefficient
#' is to be computed, based on the [stats::cor()] function. Options are:
#'    c("pearson" "kendall", "spearman").
#' @param plotly logical: when `FALSE` (default), the `ggplot2` plot will be 
#' returned. `TRUE` option returns the `plotly` version of the plot.
#' @param savePlot logical: when `FALSE` (default), the plot will not be saved.
#' If `TRUE` the plot will be saved, and requires the `filePlot` argument.
#' @param filePlot file name where the plot will be saved. For more information,
#' please consult the [ggplot2::ggsave()] function.
#'
#' @return returns an object of `ggplot` or `plotly` class.
#'
#' @examples
#' ## Targets file
#' targetspath <- system.file("extdata", "targets.txt", 
#' package = "systemPipeR")
#' targets <- read.delim(targetspath, comment = "#")
#' cmp <- systemPipeR::readComp(file = targetspath, 
#' format = "matrix", delim = "-")
#' ## Count table file
#' countMatrixPath <- system.file("extdata", "countDFeByg.xls",
#'  package = "systemPipeR")
#' countMatrix <- read.delim(countMatrixPath, row.names = 1)
#' ## Plot
#' exploredds <- exploreDDS(countMatrix, targets,
#'     cmp = cmp[[1]],
#'     preFilter = NULL, transformationMethod = "rlog"
#' )
#' hclustplot(exploredds, method = "spearman")
#' hclustplot(exploredds, method = "spearman", savePlot = TRUE, 
#' filePlot = "cor.pdf")
#' @export
#' @importFrom ape as.phylo
#' @importFrom ggplot2 coord_cartesian margin ggsave
#' @importFrom ggtree ggtree geom_tiplab theme_tree
#' @importFrom plotly ggplotly
#' @importFrom stats cor hclust dist
#' @importFrom SummarizedExperiment assay
hclustplot <- function(exploredds, method = "spearman", plotly = FALSE,
                       savePlot = FALSE, filePlot = NULL) {
    ## Validations
    if (all(!methods::is(exploredds) == "DESeqTransform")) stop("'exploredds' needs to be assignes an object of class 'DESeqTransform'. For more information check 'help(exploreDDS)'.")
    ## cor() computes the correlation coefficient
    d <- stats::cor(SummarizedExperiment::assay(exploredds), method = method)
    ## Hierarchical cluster analysis
    hc <- stats::hclust(stats::dist(1 - d))
    ## plot phylogenetic trees
    plot <- ggtree::ggtree(ape::as.phylo(hc), color = "blue") +
        ggtree::geom_tiplab() +
        ggplot2::coord_cartesian(clip = "off") +
        ggtree::theme_tree(plot.margin = ggplot2::margin(6, 60, 6, 6))
    if (savePlot == TRUE) {
        ggplot2::ggsave(filePlot, scale = 0.8)
    }
    ## Return
    if (plotly == TRUE) {
        return(plotly::ggplotly(plot))
    }
    return(plot)
}