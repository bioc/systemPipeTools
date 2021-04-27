################
## heatMaplot ##
################
#' @title Hierarchical Clustering HeatMap (heatMaplot)
#' @description This function performs hierarchical clustering on the
#' transformed expression matrix generated with the DESeq2 package. It uses,
#' by default, a Pearson correlation-based distance measure and complete linkage
#' for cluster join.
#'
#' @param exploredds object of class [DESeq2::DESeqTransform()].
#' @param clust select the data to apply the distance matrix computation.
#' If `samples` selected, it will be applied the [stats::dist()] function to the
#' transformed count matrix to get sample-to-sample distances. If `ind`, it is
#' necessary to provide the list of differentially expressed genes,
#' for the `exploredds` subset.
#' @param DEGlist List of up or down regulated gene/transcript identifiers 
#' meeting the chosen filter settings for all comparisons defined in data
#' frames `pval` and `log2FC`.
#' @param plotly logical: when `FALSE` (default), the `ggplot2` plot will be
#' returned. `TRUE` option returns the `plotly` version of the plot.
#' @param savePlot logical: when `FALSE` (default), the plot will not be saved.
#' If `TRUE` the plot will be saved, and requires the `filePlot` argument.
#' @param filePlot file name where the plot will be saved. For more information,
#' please consult the [ggplot2::ggsave()] function.
#' @param ... additional parameters for the [pheatmap::pheatmap()] function.
#'
#' @return returns an object of `pheatmap` or `plotly` class.
#'
#' @examples
#' ### Load data
#' targetspath <- system.file("extdata", "targets.txt", 
#' package = "systemPipeR")
#' targets <- read.delim(targetspath, comment = "#")
#' cmp <- systemPipeR::readComp(file = targetspath, 
#' format = "matrix", delim = "-")
#' countMatrixPath <- system.file("extdata", "countDFeByg.xls", 
#' package = "systemPipeR")
#' countMatrix <- read.delim(countMatrixPath, row.names = 1)
#' ## Samples plot
#' exploredds <- exploreDDS(countMatrix, targets,
#'     cmp = cmp[[1]],
#'     preFilter = NULL, transformationMethod = "rlog"
#' )
#' heatMaplot(exploredds, clust = "samples", plotly = TRUE)
#' ## Individuals genes identified in DEG analysis
#' ### DEG analysis with `systemPipeR`
#' degseqDF <- systemPipeR::run_DESeq2(
#'     countDF = countMatrix,
#'     targets = targets, cmp = cmp[[1]], independent = FALSE
#' )
#' DEG_list <- systemPipeR::filterDEGs(
#'     degDF = degseqDF,
#'     filter = c(Fold = 2, FDR = 10)
#' )
#' ### Plot
#' heatMaplot(exploredds,
#'     clust = "ind",
#'     DEGlist = unique(as.character(unlist(DEG_list[[1]])))
#' )
#' @export
#' @importFrom ggplot2 ggsave
#' @importFrom pheatmap pheatmap
#' @importFrom plotly plot_ly
#' @importFrom stats dist
#' @importFrom SummarizedExperiment assay
heatMaplot <- function(exploredds, clust, DEGlist = NULL, plotly = FALSE,
                       savePlot = FALSE, filePlot = NULL, ...) {
    ## Validations
    if (all(!methods::is(exploredds) == "DESeqTransform")) stop("'exploredds' needs to be assignes an object of class 'DESeqTransform'. For more information check 'help(exploreDDS)'")
    anno <- as.data.frame(exploredds$condition)
    colnames(anno) <- "Condition"
    ## sample-to-sample distances
    if (clust == "samples") {
        sampleDists <- stats::dist(t(SummarizedExperiment::assay(exploredds)))
        sampleDistMatrix <- as.matrix(sampleDists)
        rownames(anno) <- colnames(sampleDistMatrix)
        if (plotly == FALSE) {
            pheatPlot <- pheatmap::pheatmap(sampleDistMatrix,
                                            clustering_distance_rows = sampleDists,
                                            clustering_distance_cols = sampleDists, annotation_col = anno
            )
        } else if (plotly == TRUE) {
            plot <- plotly::plot_ly(
                x = colnames(sampleDistMatrix), y = rownames(sampleDistMatrix),
                z = sampleDistMatrix, type = "heatmap"
            )
        }
    } else if (clust == "ind") {
        ## Hierarchical clustering on the transformed expression matrix 
        ## subsetted by the DEGs identified in differential expression analysis.
        if (any(is.null(DEGlist) | !is.character(DEGlist))) stop("Provide a character vector with the gene names identified in differential expression analysis.")
        dist <- SummarizedExperiment::assay(exploredds)[DEGlist, ]
        rownames(anno) <- colnames(dist)
        if (plotly == FALSE) {
            pheatPlot <- pheatmap::pheatmap(dist,
                                            scale = "row", clustering_distance_rows = "correlation",
                                            clustering_distance_cols = "correlation", annotation_col = anno
            )
        } else if (plotly == TRUE) {
            plot <- plotly::plot_ly(
                x = colnames(dist), y = rownames(dist), z = dist,
                type = "heatmap"
            )
        }
    } else {
        stop("Supported clust include 'samples' and 'ind'")
    }
    if (savePlot == TRUE) {
        ggplot2::ggsave(plot = pheatPlot, filename = filePlot)
    }
    ## Return
    if (plotly == TRUE) {
        return(plot)
    }
    return(pheatPlot)
}