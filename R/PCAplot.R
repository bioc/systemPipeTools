#############
## PCAplot ##
#############
#' @title PCAplot
#' @description This function plots a Principal Component Analysis (PCA) from
#' transformed expression matrix. This plot shows samples variation based on the
#' expression values and identifies batch effects.
#'
#' @param exploredds object of class [DESeq2::DESeqTransform()].
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
#' targetspath <- system.file("extdata", "targets.txt", package = "systemPipeR")
#' targets <- read.delim(targetspath, comment = "#")
#' cmp <- systemPipeR::readComp(file = targetspath, format = "matrix", 
#' delim = "-")
#' ## Count table file
#' countMatrixPath <- system.file("extdata", "countDFeByg.xls", 
#' package = "systemPipeR")
#' countMatrix <- read.delim(countMatrixPath, row.names = 1)
#' ## Plot
#' exploredds <- exploreDDS(countMatrix, targets, cmp = cmp[[1]], 
#' preFilter = NULL, transformationMethod = "rlog")
#' PCAplot(exploredds, plotly = TRUE)
#' @export
#' @importFrom DESeq2 DESeqTransform plotPCA
#' @importFrom ggplot2 ggplot aes_string geom_point xlab ylab coord_fixed 
#' ggtitle ggsave
#' @importFrom plotly ggplotly
PCAplot <- function(exploredds, plotly = FALSE, savePlot = FALSE,
                    filePlot = NULL) {
    ## Validations
    if (!inherits(exploredds, "DESeqTransform")) {
        stop("'exploredds' needs to be assignes an object of class 
                'DESeqTransform'. For more information check 'help(exploreDDS)'"
             )
    }
    ## Plot
    pcaData <- DESeq2::plotPCA(exploredds,
        intgroup = "condition",
        returnData = TRUE
    )
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    Sample <- exploredds$condition
    plot <- ggplot2::ggplot(pcaData, ggplot2::aes_string("PC1", "PC2",
        color = Sample
    )) +
        ggplot2::geom_point(size = 3) +
        ggplot2::xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ggplot2::ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        ggplot2::coord_fixed() +
        ggplot2::ggtitle("Principal Component Analysis (PCA)")
    ## Save plot
    if (savePlot == TRUE) {
        if (is.null(filePlot)) {
            stop("Argument 'filePlot' is missing, please provide file name.")
        }
        ggplot2::ggsave(plot = plot, filename = filePlot)
    }
    ## Return
    if (plotly == TRUE) {
        return(plotly::ggplotly(plot))
    }
    return(plot)
}
