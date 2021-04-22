###############
## tSNEplot ##
###############
#' @title t-Distributed Stochastic Neighbor embedding with tSNEplot
#' @description This function computes and plots t-Distributed Stochastic Neighbor
#' embedding (t-SNE) analysis for unsupervised nonlinear dimensionality reduction
#' of count expression matrix. Internally, it is applied the [Rtsne::Rtsne()]
#' function, using the exact t-SNE computing with `theta=0.0`.
#'
#' @param countMatrix `date.frame` or `matrix` containing raw read counts.
#' @param targets targets `data.frame`.
#' @param plotly logical: when `FALSE` (default), the `ggplot2` plot will be returned.
#' `TRUE` option returns the `plotly` version of the plot.
#' @param savePlot logical: when `FALSE` (default), the plot will not be saved.
#' If `TRUE` the plot will be saved, and requires the `filePlot` argument.
#' @param filePlot file name where the plot will be saved. For more information, please consult the
#' [ggplot2::ggsave()] function.
#' @param ... additional parameters for the [Rtsne::Rtsne()] function.
#'
#' @return returns an object of `ggplot` or `plotly` class.
#'
#' @examples
#' targetspath <- system.file("extdata", "targets.txt", package = "systemPipeR")
#' targets <- read.delim(targetspath, comment = "#")
#' cmp <- systemPipeR::readComp(file = targetspath, format = "matrix", delim = "-")
#' countMatrixPath <- system.file("extdata", "countDFeByg.xls", package = "systemPipeR")
#' countMatrix <- read.delim(countMatrixPath, row.names = 1)
#' set.seed(42)
#' tSNEplot(countMatrix, targets, perplexity = 5)
#' @export
#' @importFrom ggplot2 ggplot aes aes_string geom_point ggtitle ggsave
#' @importFrom plotly ggplotly
#' @importFrom Rtsne Rtsne
tSNEplot <- function(countMatrix, targets, plotly = FALSE, savePlot = FALSE,
                     filePlot = NULL, ...) {
    ## Validations
    if (is.data.frame(countMatrix)) {
        countMatrix <- as.matrix(countMatrix)
    } else if (is.matrix(countMatrix)) {
        countMatrix <- countMatrix
    } else {
        stop("countMatrix needs to be assigned an object of class 'data.frame' OR 'matrix'")
    }
    if (!is.data.frame(targets)) stop("targets needs to be assignes an object of
                                                                        class 'data.frame'")
    ## data manipulation
    ## removes duplicates and transpose matrix, samples perspective
    countDF_uni <- t(unique(countMatrix))
    tsne_out <- Rtsne::Rtsne(countDF_uni, dims = 2, theta = 0.0, ...)
    targets <- data.frame(targets)
    Sample <- targets$Factor
    plotdata <- data.frame(tsne_x = tsne_out$Y[, 1], tsne_y = tsne_out$Y[, 2])
    ## Plot
    plot <- ggplot2::ggplot(plotdata, ggplot2::aes_string(
        x = "tsne_x",
        y = "tsne_y"
    )) +
        ggplot2::geom_point(size = 3, ggplot2::aes(color = Sample)) +
        ggplot2::ggtitle("t-SNE")
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
