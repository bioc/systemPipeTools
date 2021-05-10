#############
## GLMplot ##
#############
#' @title Dimension Reduction with GLMplot
#' @description This function computes and plots generalized principal 
#' components analysis for dimension reduction of count expression matrix.
#'
#' @param exploredds object of class [DESeq2::DESeqDataSet()], generated from 
#' `exploreDDS` function.
#' @param L desired number of latent dimensions (positive integer).
#' @param plotly logical: when `FALSE` (default), the `ggplot2` plot will be 
#' returned. `TRUE` option returns the `plotly` version of the plot.
#' @param savePlot logical: when `FALSE` (default), the plot will not be saved.
#' If `TRUE` the plot will be saved, and requires the `filePlot` argument.
#' @param filePlot file name where the plot will be saved. For more information,
#' please consult the [ggplot2::ggsave()] function.
#' @param ... additional parameters for the [glmpca::glmpca()] function.
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
#' preFilter = NULL, transformationMethod = "raw")
#' GLMplot(exploredds, plotly = FALSE)
#' @export
#' @importFrom ggplot2 ggplot aes_string geom_point aes coord_fixed ggtitle 
#' ggsave
#' @importFrom glmpca glmpca
#' @importFrom plotly ggplotly
#' @importFrom SummarizedExperiment assay
#' @keywords visualization
#' @references 
#'   F. William Townes and Kelly Street (2020). glmpca: Dimension Reduction of 
#'   Non-Normally Distributed Data. R package version 0.2.0. 
#'   <https://CRAN.R-project.org/package=glmpca>
GLMplot <- function(exploredds, L = 2, plotly = FALSE, savePlot = FALSE,
                    filePlot = NULL, ...) {
    ## Add validation, need to be counts reads
    if (inherits(exploredds, "DESeqDataSet")) {
        count_mat <- SummarizedExperiment::assay(exploredds)
    } else if (!inherits(exploredds, "DESeqDataSet")) {
          stop("'exploredds' needs to be assignes an object of class 
             'DESeqDataSet'. For more information check 
             'help(exploreDDS)', and select the transformationMethod='raw'")
      }
    ## glmpca is performed on raw counts
    nozero <- count_mat[which(rowSums(count_mat) > 0), ]
    gpca <- glmpca::glmpca(nozero, L = L, ...)
    gpca.dat <- gpca$factors
    gpca.dat$condition <- exploredds$condition
    Samples <- as.character(exploredds$condition)
    plot <- ggplot2::ggplot(gpca.dat, ggplot2::aes_string("dim1", "dim2")) +
        ggplot2::geom_point(size = 3, ggplot2::aes(color = Samples)) +
        ggplot2::coord_fixed() +
        ggplot2::ggtitle("Generalized PCA (GLM-PCA)")
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
