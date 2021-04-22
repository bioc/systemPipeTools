####################
## exploreDDSplot ##
####################
#' @title exploreDDSplot
#' @description Scatterplot of transformed counts from 
#' two samples or grid of all samples
#'
#' @param countMatrix `date.frame` or `matrix` containing raw read counts
#' @param targets targets `data.frame`
#' @param cmp `character matrix` where comparisons are defined in two columns.
#' This matrix should be generated with the [systemPipeR::readComp()] function
#' from the targets file. Values used for comparisons need to match those in the
#'    `Factor` column of the targets file.
#' @param preFilter allows removing rows in which there are very few reads.
#' Accepts a numeric value with the minimum of total reads to keep.
#' Default is `NULL`.
#' @param samples a `character vector` of two samples or `ALL` samples in the
#' dataset. Could be specified the `SampleName` column name of the targets
#' file or the respective numeric values. Also, if set as `ALL`,
#' a correlation matrix it will be plot.
#' @param blind logical, whether to blind the transformation to the experimental
#' design (see varianceStabilizingTransformation), from [DESeq2::vst()]
#' or [DESeq2::rlog()].
#' @param scattermatrix if `samples` set as `ALL`, requires to assign `TRUE`
#' to build    a correlation matrix and plot the correlogram of all the samples.
#' @param plotly logical: when `FALSE` (default), the `ggplot2` plot will be
#' returned. `TRUE` returns the `plotly` version of the plot.
#' @param savePlot logical: when `FALSE` (default), the plot will not be saved.
#' If `TRUE` the plot will be saved, and requires the `filePlot` argument.
#' @param filePlot file name where the plot will be saved. For more information,
#'    please consult the [ggplot2::ggsave()] function.
#'
#' @return returns an object of `ggplot2 plot`.

#' @examples
#' library(systemPipeR)
#' ## Targets file
#' targetspath <- system.file("extdata", "targets.txt", package = "systemPipeR")
#' targets <- read.delim(targetspath, comment = "#")
#' cmp <- systemPipeR::readComp(
#'     file = targetspath,
#'     format = "matrix", delim = "-"
#' )
#' ## Count table file
#' countMatrixPath <- system.file("extdata", "countDFeByg.xls",
#'     package = "systemPipeR"
#' )
#' countMatrix <- read.delim(countMatrixPath, row.names = 1)
#' ## Plot
#' exploreDDSplot(countMatrix, targets,
#'     cmp = cmp[[1]], preFilter = NULL,
#'     samples = c(3, 4)
#' )
#' exploreDDSplot(countMatrix, targets,
#'     cmp = cmp[[1]], samples = c("M1A", "M1B"), save = TRUE,
#'     filePlot = "transf_deseq2.pdf"
#' )
#' ## Plot Correlogram
#' exploreDDSplot(countMatrix, targets,
#'     cmp = cmp[[1]], preFilter = NULL,
#'     samples = c("M1A", "M1B"), scattermatrix = TRUE
#' )
#' @export exploreDDSplot
#' @importFrom DESeq2 estimateSizeFactors counts
#' @importFrom dplyr bind_rows as_tibble mutate group_by do
#' @importFrom GGally ggpairs
#' @importFrom ggplot2 aes aes_string ggplot geom_hex coord_fixed facet_grid
#' @importFrom ggplot2 .data xlab ylab ggsave
#' @importFrom plotly plot_ly subplot
#' @importFrom SummarizedExperiment assay
#' @importFrom magrittr "%>%"
exploreDDSplot <- function(countMatrix, targets, cmp = cmp[[1]],
                           preFilter = NULL, samples, blind = TRUE,
                           scattermatrix = FALSE, plotly = FALSE,
                           savePlot = FALSE, filePlot = NULL) {
    ## Validations
    SampleName <- targets$SampleName
    names(SampleName) <- targets$SampleName
    if (is.numeric(samples)) {
        samples <- SampleName[samples]
        if (!all(samples %in% SampleName)) {
            stop(paste("samples position can be
                       assigned from the following options",
                       paste0(1:length(SampleName),
                              collapse = ", "
                       ),
                       sep = " "
            ))
        }
    } else if (is.character(samples)) {
        if (all(samples == "ALL")) {
            samples <- SampleName
            if (!scattermatrix == "TRUE") stop("'scattermatrix' argument needs
                                               to set as TRUE in the case of ALL 
                                               the samples selected.")
        } else {
            samples <- SampleName[samples]
            if (!all(samples %in% SampleName)) {
                stop(paste("samples names can be
     assigned from the following
     options",
     paste0((SampleName),
            collapse = ", "
     ),
     sep = " "
                ))
            }
        }
    }
    transformation <- . <- NULL
    ## Calculate the data transformations
    suppressWarnings({
        vst <- exploreDDS(countMatrix, targets,
                          cmp = cmp, preFilter = preFilter,
                          transformationMethod = "vst", blind = blind
        )
        rlog <- exploreDDS(countMatrix, targets,
                           cmp = cmp, preFilter = preFilter,
                           transformationMethod = "rlog", blind = blind
        )
        dss <- exploreDDS(countMatrix, targets,
                          cmp = cmp, preFilter = preFilter,
                          transformationMethod = "raw"
        )
        dss <- DESeq2::estimateSizeFactors(dss)
    })
    ## create dataframe with transformed values
    transform_df <- dplyr::bind_rows(
        dplyr::as_tibble(log2(DESeq2::counts(dss, normalized = TRUE)[, samples] + 1)) %>%
            dplyr::mutate(transformation = "log2(x + 1)"),
        dplyr::as_tibble(SummarizedExperiment::assay(vst)[, samples]) %>% dplyr::mutate(transformation = "vst"),
        dplyr::as_tibble(SummarizedExperiment::assay(rlog)[, samples]) %>% dplyr::mutate(transformation = "rlog")
    )
    names <- colnames(transform_df)[1:2]
    lvl <- levels(factor(transform_df$transformation))
    ## plot
    if (scattermatrix == TRUE) {
        plot <- GGally::ggpairs(transform_df,
                                title = "Scatterplot of transformed counts",
                                ggplot2::aes_string(colour = "transformation")
        )
    } else {
        plot <- ggplot2::ggplot(transform_df, ggplot2::aes(
            x = .data[[names[1]]],
            y = .data[[names[2]]]
        )) +
            ggplot2::geom_hex(bins = 80) +
            ggplot2::coord_fixed() +
            ggplot2::facet_grid(. ~ transformation) +
            ggplot2::xlab(names[1]) +
            ggplot2::ylab(names[2])
    }
    if (savePlot == TRUE) {
        ggplot2::ggsave(filePlot, scale = 0.8)
    }
    ## Return
    if (plotly == TRUE) {
        plot <- transform_df %>%
            dplyr::group_by(transformation) %>%
            dplyr::do(p = plotly::plot_ly(.,
                                          x = .data[[names[1]]], y = .data[[names[2]]], 
                                          color = ~transformation, type = "scatter",
                                          name = ~transformation, showlegend = TRUE, 
                                          legendgroup = ~transformation
            )) %>%
            plotly::subplot(nrows = 1, shareX = TRUE, shareY = TRUE)
    }
    return(plot)
}
