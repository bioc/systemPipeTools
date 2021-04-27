#############################
## showDT function ##
#############################
#' @title Create an HTML table using DT package with fixed columns
#' @description Create an HTML table using DT package with fixed columns
#'
#' @param data data object (either a matrix or a data frame).
#' @param ... Additional arguments used by dDT::atatable() function.
#' 
#' @return returns an object of `datatables` and `htmlwidget`.
#' @importFrom DT datatable
#' @export
#'
#' @examples
#' showDT(iris)
showDT <- function(data, ...) {
    DT::datatable(
        data,
        extensions = c("FixedColumns", "Scroller"),
        options = list(
            scrollX = TRUE,
            fixedColumns = TRUE,
            deferRender = TRUE,
            scrollY = 200,
            scroller = TRUE
        )
    )
}
## Usage:
# targetspath <- system.file("extdata", "targets.txt", package="systemPipeR") 
# targets <- read.delim(targetspath, comment.char = "#")
# showDT(targets)