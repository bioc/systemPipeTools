library(SPtools)
# context("showDT")
test_that("check class", {
    x <- showDT(iris)
    expect_s3_class(x, "datatables")
})
    
# context("plots")
test_that("check output", {
targetspath <- system.file("extdata", "targets.txt", package = "systemPipeR")
targets <- read.delim(targetspath, comment = "#")[13:16,]
cmp <- systemPipeR::readComp(file = targetspath, format = "matrix", delim = "-")
cmp$CMPset1 <- cmp$CMPset1[7:7, 1:2]
## Count table file
countMatrixPath <- system.file("extdata", "countDFeByg.xls", 
                               package = "systemPipeR")
countMatrix <- read.delim(countMatrixPath, row.names = 1)[,13:16]
## exploreDDS
expect_warning(exploredds <- exploreDDS(countMatrix, targets, cmp = cmp[[1]], 
                                        preFilter = NULL, 
                                        transformationMethod = "rlog"))
## Plot exploreDDSplot
p1 <- exploreDDSplot(countMatrix, targets, cmp=cmp[[1]], preFilter=NULL, 
                     samples=c(1,2))
expect_s3_class(p1, "ggplot")
## Plot hclustplot
p2 <- hclustplot(exploredds, method = "spearman")
expect_s3_class(p2, "ggplot")
## Plot heatMaplot Samples
p3 <- heatMaplot(exploredds, clust="samples", plotly = TRUE)
expect_s3_class(p3, "plotly")

## Plot PCAplot
p4 <- PCAplot(exploredds, plotly = TRUE)
expect_s3_class(p4, "plotly")
## Plot MDSplot
p5 <- MDSplot(exploredds, plotly = FALSE)
expect_s3_class(p5, "ggplot")
## Individuals genes identified in DEG analysis
### DEG analysis with `systemPipeR`
expect_warning(degseqDF <- systemPipeR::run_DESeq2(countDF = countMatrix, 
                                                   targets = targets, 
                                    cmp = cmp[[1]], independent = FALSE))
DEG_list <- systemPipeR::filterDEGs(degDF = degseqDF, filter = c(Fold = 0, 
                                                                 FDR = 40))
## Plot heatMaplot
p6 <- heatMaplot(exploredds, clust="ind", 
                 DEGlist = unique(as.character(unlist(DEG_list[[1]]))))
expect_s3_class(p6, "pheatmap")
## Plot tSNEplot
p7 <- tSNEplot(countMatrix, targets, perplexity = 1)
expect_s3_class(p7, "ggplot")
## Plot MAplot
p8 <- MAplot(degseqDF, comparison = "M12-A12", filter = c(Fold = 1, FDR = 20), 
             genes = "ATCG00280")
expect_s3_class(p8, "ggplot")
## Plot volcanoplot
p9 <- volcanoplot(degseqDF, comparison = "M12-A12", 
                  filter = c(Fold = 1, FDR = 20), genes = "ATCG00280")
expect_s3_class(p9, "ggplot")
## Plot GLMplot
expect_warning(exploredds_raw <- exploreDDS(countMatrix, targets, cmp=cmp[[1]], 
                             preFilter=NULL, transformationMethod="raw"), 
               "some variables in design formula are characters")
p10 <- GLMplot(exploredds_raw, plotly = FALSE)
expect_s3_class(p10, "ggplot")
})