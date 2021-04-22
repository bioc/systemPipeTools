# systemPipeTools: Tools for data visualization

<!-- badges: start -->
![R-CMD-check](https://github.com/systemPipeR/systemPipeTools/workflows/R-CMD-check/badge.svg)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![codecov](https://codecov.io/gh/systemPipeR/systemPipeTools/branch/master/graph/badge.svg?token=PwWVN6tTh3)](https://codecov.io/gh/systemPipeR/systemPipeTools)
<!-- badges: end -->

The systemPipeTools package provides a set of utilities for High Throughput Sequence Data
Visualization. This package is designed to extend and provide a visualization
and utilities interface for the `systemPipeR` package. 

### Installation

Get the released version from Bioconductor:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager") }
BiocManager::install("systemPipeTools"")
```
Or the development version from GitHub:
```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager") }
BiocManager::install("systemPipeR/systemPipeTools"")
```