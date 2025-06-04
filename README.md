# variantCell 0.1.0 - alpha (in development)

[![Documentation](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://cchmc.github.io/variantCell/) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![Status: In Development](https://img.shields.io/badge/Status-In%20Development-blue)]()

A tool for analyzing single-cell SNP data with focus on organ transplant.

## Overview

variantCell is an R-based computational tool designed to analyze single-cell SNP data with a particular focus on transplant settings. It works with widely-used single-cell genomics workflows to identify, characterize, and visualize SNPs at cellular resolution, enabling researchers to distinguish donor from recipient cells and quantify SNPs that may influence transplant outcomes.

The package provides a comprehensive framework for processing CellSNP and Vireo output, normalizing SNP counts, and performing differential SNP analysis across cell types or clinical conditions. With built-in visualization capabilities, users can explore variant distributions within genes of interest and create SNP heatmaps across cell clusters.

Unlike existing approaches that analyze variants at the bulk tissue level, variantCell leverages single-cell resolution to reveal cell type-specific genetic signals and heterogeneity within transplanted tissues. This approach offers unique insights into the genetic factors driving transplant rejection, immune response variation, and cellular chimerism in complex tissues, with applications extending to cancer, development, and other areas where genetic mosaicism is biologically significant.

## Core Dependencies

R6, data.table, Matrix, ggplot2, cowplot, GenomicRanges, IRanges, AnnotationHub, matrixStats, ensembldb, circlize, ComplexHeatmap

## Optional Dependencies

Seurat: If using Seurat objects as input (not required if using data frames), SingleCellExperiment: If using SingleCellExperiment objects as input Parallel processing findDESNPs function: parallel, doParallel, foreach

## Getting Started

**Documentation:** <https://cchmc.github.io/variantCell/>

To add samples to the SNP database, the package requires: - An output directory from CellSNP along with cell metadata (as a Seurat object, dataframe, or SingleCellExperiment) - Cell prefixes can be optionally specified when adding sample data to match cell IDs to integrated data

For transplant mode, you'll also need: - A Vireo output directory - Specification of which donors are 'Recipient' and 'Donor' (can be determined using the `process_vireo` helper function)

After adding samples, a SNP database can be built using the `buildSNPDatabase` function, enabling:

-   Finding SNPs present or absent in specific groups

-   Performing cell (or group)-level differential SNP expression analysis (uses parallel processing)

SNPs can be visualized:

-   Plotting SNPs for specific genes

-   Creating heatmaps of SNPs within cell clusters / groups

## Installation

Install from GitHub

Check if devtools is installed, install if needed

`if (!requireNamespace("devtools", quietly = TRUE)) { install.packages("devtools") }`

`devtools::install_github("cchmc/variantCell")`

## License

This project is licensed under the MIT License.

## Acknowledgements

Funding: Heart Institute, Cincinnati Children's Hospital Medical Center

## Citation

If you use this software, please cite:

Andrew Potter, Don Hayes (in preparation). variantCell: A tool for single-cell SNP analysis. GitHub: <https://cchmc.github.io/variantCell/>
