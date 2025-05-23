#' variantCell: Single-Cell SNP Analysis for Transplant Studies
#' variantCell: A tool for analyzing single-cell SNP data with focus on organ transplant
#'
#' @description
#' A tool for analyzing single-cell SNP data with focus on organ transplant.
#' variantCell processes CellSNP and Vireo output to identify, characterize, and
#' visualize SNPs at cellular resolution, enabling researchers to distinguish donor from
#' recipient cells and quantify SNPs that may influence transplant outcomes.
#'
#' @section Core Features:
#' \itemize{
#'   \item Processing CellSNP and Vireo output data
#'   \item SNP identification and characterization across cell types
#'   \item Differential SNP analysis
#'   \item Visualization of SNP distributions and patterns
#'   \item Integration with single-cell metadata
#' }
#'
#' @import R6
#' @import ggplot2
#' @import cowplot
#' @import GenomicRanges
#' @import AnnotationHub
#' @import ensembldb
#' @import circlize
#' @import ComplexHeatmap
#'
#' @importFrom IRanges IRanges
#' @importFrom data.table fread rbindlist setkey as.data.table data.table
#' @importFrom Matrix Matrix rowSums colSums Diagonal sparseMatrix t nnzero
#' @importFrom matrixStats rowVars rowMeans2 rowSds colVars colMeans2
#' @importFrom grid unit gpar
#' @importFrom stats p.adjust wilcox.test t.test fisher.test chisq.test median
#' @importFrom utils head tail read.delim write.csv
#' @importFrom grDevices colorRampPalette
#' @importFrom methods is as new setGeneric setMethod
#'
#' @keywords internal
"_PACKAGE"
