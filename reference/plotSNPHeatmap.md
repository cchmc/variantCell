# plotSNPHeatmap: Create a heatmap visualization of SNP expression across cell groups

This function generates a heatmap visualization of SNP expression data
across different cell groups, with optional splitting by an additional
variable. It aggregates SNP read counts within each group and can filter
based on alternative allele frequency thresholds. The resulting heatmap
shows expression patterns with gene-based row annotation and clustering
options.

## Arguments

- genes:

  Character vector of gene names to include in the heatmap. Can be NULL
  if snp_indices is provided instead.

- snp_indices:

  Integer vector of specific SNP indices to include in the heatmap. Can
  be NULL if genes is provided instead.

- group.by:

  Character. Column name in metadata to use for primary grouping of
  cells.

- split.by:

  Character, optional. Column name in metadata to use for secondary
  grouping/splitting.

- min_alt_frac:

  Numeric, 0 to 1. Minimum alternative allele fraction required for a
  SNP to be counted as expressed in a cell. Default is 0.2.

- scale_data:

  Logical. Whether to scale data by row for visualization. Default is
  TRUE.

- max_scale:

  Numeric. Maximum value for scaled data (values will be capped at
  ±max_scale). Default is 2.

- cluster_rows:

  Logical. Whether to cluster rows (SNPs) in the heatmap. Default is
  TRUE.

- cluster_cols:

  Logical. Whether to cluster columns (cell groups) in the heatmap.
  Default is TRUE.

- show_rownames:

  Logical. Whether to display row names (SNP identifiers) in the
  heatmap. Default is TRUE.

- show_colnames:

  Logical. Whether to display column names (group identifiers) in the
  heatmap. Default is TRUE.

- fontsize_row:

  Numeric. Font size for row names. Default is 8.

- fontsize_col:

  Numeric. Font size for column names. Default is 8.

- exclude_empty:

  Logical. Whether to exclude rows and columns with no expression data.
  Default is TRUE.

- normalize_by_cells:

  Logical. Whether to normalize expression values by total cell count in
  each group (TRUE) or only by expressing cells (FALSE). Default is
  TRUE.

- data_out:

  Logical. Whether to return the underlying data matrices instead of the
  heatmap object. Default is FALSE.

- use_rs_ids:

  Logical. Whether to use rs# identifiers for row labels when available.
  Default is TRUE. Falls back to chromosome:position format if rs# not
  available.

- rs_id_format:

  Character. Format for rs# display: "rs_only" (just rs#), "rs_with_pos"
  (rs# with position), or "mixed" (rs# when available, otherwise
  position). Default is "mixed".

## Value

If data_out is FALSE (default), returns a ComplexHeatmap object that can
be directly plotted or further customized. If data_out is TRUE, returns
a list containing:

- raw_matrix: Matrix of raw expression values

- scaled_matrix: Matrix of scaled expression values

- cell_counts: Matrix of total cell counts per group

- expr_cell_counts: Matrix of expressing cell counts per group

- snp_info: Data frame with SNP identifiers, gene names, and feature
  types

## Details

This function calculates the mean expression of SNPs across cell groups,
with filtering based on minimum alternative allele frequency. For each
SNP in each group, it computes:

1.  The number of cells with the SNP

2.  The number of cells expressing the SNP above the alt_frac threshold

3.  The mean expression value (normalized by total cells or expressing
    cells)

The resulting heatmap includes annotation for gene names and feature
types, with rows grouped by gene. The heatmap uses a blue-to-red color
scale for scaled expression values.

## See also

[`aggregateByGroup`](https://github.com/cchmc/variantCell/reference/aggregateByGroup.md)
for aggregating SNP data by groups
[`findSNPsByGroup`](https://github.com/cchmc/variantCell/reference/findSNPsByGroup.md)
for identifying differential SNPs between groups
[`plotSNPs`](https://github.com/cchmc/variantCell/reference/plotSNPs.md)
for visualizing SNP distribution along genomic regions

## Examples

``` r
# Basic usage with default parameters
if (FALSE) { # \dontrun{
project$plotSNPHeatmap(genes = "BRCA1", group.by = "cell_type")

# Plot multiple genes with custom settings
project$plotSNPHeatmap(
  genes = c("TP53", "KRAS", "EGFR"),
  group.by = "cell_type",
  split.by = "patient",
  min_alt_frac = 0.1,
  cluster_rows = FALSE
)

# Return the underlying data for custom processing
snp_data <- project$plotSNPHeatmap(
  genes = "APOE",
  group.by = "condition",
  data_out = TRUE
)
} # }
```
