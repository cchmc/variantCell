# plotSNPs "Visualize SNP Distribution for a Gene"

Creates genomic visualizations of SNP distribution for a specified gene,
showing alternative allele frequencies across different cell
populations. This function generates interactive plots displaying SNPs
along the gene's coordinates with various grouping options.

## Arguments

- gene:

  Character. Name of the gene to visualize.

- group.by:

  Character, optional. Primary grouping variable from metadata. If NULL,
  all cells are treated as one group.

- split.by:

  Character, optional. Secondary grouping variable for within-group
  comparisons. Useful for donor/recipient or condition comparisons
  within cell types.

- idents:

  Vector, optional. Specific identity values to include. Must be values
  from group.by column.

- min_depth:

  Integer. Minimum read depth required to include a SNP in
  visualization.

- min_cells:

  Integer. Minimum cells per group required to include a group in
  visualization.

- min_alt_frac:

  Numeric between 0 and 1. Minimum alternative allele fraction required
  to include a SNP. Set to 0 to show all SNPs regardless of alt
  fraction.

- flank_size:

  Integer. Size of flanking regions (in bp) to include around the gene.

- plot_density:

  Logical. Whether to include density distribution plots below the main
  visualization.

- data_out:

  Logical. If TRUE, returns a data frame with SNP data instead of plots.

- use_normalized:

  Logical. Whether to use normalized read depth values.

- color_scheme:

  Vector. Named vector with "low" and "high" colors for alt fraction
  gradient.

- point_size_range:

  Numeric vector of length 2. Range of point sizes (minimum and maximum)
  for depth representation.

- include_rs_ids:

  Logical. Applies if data_out is set to true

## Value

If data_out = FALSE (default): Plot grid with main SNP visualization and
optional density plots. If data_out = TRUE: Data frame containing SNP
information for the gene including positions, alternative allele
frequencies, depths, and other metrics for each group.

## Details

This function creates a comprehensive visualization of SNPs within a
gene region, showing:

1.  A main plot with:

    - Gene structure with exons displayed as black rectangles

    - SNP positions represented as points along genomic coordinates

    - Alternative allele frequencies encoded by color

    - Read depth encoded by point size

    - Groups and splits organized vertically

2.  Optional density plots showing:

    - Alternative allele fraction distribution

    - Read depth distribution (log10 scale)

The function applies filters to ensure only high-quality data is
displayed. SNPs with low depth or low alternative allele frequency can
be excluded. Groups with insufficient cells are also not displayed.

## Note

- Requires gene annotation data from AnnotationHub (automatically
  accessed)

- May take longer to generate for genes with many SNPs or when using
  many groups

- Setting min_alt_frac = 0 shows all SNPs but may make plots very busy

- The data_out parameter is useful for further custom analysis or
  visualization

## See also

[`findDESNPs`](https://github.com/cchmc/variantCell/reference/findDESNPs.md)
for identifying differentially expressed SNPs
[`plotSNPHeatmap`](https://github.com/cchmc/variantCell/reference/plotSNPHeatmap.md)
for heatmap visualization of multiple genes

## Examples

``` r
# Basic visualization of a gene


if (FALSE) { # \dontrun{

# Split by cell type and donor type
project$plotSNPs(
  gene = "IL7R",
  group.by = "cell_type",
  split.by = "donor_type",
  min_alt_frac = 0.1
)

# Focus on specific cell types
project$plotSNPs(
  gene = "PDCD1",
  group.by = "cell_type",
  idents = c("CD4_T", "CD8_T", "Treg"),
  flank_size = 10000
)

# Return data frame instead of plot for custom visualization
snp_data <- project$plotSNPs(
  gene = "HLA-DRB1",
  group.by = "condition",
  data_out = TRUE
)
} # }
```
