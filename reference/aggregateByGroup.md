# aggregateByGroup: Group-Level SNP Data Aggregation

Aggregates single-cell SNP data into group-level summaries based on a
specified metadata column. This function collapses individual cell SNP
counts into group-level matrices, which can be used for group-level
differential SNP analyses. The function supports both transplant and
non-transplant modes, donor type filtering, and normalized expression
values.

## Arguments

- group_by:

  Character. Column name in metadata to use for grouping cells. Must be
  present in cell_metadata.

- sample_column:

  Character. Column name in metadata to use for sample stratification.
  If NULL, operates in group-only mode. If specified, creates
  sample-aware aggregation that preserves sample-level information for
  statistical testing. Default: "sample_id".

- donor_type:

  Character, optional. Specific donor type to analyze (e.g., "Donor" or
  "Recipient"). If NULL, uses all cells. Ignored in non-transplant mode.

- min_cells_per_group:

  Integer. Minimum number of cells required for a group to be included
  in analysis. Groups with fewer cells are marked as
  "filtered_low_cells" in the metadata.

- use_normalized:

  Logical. Whether to include normalized depth counts in the output
  (TRUE) or only use raw counts (FALSE).

## Value

A list containing:

- ad_matrix:

  Aggregated alternative allele counts matrix (SNPs x Groups)

- dp_matrix:

  Aggregated depth matrix (SNPs x Groups)

- dp_matrix_normalized:

  Aggregated normalized depth matrix (SNPs x Groups), if available and
  requested

- metadata:

  Data frame with group-level metadata and QC metrics

- group_by:

  The metadata column used for grouping

- parameters:

  List of parameters used for aggregation

- snp_info:

  Data frame with SNP information

- snp_annotations:

  Data frame with SNP annotations

## Details

This function works by:

1.  Filtering cells based on donor_type if specified (e.g., only use
    Donor cells)

2.  Identifying unique values in the grouping column (e.g., cell_type)

3.  Summing alternative allele counts and depth counts across all cells
    in each group

4.  Creating group-level metadata with cell counts and quality metrics

5.  Filtering groups with fewer cells than the specified threshold

The function automatically detects non-transplant mode (single donor
type) and adjusts its behavior accordingly. It also checks for
normalized counts and includes them in the output if available and
requested.

## Note

- This function is typically used as a preprocessing step before
  [`findSNPsByGroup()`](https://github.com/cchmc/variantCell/reference/findSNPsByGroup.md)

- The aggregated matrices no longer contain cell-level information; all
  counts are summed across cells in each group

- For transplant data, it's often useful to analyze donor and recipient
  cells separately by specifying the donor_type parameter

- Groups with fewer cells than min_cells_per_group are marked as
  "filtered_low_cells" in the metadata but are still included in the
  output matrices

## Examples

``` r

if (FALSE) { # \dontrun{
# Basic usage - aggregate by cell type (group-only mode)
collapsed <- project$aggregateByGroup(
  group_by = "cell_type",
  sample_column = NULL,
  use_normalized = TRUE
)

# Sample-stratified aggregation for statistical testing
stratified_agg <- project$aggregateByGroup(
  group_by = "disease_status",
  sample_column = "patient_id",
  donor_type = "Recipient",
  use_normalized = TRUE
)

# Analyze only donor cells with stricter filtering
donor_agg <- project$aggregateByGroup(
  group_by = "cell_type",
  sample_column = "sample_id",
  donor_type = "Donor",
  min_cells_per_group = 5
)
} # }
```
