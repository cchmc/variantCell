# subsetVariantCell: Subset Cells Based on Metadata Values

Filters cells in the project based on values in a specified metadata
column. This function can either create a new variantCell object with
the subset data (default) or modify the current object in-place.

## Arguments

- column:

  Character. Name of the metadata column to filter on.

- values:

  Vector. Values to include or exclude (depending on `invert`
  parameter).

- invert:

  Logical. If FALSE (default), keeps cells matching the values; if TRUE,
  keeps cells NOT matching the values.

- copy:

  Logical. If TRUE (default), returns a new variantCell object with the
  subset; if FALSE, modifies the current object in-place.

## Value

If copy=TRUE: Returns a new variantCell object containing only the
subset data. If copy=FALSE: Returns the modified object invisibly (for
method chaining).

## Details

This function subsets all project data components, including:

- Alternative allele (AD) matrix

- Depth (DP) matrix

- Normalized depth matrix (if available)

- Cell metadata

The function also handles sample management, removing samples that no
longer have any cells after filtering. This ensures data consistency
throughout the object.

With copy=TRUE (default), the original object remains untouched and a
new object with only the selected data is returned, allowing for
exploration of subsets without risk of data loss. With copy=FALSE, the
subsetting operation permanently modifies the original object, which is
more memory-efficient but irreversible.

## Note

- When copy=TRUE, this function performs a deep clone which may use
  significant memory for large datasets

- In most bioinformatics workflows, it's recommended to keep the
  original data intact and work with copies to enable different analysis
  branches

- The function prints a summary of changes for verification

## See also

[`aggregateByGroup`](https://github.com/cchmc/variantCell/reference/aggregateByGroup.md)
for grouping cells after subsetting

## Examples

``` r
if (FALSE) { # \dontrun{
# Create a new variantCell object with only T cells (default behavior)
t_cell_project <- project$subsetVariantCell("cell_type", c("CD4", "CD8", "Treg"))

# Modify the original project in-place (more memory efficient but irreversible)
project$subsetVariantCell("patient_id", "Patient1", copy = FALSE)

# Create a subset excluding cells from a specific condition
no_acr_project <- project$subsetVariantCell("condition", "ACR", invert = TRUE)

# Method chaining example (when using copy=FALSE)
results <- project$subsetVariantCell("donor_type", "Donor", copy = FALSE)$findDESNPs(
  ident.1 = "T_cells",
  ident.2 = "B_cells"
)
} # }
```
