# setProjectIdentity: Set Project-Wide Cell Identity Variable

Sets the metadata column to use as the primary cell identity variable
for all downstream analyses. This function establishes which cell
grouping will be used in functions like findDESNPs() and other SNP
analysis methods.

## Arguments

- ident_col:

  Character. Name of the column in cell_metadata to use as cell
  identity. Must be an existing column in the metadata.

## Value

Returns the object invisibly (for method chaining).

## Details

This function validates that the specified column exists in the cell
metadata, then sets it as the active identity for all subsequent
analyses. It also prints a summary of the unique identities found and
the number of cells in each category.

The identity column is crucial for many analysis functions as it defines
how cells are grouped when comparing SNP expression or presence between
populations.

Common identity variables include cell type annotations, cluster IDs,
condition labels, or any other categorical metadata that meaningfully
separates cell populations.

## Note

- This function must be called before using analysis methods that rely
  on cell identities

- The function prints a summary of the identities found, which can be
  useful for verification

- The previously set identity (if any) is overwritten by this function

## See also

[`getCurrentIdentity`](https://github.com/cchmc/variantCell/reference/getCurrentIdentity.md)
for checking the currently active identity
[`findDESNPs`](https://github.com/cchmc/variantCell/reference/findDESNPs.md)
and
[`findSNPsByGroup`](https://github.com/cchmc/variantCell/reference/findSNPsByGroup.md)
which use the set identity

## Examples

``` r

if (FALSE) { # \dontrun{

# Set cell type as the active identity

project$setProjectIdentity("cell_type")

# Use cluster IDs instead
project$setProjectIdentity("seurat_clusters")

# Set disease status as identity for case-control comparisons
project$setProjectIdentity("disease_status")

# Method chaining example
results <- project$setProjectIdentity("cell_type")$findDESNPs(
  ident.1 = "T_cells",
  donor_type = "Donor"
)
} # }
```
