# downsampleVariant: Downsample Cells by Group to Balance Cell Numbers

Reduces the number of cells in the variantCell object by downsampling
each group to a maximum number of cells. This function is useful for
balancing cell numbers across groups, reducing computational burden, and
mitigating the effects of groups with very different cell counts on
downstream analyses.

## Arguments

- max_cells:

  Integer. Maximum number of cells to keep from each group. Groups with
  fewer cells than this threshold will retain all their cells. Default:
  1000.

- group_by:

  Character, optional. Metadata column to use for grouping cells. If
  NULL, uses the current project identity set by setProjectIdentity().
  Default: NULL.

- seed:

  Integer. Random seed for reproducible downsampling. Default: 42.

## Value

Returns the modified object invisibly (for method chaining).

## Details

The function performs downsampling by:

1.  Grouping cells based on the specified metadata column

2.  For each group, if cell count exceeds max_cells, randomly selecting
    max_cells cells to keep

3.  Updating all matrices and metadata to include only the selected
    cells

4.  Maintaining consistency across all data structures in the object

This operation modifies the object in-place, permanently removing cells
that aren't selected. It's particularly useful when working with
imbalanced datasets, where some cell types or conditions have many more
cells than others, which could bias analytical results.

The function automatically handles updates to all relevant data
structures, including:

- Alternative allele (AD) matrix

- Depth (DP) matrix

- Normalized depth matrix (if available)

- Cell metadata

- Sample-level information

## Note

- This function modifies the object in-place (no copy is created)

- Downsampling is performed randomly for each group

- The seed parameter ensures reproducibility of random sampling

- Groups with fewer cells than max_cells will keep all their cells

- If after downsampling a sample has no remaining cells, it will be
  removed from the object

- A detailed summary of the downsampling is printed to the console

## See also

[`setProjectIdentity`](https://github.com/cchmc/variantCell/reference/setProjectIdentity.md)
for setting the grouping identity
[`subsetVariantCell`](https://github.com/cchmc/variantCell/reference/subsetVariantCell.md)
for other filtering operations

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage - downsample to 500 cells per cell type
project$setProjectIdentity("cell_type")
project$downsampleVariant(max_cells = 500)

# Downsample by a different grouping variable
project$downsampleVariant(
  max_cells = 200,
  group_by = "condition",
  seed = 123  # Use different seed for different random selection
)

# Use with method chaining
results <- project$downsampleVariant(max_cells = 300)$findDESNPs(
  ident.1 = "T_cells",
  ident.2 = "B_cells"
)
} # }
```
