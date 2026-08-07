# getCurrentIdentity: Get Current Project-Wide Cell Identity Information

Retrieves and displays information about the currently active cell
identity variable. This function returns details about which metadata
column is being used for cell grouping and provides a summary of the
cell distribution across the different identity categories.

## Value

Invisibly returns a list containing:

- identity:

  Character. Name of the current identity column.

- distribution:

  Table. Distribution of cells across identity categories.

- total_cells:

  Integer. Total number of cells in the dataset.

If no identity is set, returns NULL.

## Details

This function checks whether a project-wide identity has been set using
[`setProjectIdentity()`](https://github.com/cchmc/variantCell/reference/setProjectIdentity.md).
If an identity is active, it prints the name of the identity column and
displays a summary of how many cells belong to each identity category.

The function is useful for:

- Verifying which grouping variable is currently active

- Checking the cell distribution across groups before analysis

- Confirming that identity assignments are as expected

- Debugging when analysis results seem unexpected

## Note

- If no identity has been set, the function returns NULL and displays a
  notification

- The function both prints information to the console and returns data
  invisibly

- The returned list can be captured and used programmatically if needed

## See also

[`setProjectIdentity`](https://github.com/cchmc/variantCell/reference/setProjectIdentity.md)
for setting the active identity

## Examples

``` r

if (FALSE) { # \dontrun{

# Check current identity
project$getCurrentIdentity()

# Capture the return value for programmatic use
id_info <- project$getCurrentIdentity()
if(!is.null(id_info)) {
  # Find the most abundant cell type
  most_common <- names(which.max(id_info$distribution))
  cat("Most common cell type:", most_common)
}

# Use in a workflow
project$setProjectIdentity("cell_type")
project$getCurrentIdentity()  # Verify it worked
project$findDESNPs(ident.1 = "T_cells")
} # }
```
