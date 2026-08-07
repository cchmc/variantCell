# process_tsv: Process a TSV file with prefix addition

Reads a TSV file, adds a prefix to a specified column, and sets that
column as row names in the resulting data frame. This function is
typically used to preprocess sample metadata files and ensure cell
identifiers match across different data sources.

## Arguments

- file_path:

  Character. Path to the TSV file to be read.

- prefix_column:

  Character. Name of the column to which the prefix should be added and
  which will be used as row names.

- prefix_text:

  Character. Text to prepend to each value in the specified column.

## Value

A data frame where:

- The specified column has been prefixed with the given text

- The modified column has been set as row names

- The column itself has been removed from the data frame

## Details

This function is particularly useful when processing Vireo donor
assignment files that need to be integrated with other single-cell data
(like Seurat objects), where cell identifiers need to match exactly for
proper merging.

## Note

The function assumes that the file exists and is formatted as a proper
TSV file. No validation is performed to check if the column specified by
`prefix_column` exists in the TSV file.

## Examples

``` r
if (FALSE) { # \dontrun{
# Process a donor assignment file for Vireo
donor_meta <- process_tsv(
  file_path = "path/to/donor_ids.tsv",
  prefix_column = "cell",
  prefix_text = "Patient1_Sample3_"
)

# Now donor_meta has the prefixed cell IDs as row names
# and can be easily matched with Seurat object cell barcodes
} # }
```
