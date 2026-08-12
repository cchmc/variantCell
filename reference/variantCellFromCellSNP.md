# variantCellFromCellSNP: Build a variantCell Object from cellsnp-lite Output

Convenience wrapper: calls
[`buildCellSNPDatabase`](https://github.com/cchmc/variantCell/reference/buildCellSNPDatabase.md)
and returns a variantCell object with the database attached and an
identity optionally set.

## Usage

``` r
variantCellFromCellSNP(..., identity = NULL)
```

## Arguments

- ...:

  Passed to
  [`buildCellSNPDatabase`](https://github.com/cchmc/variantCell/reference/buildCellSNPDatabase.md).

- identity:

  Character or NULL. Metadata column to set as the project identity.

## Value

A variantCell object with the database attached.

## See also

[`buildCellSNPDatabase`](https://github.com/cchmc/variantCell/reference/buildCellSNPDatabase.md)

## Examples

``` r
if (FALSE) { # \dontrun{
project <- variantCellFromCellSNP(
  dirs       = file.path("editing_pilot/run", samples),
  sample_ids = samples,
  cell_metadata = seurat_meta,
  identity   = "donor_type")

idx <- project$computeAlleleFractionIndex(within = "cell_type")
} # }
```
