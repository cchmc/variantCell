# buildCellSNPDatabase: Import cellsnp-lite Output as a SNP Database

Reads one or more cellsnp-lite output directories and merges them into
the `snp_database` structure used by variantCell, taking the union of
sites across samples.

This is the import path for runs against a non-germline region list -
RNA editing sites, chrM, a candidate panel - which have no vireo output
because there is no genotype to deconvolve. For the germline genotyping
workflow use
[`addSampleData()`](https://github.com/cchmc/variantCell/reference/addSampleData.md)
plus
[`buildSNPDatabase()`](https://github.com/cchmc/variantCell/reference/buildSNPDatabase.md).

## Usage

``` r
buildCellSNPDatabase(dirs, sample_ids, cell_metadata = NULL,
  barcode_prefix = NULL, min_total_dp = 0, keep_oth = TRUE,
  require_metadata = FALSE, verbose = TRUE)
```

## Arguments

- dirs:

  Character vector of cellsnp-lite output directories.

- sample_ids:

  Character vector of sample names, same length as `dirs`.

- cell_metadata:

  Optional data frame of per-cell metadata to join. Must contain a
  `cell_id` column, or have cell IDs as rownames, matching the
  constructed cell IDs. Columns such as `donor_type` and `cell_type`
  come from here.

- barcode_prefix:

  Character vector, or NULL. Prefix prepended to each sample's barcodes
  to build cell IDs. Defaults to `paste0(sample_ids, "_")`, and must
  reproduce the cell naming used by whatever object supplies
  `cell_metadata`.

- min_total_dp:

  Numeric. Drop sites whose summed depth across all cells and samples
  falls below this. 0 keeps everything.

- keep_oth:

  Logical. Retain the OTH matrix. Default TRUE.

- require_metadata:

  Logical. If TRUE, cells absent from `cell_metadata` are dropped. If
  FALSE (default) they are kept with NA metadata.

- verbose:

  Logical. Print progress.

## Value

A list with `ad_matrix`, `dp_matrix`, `oth_matrix` (or NULL),
`dp_matrix_normalized` (NULL), `snp_info`, `cell_metadata` and
`qc_report`. The QC report carries the pooled alt fraction, the per-base
error floor and their ratio.

## Why OTH is retained

At an A\>G site the two remaining bases can only be sequencing error, so
OTH/DP is an internal, position-matched error floor for the AD/DP
signal. Comparing the two is what separates a real editing signal from a
base-calling artifact. On the 2026-08-12 editing pilot the ratio was
7.4x (TBX3) and 9.5x (CLAD-2).

## See also

[`variantCellFromCellSNP`](https://github.com/cchmc/variantCell/reference/variantCellFromCellSNP.md),
[`computeAlleleFractionIndex`](https://github.com/cchmc/variantCell/reference/computeAlleleFractionIndex.md)

## Examples

``` r
if (FALSE) { # \dontrun{
db <- buildCellSNPDatabase(
  dirs       = file.path("editing_pilot/run", c("TBX3", "CLAD-2")),
  sample_ids = c("TBX3", "CLAD-2"),
  cell_metadata = seurat_meta)

db$qc_report$signal_to_error
} # }
```
