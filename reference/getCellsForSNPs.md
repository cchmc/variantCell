# getCellsForSNPs: Extract Cell IDs Based on SNP Expression Criteria

Identifies cells that express specific SNPs at defined thresholds. This
function extracts cell IDs based on alternative allele fraction and
depth criteria for specified SNPs. Useful for annotating cells in Seurat
or other single-cell objects based on genetic variant expression
patterns.

## Arguments

- snp_ids:

  Character vector. SNP identifiers to query. Three forms are accepted,
  tried in this order: the full `CHROM_POS_REF_ALT` key as it appears in
  `snp_info$snp_id`; an rs identifier, if the database was built with
  `add_rs_ids = TRUE`; or `chromosome:position` (e.g. "1:12345"), with
  or without a "chr" prefix.

- min_alt_frac:

  Numeric. Minimum alternative allele fraction threshold. Cells must
  have alt_frac \>= this value to be included. Use 0 to include cells
  not expressing the alternative allele. Default: 0.2.

- max_alt_frac:

  Numeric. Maximum alternative allele fraction threshold. Cells must
  have alt_frac \<= this value to be included. Use 1.0 to include all
  expressing cells. Set to 0 to get only reference allele cells.
  Default: 1.0.

- min_dp:

  Integer. Minimum depth (read coverage) required for the SNP in each
  cell. Cells with DP \< min_dp are excluded. Default: 5.

- sample_ids:

  Character vector, optional. Restrict analysis to specific samples. If
  NULL, uses all samples in the database. Default: NULL.

## Value

A named list where each element corresponds to a queried SNP:

- snp_id:

  Character vector of cell IDs meeting the criteria for this SNP

If a SNP is not found, returns an empty character vector with a warning.
The list also includes an attribute "summary" with per-SNP statistics.

## Details

This function searches through the SNP database to find cells expressing
specified SNPs within defined thresholds. It's particularly useful for:

- Identifying cells carrying specific mutations

- Finding cells not expressing certain alleles (min_alt_frac = 0,
  max_alt_frac = 0)

- Annotating cell subsets based on genetic variants for downstream
  analysis

- Quality control based on read depth requirements

An rs identifier or a `chromosome:position` can resolve to more than one
row, since a position may carry several alternative alleles. In that
case the alternative counts are pooled across alleles and read depth,
which is a property of the position rather than of an allele, is counted
once. Query the full `CHROM_POS_REF_ALT` key to restrict to a single
allele.

## Note

- Cells must meet ALL criteria (alt_frac range AND min_dp) to be
  included

- The function only searches SNPs present in the current database

- For cells not expressing the alt allele, use min_alt_frac=0,
  max_alt_frac=0

- Empty results may indicate SNP not found or no cells meeting criteria

## See also

[`findDESNPs`](https://github.com/cchmc/variantCell/reference/findDESNPs.md)
for differential SNP analysis
[`plotSNPs`](https://github.com/cchmc/variantCell/reference/plotSNPs.md)
for SNP visualization

## Examples

``` r
if (FALSE) { # \dontrun{
# Find cells expressing specific SNPs above 20% alt fraction
expressing_cells <- project$getCellsForSNPs(
  snp_ids = c("1:12345", "2:67890"),
  min_alt_frac = 0.2,
  min_dp = 5
)

# Find cells NOT expressing the alternative allele (reference only)
ref_cells <- project$getCellsForSNPs(
  snp_ids = "1:12345",
  min_alt_frac = 0,
  max_alt_frac = 0,
  min_dp = 5
)

# Find cells with high expression of alt allele (>80%)
high_alt_cells <- project$getCellsForSNPs(
  snp_ids = "rs123456",
  min_alt_frac = 0.8,
  min_dp = 10
)

# Use results to annotate Seurat object
seurat_obj$expressing_snp1 <- ifelse(
  colnames(seurat_obj) %in% expressing_cells$`1:12345`,
  "expressing", "not_expressing"
)
} # }
```
