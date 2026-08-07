# findSNPsByGroup: Group-Level SNP Presence Analysis

Identifies SNPs that are exclusively or predominantly present in one
cell group compared to another. This function analyzes alternative
allele frequencies between groups using aggregated data to detect
group-specific genetic variants.

## Arguments

- ident.1:

  Character. Primary group identity to analyze.

- ident.2:

  Character, optional. Secondary group identity to compare against. If
  NULL, compares against all other groups combined.

- aggregated_data:

  List. Output from aggregateByGroup function with required matrices and
  metadata.

- min_depth:

  Integer. Minimum total read depth required for a group to consider a
  SNP.

- min_alt_frac:

  Numeric between 0 and 1. Minimum alternative allele fraction required
  in a group for a SNP to be considered present.

- max_alt_frac_other:

  Numeric between 0 and 1. Maximum alternative allele fraction allowed
  in the other group for a SNP to be considered absent there.

- return_all:

  Logical. Whether to return all results regardless of significance.

- include_rs_ids:

  Logical. Whether to include rs# identifiers in results if available.

- presence_score_weights:

  Numeric vector of length 3. Weights for calculating presence scores in
  sample-stratified mode: c(fold_change, sample_consistency,
  depth_reliability). Default: c(0.4, 0.3, 0.3). Weights should sum to
  1.0.

## Value

List containing:

- results:

  Data frame of group-specific SNPs with metrics including genomic
  position, gene annotation, depth metrics, allele frequencies, and
  presence classification.

- summary:

  List with analysis overview including counts of SNPs present in each
  group and parameters used for filtering.

## Details

The function identifies SNPs that are present in one group but absent in
another by applying thresholds to alternative allele frequencies. For
each SNP, a presence score is calculated that quantifies the strength of
evidence for group-specific presence, considering both the frequency
difference and the read depth.

A SNP is considered "present" in a group when its alternative allele
frequency exceeds `min_alt_frac` and the read depth exceeds `min_depth`.
It is considered "absent" in the other group when its alternative allele
frequency is below `max_alt_frac_other` and the read depth exceeds
`min_depth`.

The presence score formula is: score = (alt_frac_present -
alt_frac_absent) \* (depth/min_depth) \* (1 -
alt_frac_absent/min_alt_frac)

## Note

- This function operates on pre-aggregated data from
  [`aggregateByGroup()`](https://github.com/cchmc/variantCell/reference/aggregateByGroup.md)
  rather than raw SNP data

- Non-transplant mode is automatically detected from the aggregated data
  parameters

- Results are sorted by presence score, with highest-scoring SNPs listed
  first

## See also

[`aggregateByGroup`](https://github.com/cchmc/variantCell/reference/aggregateByGroup.md)
for preparing input data
[`findDESNPs`](https://github.com/cchmc/variantCell/reference/findDESNPs.md)
for cell-level differential analysis
[`plotSNPs`](https://github.com/cchmc/variantCell/reference/plotSNPs.md)
for visualizing the identified SNPs

## Examples

``` r

if (FALSE) { # \dontrun{
# Aggregate SNP data by cell type
agg_data <- proj$aggregateByGroup(
  group_by = "cell_type",
  donor_type = "Donor",
  use_normalized = TRUE
)

# Find T cell-specific SNPs
tc_snps <- proj$findSNPsByGroup(
  ident.1 = "T_cells",
  ident.2 = "B_cells",
  aggregated_data = agg_data,
  min_depth = 20,
  min_alt_frac = 0.25,
  max_alt_frac_other = 0.05
)

# Comparing patient groups
patient_snps <- proj$findSNPsByGroup(
  ident.1 = "ACR",
  ident.2 = "No_ACR",
  aggregated_data = patient_data,
  min_alt_frac = 0.1,
  max_alt_frac_other = 0.02
)
} # }
```
