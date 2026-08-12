# computeAlleleFractionIndex: Aggregate Alt-Allele Fraction per Pseudobulk

Collapses AD/DP over a set of sites into a single alt-fraction index per
pseudobulk, where a pseudobulk is one level of `split_by` (normally
sample) crossed with one level of `group_by` and optionally `within`.

This is the Alu Editing Index construction generalised: one number per
pseudobulk rather than one test per site. Because it produces a single
statistic it carries no genome-wide multiple-testing burden, which is
what makes it usable at the group sizes typical of a transplant cohort,
where a per-site scan cannot reach significance at all.

## Arguments

- group_by:

  Character. Metadata column defining the arms being compared. Defaults
  to the current project identity.

- split_by:

  Character. Metadata column defining the unit of observation. Defaults
  to "sample_id". Do not set this to a per-cell column - the whole point
  is that cells within a library are not independent observations.

- within:

  Character or NULL. Optional further stratification, e.g. a cell type
  column, producing one index per sample x arm x cell type.

- sites:

  Optional site restriction: row indices, rs IDs, or "CHROM:POS"
  strings. NULL uses all sites in the database.

- min_dp_site:

  Numeric. Minimum pseudobulk depth for a site to contribute.

- min_cells:

  Numeric. Minimum cells for a pseudobulk to be reported.

- common_sites:

  Logical. If TRUE (default), restricts every pseudobulk within a
  `split_by` unit to sites covered in all of that unit's arms, so an
  index difference cannot be manufactured by the two arms using
  different sets of positions.

- weighted:

  Logical. If TRUE (default) the index is the read-weighted
  sum(AD)/sum(DP). If FALSE it is the unweighted mean of per-site
  fractions, which is dominated by low-coverage sites.

- verbose:

  Logical. Print progress.

## Value

A data frame with one row per pseudobulk: split, group, within (if
requested), n_cells, n_sites, total_ad, total_dp, index. Attributes
carry the parameters used and `min_attainable_p_paired`, the permutation
ceiling for a paired test over the observed number of units.

## Why the sample is the unit

Cells within a library share a genome, a capture and a sequencing run.
Treating thousands of them as independent observations is
pseudoreplication and inflates significance by orders of magnitude.
Cells are aggregated to pseudobulk first; any test then runs across
pseudobulks.

## See also

[`findDEAlleleFraction`](https://github.com/cchmc/variantCell/reference/findDEAlleleFraction.md)
for the per-site test,
[`findSNPsByGroup`](https://github.com/cchmc/variantCell/reference/findSNPsByGroup.md)
for genotype,
[`findDESNPs`](https://github.com/cchmc/variantCell/reference/findDESNPs.md)
for depth

## Examples

``` r
if (FALSE) { # \dontrun{
project$setProjectIdentity("donor_type")

# One index per sample per arm
idx <- project$computeAlleleFractionIndex()

# Stratified by cell type - the paired within-sample contrast
idx <- project$computeAlleleFractionIndex(within = "cell_type")

# Restricted to a candidate site list
idx <- project$computeAlleleFractionIndex(sites = editing_sites)
} # }
```
