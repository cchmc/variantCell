# findDEAlleleFraction: Per-Site Alt-Allele Fraction Between Groups

Tests, at each site, whether the fraction of reads carrying the
alternative base differs between two groups. This is the statistic for
RNA-level variation - editing rate, allele-specific expression, mtDNA
heteroplasmy - as distinct from
[`findSNPsByGroup()`](https://github.com/cchmc/variantCell/reference/findSNPsByGroup.md),
which tests genotype, and
[`findDESNPs()`](https://github.com/cchmc/variantCell/reference/findDESNPs.md),
which tests read depth. Nothing else in the package tests AD/DP itself.

Cells are aggregated to one pseudobulk per `split_by` unit per arm
before testing, so the unit of observation is the sample rather than the
cell.

## Arguments

- ident.1:

  Character. First group.

- ident.2:

  Character or NULL. Second group; NULL means all other cells.

- group_by:

  Character. Identity column. Defaults to the current identity.

- split_by:

  Character. Unit of observation, default "sample_id".

- sites:

  Optional site restriction: row indices, rs IDs, or "CHROM:POS".

- min_dp:

  Numeric. Minimum pseudobulk depth per arm per unit.

- min_units:

  Numeric. Minimum units per arm for a site to be tested.

- paired:

  Logical. If TRUE, units contributing both arms are tested with a
  paired signed-rank test on the per-unit fraction difference. Falls
  back to unpaired with a warning if fewer than two units contribute
  both arms.

- test:

  Character. "wilcoxon" (default, distribution-free) or "betabinom"
  (pooled likelihood ratio with method-of-moments overdispersion).

- p_adjust:

  Character. Method passed to `p.adjust`.

- verbose:

  Logical. Print progress.

## Value

A data frame of per-site results ordered by p-value: CHROM, POS, rs_id,
gene_name, af_1, af_2, af_diff, dp_1, dp_2, n_units, p_value, p_adj.
Attributes carry `min_attainable_p`, `bonferroni_threshold` and
`significance_reachable`.

## Statistical ceiling

The function reports the smallest p-value attainable given the number of
units and warns when that ceiling cannot clear Bonferroni over the
number of sites tested. This is a permutation-count limit - it depends
only on group sizes, so no amount of sequencing depth, cell number or
choice of statistic can beat it.

For a paired test over n units the floor is 2/2^n: six pairs floor at
0.031, three pairs at 0.25. Against a genome-wide editing site set of
~250,000 positions, Bonferroni sits near 2e-7, so per-site significance
is unreachable below roughly 22 pairs.

When the warning fires, the two correct responses are to restrict
`sites` to a small candidate list (recoding sites, a targeted panel), or
to use
[`computeAlleleFractionIndex`](https://github.com/cchmc/variantCell/reference/computeAlleleFractionIndex.md),
which spends a single test and so carries no genome-wide correction at
all.

## See also

[`computeAlleleFractionIndex`](https://github.com/cchmc/variantCell/reference/computeAlleleFractionIndex.md)
for the aggregate index,
[`checkGenotypeConcordance`](https://github.com/cchmc/variantCell/reference/checkGenotypeConcordance.md)
for the analogous guard on presence/absence contrasts

## Examples

``` r
if (FALSE) { # \dontrun{
project$setProjectIdentity("donor_type")

# Targeted: a candidate site list, where per-site testing can pay off
res <- project$findDEAlleleFraction("Donor", "Recipient",
                                    sites = recoding_sites)

# Genome-wide: warns that the ceiling cannot clear Bonferroni
res <- project$findDEAlleleFraction("Donor", "Recipient")
attr(res, "significance_reachable")
} # }
```
