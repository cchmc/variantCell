# Check whether two cell groups are genetically distinct

Measures pseudobulk genotype concordance between two sets of cells, and
within each set, to establish whether a presence/absence contrast
between them is a genotype comparison at all.

Germline presence/absence is a genotype test. It is informative when the
two groups are different genomes (donor versus recipient within a
sample, or two individuals), structurally empty when they share a
genome, and confounded when either group pools several genomes – in that
last case hits reflect which individuals landed in which group, not the
grouping variable.

## Usage

``` r
checkGenotypeConcordance(
  snp_database,
  cells1,
  cells2,
  split_by = NULL,
  n_sites = 30000,
  min_dp = 10,
  min_cells = 50,
  min_shared_sites = 500,
  hom_ref = 0.15,
  hom_alt = 0.85,
  same_genome_threshold = 0.9,
  different_genome_threshold = 0.8,
  verbose = TRUE
)
```

## Arguments

- snp_database:

  The `snp_database` element of a variantCell project.

- cells1, cells2:

  Cell identifiers (matching `cell_metadata$cell_id`) or integer column
  indices defining the two groups.

- split_by:

  Optional vector, one entry per cell in `cell_metadata`, used to split
  each group into sub-pseudobulks so within-group genome heterogeneity
  can be measured. Defaults to `sample_id` when present.

- n_sites:

  Number of best-covered autosomal sites to use. Default 30000.

- min_dp:

  Minimum pseudobulk depth for a site to be called. Default 10.

- min_cells:

  Minimum cells for a sub-pseudobulk to be used. Default 50.

- min_shared_sites:

  Minimum jointly callable sites for a comparison. Default 500.

- hom_ref, hom_alt:

  Alt-fraction bounds separating homozygous reference, heterozygous and
  homozygous alternative. Defaults 0.15 and 0.85.

- same_genome_threshold:

  Concordance at or above which two pseudobulks are the same individual.
  Default 0.90.

- different_genome_threshold:

  Concordance at or below which they are different individuals. Default
  0.80.

- verbose:

  Print a summary. Default TRUE.

## Value

A list with `concordance` (between the two groups), `n_shared_sites`,
`within1`/`within2` (minimum concordance among each group's own
sub-pseudobulks, NA if not splittable), `n_genomes1`/`n_genomes2`
(distinct individuals per group, obtained by single-linkage clustering
of sub-pseudobulks on concordance, so no patient column is needed),
`min_attainable_p`, `bonferroni_threshold`, `significance_reachable`,
`verdict`, and `message`. `verdict` is one of `"different_genomes"`,
`"same_genome"`, `"heterogeneous_groups"`, or `"indeterminate"`.

## Details

Only autosomal sites are used, so the result is unaffected by sex. Sites
are ranked by total depth across the whole database and the best-covered
`n_sites` retained, which keeps the aggregation cheap and the calls well
supported.

A `"heterogeneous_groups"` verdict means at least one group contains
more than one individual. That is a legitimate design – a case-control
genetic association study – but the evidence available is capped by the
number of individuals, not by sequencing. For `n1` versus `n2`
individuals the smallest attainable two-sided p for a perfectly
separating variant is `2/choose(n1 + n2, n1)`, a permutation-count
ceiling that no amount of depth, cell number or statistical refinement
can beat. The function reports that ceiling against a Bonferroni
threshold so it is clear up front whether any hit could reach
significance; when it cannot, rank hits as exploratory candidates rather
than thresholding them.

Two caveats survive at any n. Unrelated individuals differ in ancestry,
which shifts allele frequencies genome-wide independently of the
grouping variable (association studies correct for this with ancestry
principal components, which cannot be estimated at small n). And where
the grouping variable is confounded with assay chemistry, the
eligible-site denominator differs between the groups.

Calibration on the 32-sample lung transplant cohort: concordance runs
~0.99 between compartments of the same individual and ~0.61 between
donor and recipient within a sample. Pooling donor cells across patients
drops the within-group concordance to ~0.52, below the between-group
value – the groups are then less internally consistent than they are
different from each other.

## See also

[`findSNPsByGroup`](https://github.com/cchmc/variantCell/reference/findSNPsByGroup.md),
which runs this as a guard by default.

## Examples

``` r
if (FALSE) { # \dontrun{
# Valid contrast: two genomes in one library
checkGenotypeConcordance(
  project$snp_database,
  cells1 = which(cm$sample_id == "TBX60" & cm$donor_type == "Donor"),
  cells2 = which(cm$sample_id == "TBX60" & cm$donor_type == "Recipient"))
#> concordance 0.611, verdict different_genomes

# Structurally empty: same person, two biopsies
checkGenotypeConcordance(
  project$snp_database,
  cells1 = which(cm$sample_id == "TBX60" & cm$donor_type == "Recipient"),
  cells2 = which(cm$sample_id == "TBX66" & cm$donor_type == "Recipient"))
#> concordance 0.992, verdict same_genome

# Confounded: donor_type pooled across patients
checkGenotypeConcordance(
  project$snp_database,
  cells1 = which(cm$donor_type == "Donor"),
  cells2 = which(cm$donor_type == "Recipient"))
#> concordance 0.723, within-group 0.52 / 0.54, verdict heterogeneous_groups
} # }
```
