# findDESNPs: Cell-Level Differential SNP Expression Analysis

Identifies differentially expressed SNPs between cell populations by
comparing read depths and alternative allele frequencies. This function
performs comprehensive statistical analysis at the single-cell level,
with support for parallel processing to improve performance on large
datasets.

## Arguments

- ident.1:

  Character. Primary cell identity to analyze.

- ident.2:

  Character, optional. Secondary cell identity to compare against. If
  NULL, compares against all other cells.

- donor_type:

  Character, optional. Donor type to restrict analysis to ("Donor" or
  "Recipient"). If NULL, uses all cells regardless of donor type.

- use_normalized:

  Logical. Whether to use normalized depth counts (TRUE) or raw counts
  (FALSE).

- min_expr_cells:

  Integer. Minimum number of expressing cells required in each group.

- min_alt_frac:

  Numeric, 0 to 1. Minimum alternative allele fraction to consider a
  cell as expressing.

- logfc.threshold:

  Numeric. Minimum absolute log2 fold-change required to report a SNP.

- calc_p:

  Logical. Whether to calculate p-values (Wilcoxon test). Set to FALSE
  to save computation time.

- p.adjust.method:

  Character. Method for p-value adjustment, passed to p.adjust().
  Default: "BH" (Benjamini-Hochberg).

- return_all:

  Logical. Whether to return all SNPs or only significant ones.

- pseudocount:

  Numeric. Value added to expression values before log transformation.

- min.p:

  Numeric. Minimum p-value to report (prevents numerical underflow).

- debug:

  Logical. Whether to print debugging information during analysis.

- n_cores:

  Integer, optional. Number of CPU cores to use for parallel processing.
  If NULL, automatically uses detectCores() - 1.

- include_rs_ids:

  Logical. If true, includes RS IDs

- include_population_AF:

  Logical. If true, includes population allele frequencies

## Value

List containing:

- results:

  Data frame of differentially expressed SNPs with metrics including
  log2FC, expression values, cell counts, and significance statistics.

- summary:

  List with analysis overview, including counts of significant SNPs,
  up/downregulated SNPs, and parameter settings used.

## Details

The function calculates differential expression by comparing the average
expression of SNPs between two groups, normalized by the total number of
cells in each group. A cell counts as expressing a SNP if it has
positive read depth and an alternative allele fraction of at least
min_alt_frac. Alternative allele fractions are always computed from raw
AD/DP, never from normalized depth, regardless of use_normalized.

min_alt_frac and min_expr_cells select which SNPs are testable; they do
not restrict which cells contribute to the effect size.
avg_expr_group1/2 and the Wilcoxon test are both computed over all cells
in the group, so the reported log2fc and p-value describe the same
population.

Statistical testing is performed using Wilcoxon rank-sum test when
calc_p=TRUE. Multiple testing correction is applied using the specified
p.adjust.method.

Note that this tests normalized read depth at the variant site, which is
transcript abundance, not genotype. Allele identity never enters the
statistic; AD is used only to decide which SNPs are testable.

The parallel implementation distributes SNP processing across multiple
CPU cores for significantly improved performance on large datasets.

## What this measures, and when to use it

Read depth at a variant position is a real and correctly computed
measurement of the abundance of the transcript containing it. The
question is not whether it is reliable but whether it is the best
available instrument, and that depends entirely on the question.

Measured on 52,809 cells of the lung transplant cohort, against the
matched gene count matrix:

- SNP positions capture **51.3%** of all reads in the gene matrix
  (234.9M of 458.0M), so the measurement is not thin in aggregate.

- But **40.3%** of expressed genes (11,193 of 27,765) contain no covered
  SNP at all and are invisible to this view.

- The median gene retains only **20%** of its reads at SNP positions,
  though well-covered genes retain up to 90%.

- The median gene carries **13 SNPs** (99th percentile: 307), so it
  generates 13 correlated tests of one biological fact, and correction
  runs over ~739,000 sites rather than ~20,000 genes – roughly **37x
  stricter**.

**Do not use this as a genome-wide discovery scan.** Per test you have
on the order of 1.5% of a gene's evidence while paying 37x the
multiple-testing penalty, on 60% of the genes. Seurat's FindMarkers on
the gene count matrix strictly dominates it for "which genes differ
between these cell groups".

**It is appropriate for:**

- Targeted interrogation of a specific gene. If you already know which
  gene matters, the multiple-testing argument largely evaporates and a
  well-covered gene retains most of its reads.

- The coverage companion to
  [`findSNPsByGroup()`](https://github.com/cchmc/variantCell/reference/findSNPsByGroup.md):
  confirming depth was adequate in the group where a SNP was called
  absent, since otherwise "absent" only means "not detected".

- Checking coverage comparability before comparing allele fractions
  between groups with
  [`computeAlleleFractionIndex()`](https://github.com/cchmc/variantCell/reference/computeAlleleFractionIndex.md),
  since differential depth creates differential power.

Positional questions that gene counts genuinely cannot answer – 3' UTR
length and alternative polyadenylation from depth ratios between sites
within one gene, or transcription outside annotated genes – are real and
unique to this data, but **no current function extracts them**;
findDESNPs tests one site against groups, not sites against each other.

## Note

- Requires package 'parallel', 'foreach', and 'doParallel' for parallel
  processing

- Project identity must be set before using this function via
  setProjectIdentity()

- For non-transplant datasets, donor_type filtering is automatically
  disabled

## See also

[`setProjectIdentity`](https://github.com/cchmc/variantCell/reference/setProjectIdentity.md)
for setting the cell identity to use
[`findSNPsByGroup`](https://github.com/cchmc/variantCell/reference/findSNPsByGroup.md)
for group-level SNP analysis

## Examples

``` r

if (FALSE) { # \dontrun{
# Initialize a variantCell project

proj$setProjectIdentity('cell_type')

# Basic usage comparing T cells vs other cells, donor cells only
results <- proj$findDESNPs(
  ident.1 = "T_cells",
  ident.2 = NULL,
  donor_type = "Donor",
  min_expr_cells = 5,
  logfc.threshold = 0.25
)

# Without p-value calculation for faster processing
fast_results <- proj$findDESNPs(
  ident.1 = "CD4",
  ident.2 = "CD8",
  calc_p = FALSE,
  n_cores = 8
)

# Access results
head(results$results)
results$summary
} # }
```
