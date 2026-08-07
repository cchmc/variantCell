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

  Numeric 0-1. Minimum alternative allele fraction to consider a cell as
  expressing.

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
cells in each group. For each SNP, cells are only considered as
expressing if they have a minimum alternative allele fraction
(min_alt_frac) and positive read depth.

Statistical testing is performed using Wilcoxon rank-sum test when
calc_p=TRUE. Multiple testing correction is applied using the specified
p.adjust.method.

The parallel implementation distributes SNP processing across multiple
CPU cores for significantly improved performance on large datasets.

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
