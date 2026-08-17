# variantCell Package - Claude Code Documentation

## Package Overview

The variantCell package is an R package for analyzing single-cell
genetic variants, particularly focused on donor/recipient cell
identification in transplant genomics. The package provides
comprehensive tools for processing cellSNP data, integrating with vireo
donor assignments, and performing differential SNP analysis.

## Core Architecture

### R6 Class Structure

The package is built around an R6 class (`variantCell`) that manages: -
SNP databases from multiple samples - Cell metadata integration -
Normalization and analysis workflows - Population allele frequency
comparisons

### Key Components

#### 1. Data Processing (`R/02-process-vireo.R`, `R/03-snp-database.R`)

- **Vireo Integration**: Functions to process donor assignment data from
  vireo
  - [`process_vireo_seurat()`](https://github.com/cchmc/variantCell/reference/process_vireo_seurat.md):
    Integrate with Seurat objects
  - [`process_vireo_sce()`](https://github.com/cchmc/variantCell/reference/process_vireo_sce.md):
    Integrate with SingleCellExperiment objects
  - [`process_vireo_dataframe()`](https://github.com/cchmc/variantCell/reference/process_vireo_dataframe.md):
    Integrate with data frames
- **Sample Management**:
  [`addSampleData()`](https://github.com/cchmc/variantCell/reference/addSampleData.md)
  function for processing individual samples
- **Database Building**:
  [`buildSNPDatabase()`](https://github.com/cchmc/variantCell/reference/buildSNPDatabase.md)
  for unified SNP database creation

#### 2. SNP Analysis (`R/04-DE-SNP.R`)

- **Differential Analysis**:
  [`findDESNPs()`](https://github.com/cchmc/variantCell/reference/findDESNPs.md)
  for identifying differentially expressed SNPs
- **Group Analysis**:
  [`findSNPsByGroup()`](https://github.com/cchmc/variantCell/reference/findSNPsByGroup.md)
  with GroupOnly and SampleStratified modes
- **Population AF Support**: Integration with population allele
  frequency data from VCF files
- **Fold Enrichment**: Log2 fold change calculations comparing observed
  vs population AF

#### 3. Visualization (`R/05-plots.R`)

- **Heatmaps**:
  [`plotSNPHeatmap()`](https://github.com/cchmc/variantCell/reference/plotSNPHeatmap.md)
  for visualizing SNP patterns
- **SNP Plots**:
  [`plotSNPs()`](https://github.com/cchmc/variantCell/reference/plotSNPs.md)
  with population AF comparison support

#### 4. Utility Functions (`R/06-utils.R`)

- **Identity Management**:
  [`setProjectIdentity()`](https://github.com/cchmc/variantCell/reference/setProjectIdentity.md)
  and
  [`getCurrentIdentity()`](https://github.com/cchmc/variantCell/reference/getCurrentIdentity.md)
- **Cell Filtering**:
  [`subsetVariantCell()`](https://github.com/cchmc/variantCell/reference/subsetVariantCell.md)
  for cell subset operations
- **Downsampling**:
  [`downsampleVariant()`](https://github.com/cchmc/variantCell/reference/downsampleVariant.md)
  for balancing cell numbers
- **Cell Extraction**:
  [`getCellsForSNPs()`](https://github.com/cchmc/variantCell/reference/getCellsForSNPs.md)
  for identifying cells expressing specific SNPs

## Major Bug Fix - Memory Management for Large Samples (07/18/25)

### Problem Identified

Large samples (\>10,000 cells) with high SNP density were causing memory
failures during alt fraction calculation. The issue was particularly
evident in: - CLAD2: 324,330 SNPs × 13,124 cells → 23.3M non-zero
entries - CLAD4: 292,648 SNPs × 10,979 cells → 17.8M non-zero entries

### Root Cause

The sparse matrix operation
`alt_fractions[valid_mask] <- ad_matched[valid_mask] / dp_matched[valid_mask]`
was attempting to process millions of entries simultaneously, causing
memory exhaustion and silent failures.

### Solution Implemented

Added intelligent downsampling based on matrix density rather than fixed
cell thresholds:

1.  **New Parameter**: `max_nonzero_entries = 12000000` in
    [`addSampleData()`](https://github.com/cchmc/variantCell/reference/addSampleData.md)
2.  **Dynamic Calculation**: Automatically calculates optimal cell count
    to stay under threshold
3.  **Preserves SNP Coverage**: Preferentially downsamples cells over
    SNPs
4.  **User Configurable**: Threshold can be adjusted per analysis needs

``` r

# Usage examples:
project$addSampleData(...) # Uses default 12M threshold
project$addSampleData(..., max_nonzero_entries = 15000000) # More permissive
project$addSampleData(..., max_nonzero_entries = 8000000)  # More conservative
```

### Implementation Details

The fix is located in `R/03-snp-database.R` around lines 260-305:

``` r

# Check matrix density before processing
valid_mask <- dp_matched > 0
n_nonzero <- sum(valid_mask)

if(n_nonzero > max_nonzero_entries) {
  # Calculate target downsampling
  target_cells <- min(n_cells, round(max_nonzero_entries * n_cells / n_nonzero))
  target_cells <- max(target_cells, 1000)  # Minimum safety threshold
  
  # Downsample matrices and update all related objects
  # ... (downsampling logic)
}
```

### Results

- **CLAD2**: 23.3M → ~6,100 cells → ~12M entries ✓
- **CLAD4**: 17.8M → ~7,400 cells → ~12M entries ✓
- **Processing Time**: Reduced from hours to minutes
- **Memory Usage**: Stays within manageable limits
- **Statistical Power**: Maintains thousands of cells for robust
  analysis

## New Feature - Cell Extraction by SNP Expression (07/21/25)

### Function Added: `getCellsForSNPs()`

A new utility function was added to extract cell IDs based on SNP
expression criteria, enabling precise cell annotation in downstream
analysis tools like Seurat.

#### Key Features:

1.  **Flexible SNP Identification**: Accepts both chromosome:position
    format (e.g., “1:12345”) and rs IDs
2.  **Dual Threshold Control**:
    - `min_alt_frac` and `max_alt_frac` for alternative allele fraction
      range
    - `min_dp` for minimum read depth requirement
3.  **Quality Filtering**: Only includes cells meeting both expression
    and depth criteria
4.  **Sample Restriction**: Optional filtering to specific samples
5.  **Comprehensive Output**: Returns cell IDs plus summary statistics

#### Usage Examples:

``` r

# Find cells expressing SNPs above 20% alt fraction
expressing_cells <- project$getCellsForSNPs(
  snp_ids = c("rs12345", "2:67890"),
  min_alt_frac = 0.2,
  min_dp = 5
)

# Find cells NOT expressing the alternative allele (reference only)
ref_cells <- project$getCellsForSNPs(
  snp_ids = "rs10314",
  min_alt_frac = 0,
  max_alt_frac = 0,
  min_dp = 5
)

# Use results to annotate Seurat object
seurat_obj$snp_status <- ifelse(
  colnames(seurat_obj) %in% expressing_cells$rs12345,
  "expressing", "not_expressing"
)
```

#### Output Structure:

- **Named list**: Each SNP returns a vector of qualifying cell IDs
- **Summary statistics**: Attached as attributes showing:
  - Total cells with data (any reads for the SNP)
  - Cells meeting criteria (passing both alt_frac and depth filters)
  - Mean alternative allele fraction and depth
  - Success rate percentage

#### Quality Control:

The function distinguishes between: - **`total_cells_with_data`**: Cells
with any reads for the SNP (DP \> 0) - **`cells_meeting_criteria`**:
Cells passing both expression and minimum depth thresholds

This ensures reliable, high-quality cell annotations by filtering out
low-coverage measurements that might yield unreliable alternative allele
fractions.

#### Implementation Location:

The function is located in `R/06-utils.R` and integrates seamlessly with
the existing SNP database structure, properly handling the `snp_info`
metadata and `cell_metadata$cell_id` columns.

## Removed - prioritizeSNPs (08/12/26)

`prioritizeSNPs()` and its 13 helpers were deleted: **1,212 lines, 37%
of `04-DE-SNP.R` and 14.5% of the package.**

Reasons, in order of weight:

- **Four defects found 08/06, three of which produced silently wrong
  output** — scores attached to the wrong SNPs via a sorting
  [`merge()`](https://rdrr.io/r/base/merge.html), unannotated SNPs
  collapsed onto one score through
  [`match()`](https://rdrr.io/r/base/match.html) on NA, and
  `set.seed(42)` hijacking the caller’s global RNG. Every result it
  produced before that date was invalid.
- **Never had a test file.**
- The “ensemble ML regression” was a fixed weighted sum of hand-built
  features with no training, no held-out data and no validation.
- The LLM clinical assessment issued billable external API calls and
  emitted confident clinical claims — druggability scores, therapeutic
  potential — with no ground truth. Not defensible in a manuscript.

`ellmer` and `jsonlite` dropped from Suggests; `man/prioritizeSNPs.Rd`
and the NAMESPACE export removed. Suite still 102 assertions,
`R CMD check` Status: OK.

The prior “Critical Fixes to prioritizeSNPs” section is retained below
only as the record of why this was cut.

## Critical Fixes to findDESNPs (08/11/26)

Three defects, all active under the **default** `use_normalized = TRUE`.

### 1. Alt fractions were divided by log-normalized depth

`normalizeSnpCounts()` returns `log1p(DP * 10000 / colSums(DP))`, and
there is **no normalized AD matrix** anywhere in the package. The
per-SNP loop pulled DP from `dp_matrix_normalized` and AD from the raw
`ad_matrix`, then computed `alt_frac <- ad/dp`. That is a raw count over
a log-scaled depth, not a fraction — it reached **10.7** on the real
32-sample database.

Because `dp_norm` sits around 2–4 while `ad` is a raw count,
`>= min_alt_frac` was satisfied by almost any cell carrying a single alt
read, so **`min_alt_frac` was effectively disabled**. The
`mean_alt_frac1/2` output columns, which users read as allele fractions,
were wrong by the same factor.

Alt fractions and the expressing-cell gate are now always computed from
raw AD/DP regardless of `use_normalized`.

### 2. Effect size and p-value described different populations

`log2fc` summed depth over gated cells only; `wilcox.test(dp1, dp2)` ran
on all cells including zeros. `avg_expr1/2` now run over all cells in
the group, so the two agree. `min_alt_frac`/`min_expr_cells` select
which SNPs are *testable*, not which cells contribute to the effect size
— the Seurat convention.

### 3. A single qualifying SNP crashed the function

`dp_matrix_subset[qualifying_snp_indices, all_group_indices]` dropped to
a vector whenever exactly one SNP passed pre-filtering, and every
downstream `[i,]` failed with `incorrect number of dimensions`.
`drop = FALSE` added throughout.

### Scope

The SNP-level **pre-filter** used raw AD/DP and was correct, so SNP
*selection* was never affected. `findSNPsByGroup_GroupOnly` and
`_SampleStratified` use raw matrices throughout and are unaffected.
`findDESNPs` appears in no analysis script in `History/`, so **no
existing result is affected**.

Also removed a per-SNP [`gc()`](https://rdrr.io/r/base/gc.html) inside
the parallel worker loop.

### Tests

First `tests/` directory in the package:
`tests/testthat/test-findDESNPs.R`, 12 assertions across 5 tests
including serial/parallel agreement. `testthat` was already declared in
Suggests with edition 3.

### What findDESNPs is actually for

It tests **normalized read depth at a variant site**, which is
transcript abundance, not genotype. Documented accordingly: it is the
*coverage companion* to
[`findSNPsByGroup()`](https://github.com/cchmc/variantCell/reference/findSNPsByGroup.md)
— use it to confirm depth was adequate in the group where a SNP was
called absent, rather than as a variant-discovery method.

## New Feature - Genotype Concordance Guard (08/11/26)

### Function Added: `checkGenotypeConcordance()` (`R/08-genotype-qc.R`)

Germline presence/absence is a **genotype test**. It is informative
between different genomes, structurally empty within one genome, and
confounded when a group pools several. Nothing in a cell-identity label
says which situation the software is in, so this measures it from the
data before any test runs.

[`findSNPsByGroup()`](https://github.com/cchmc/variantCell/reference/findSNPsByGroup.md)
now calls it by default (`check_genotype = TRUE`), warns on the two bad
verdicts, and returns the verdict as `$genotype_check` so a flagged
contrast is still recognisable as flagged when the object is read back
later.

### Verdicts, calibrated on this cohort

| contrast | concordance | within-group | verdict |
|----|----|----|----|
| TBX60 Donor vs Recipient (one library) | 0.611 | — | `different_genomes` |
| TBX60 vs TBX66 Recipient (same patient) | 0.992 | — | `same_genome` |
| all Donor vs all Recipient (pooled) | 0.723 | 0.52 / 0.54 | `heterogeneous_groups` |

**The pooled row is the important one.** It is what
`findSNPsByGroup_GroupOnly()` does when the identity is `donor_type`,
and the within-group concordance (0.52) is **lower than the
between-group value (0.723)** — the groups are less internally
consistent than they are different from each other. A “Donor-specific”
SNP found that way has to be present in every donor and absent in every
recipient, which is an allele-frequency comparison between two arbitrary
sets of unrelated people.

**The only clean donor/recipient contrast is within a single sample**:
one donor genome against one recipient genome in one library. Neither
existing mode enforces this — `_GroupOnly` pools across samples, and
`_SampleStratified` requires consistency across unrelated pairs
(`sample_consistency = 0.3`).

### Cross-patient contrasts are valid but combinatorially capped

A pooled contrast across patients (ACR patients vs CLAD patients) is
**not invalid** — it is a case-control genetic association study, which
is a real design. What it cannot do is reach significance at this n, and
the reason is arithmetic rather than statistical:

    ACR (5 patients) vs CLAD (4)
      best attainable two-sided p, perfect separation : 0.0159   [= 2/C(9,5)]
      Bonferroni over 738,596 sites                   : 6.77e-08
      short by                                        : 234,000x

This is a **permutation-count ceiling**, the same failure that closed
the Patient 7 time series: it depends only on group sizes, so no
sequencing depth, cell count or better statistic can beat it. A
perfectly separating variant needs ~14 vs 14 patients merely to clear
Bonferroni, and real association studies run to hundreds because effects
are never perfectly separating.

So the guard **reports the ceiling instead of discouraging the
analysis**, and recommends ranking hits as exploratory candidates rather
than thresholding them. Two caveats survive at any n and ranking does
not fix them: unrelated individuals differ in **ancestry**, shifting
allele frequencies genome-wide independently of the grouping variable
(real GWAS corrects with ancestry PCs, unestimable at n=9); and here ACR
is chemistry-mixed while CLAD is 100% 3′, so the eligible-site
denominator differs between the arms.

`n_genomes1`/`n_genomes2` are obtained by single-linkage clustering of
sub-pseudobulks on concordance, **not** from a patient column — on the
real ACR-vs-CLAD contrast this recovered 5 vs 4 with no metadata at all,
which also means a mislabelled sample cannot inflate the count.

### Method

Best-covered autosomal sites only (sex chromosomes excluded, so a sex
difference cannot masquerade as genotype distance), pseudobulk alt
fractions binned at 0.15/0.85 into 0/1/2, pairwise concordance over
jointly callable sites. Each group is additionally split by `sample_id`
into sub-pseudobulks, so within-group heterogeneity is measured in the
same pass.

The guard defaults to 10,000 sites rather than 30,000: verdicts and
values are unchanged (0.612 / 0.990 / 0.718) at a third of the cost.

Validated against `History/2026-08-10_genotype-identity.R`, which used
the same method to establish 0.99-within / 0.63-between across all 64
sample x compartment pseudobulks. Tests in
`tests/testthat/test-checkGenotypeConcordance.R` (16 assertions,
deterministic fixtures, no RNG).

## New Feature - Donor/Recipient Inference (08/06/26)

### Function Added: `inferDonorType()` (`R/07-donor-inference.R`)

Infers which Vireo genotype cluster is the transplanted organ by
comparing the **structural-cell fraction** (Epithelial + Endothelial +
Stromal) between clusters. The graft supplies structural tissue; the
recipient supplies infiltrating immune cells.

``` r

inf <- inferDonorType(merged_obj, vireo_dir = "Vireo/cellranger8")
inf$mapping        # per-sample audit: fractions, counts, margin, call
inf$donor_types    # ready for addSampleData(donor_type = ...)
inf$object         # input object + donor_type / vireo_donor columns
```

**Validated 22/22** against an independently established mapping on the
lung transplant cohort. On the CellRanger 8 run it called 32/32 samples
with 0 ambiguous (min margin 0.500, median 0.939); Donor cells were 91%
structural and Recipient cells 99% immune.

### Why this function exists

**Vireo’s `donor0`/`donor1` labels are arbitrary and assigned
independently on every run.** They are not stable across re-processing.
Re-running this cohort under CellRanger 8 swapped the labels on **16 of
22** previously-mapped samples — CLAD1 went from
`donor0=2358/donor1=4020` to `donor0=4387/donor1=2477`, the same two
populations with opposite labels.

A donor/recipient mapping recorded against one Vireo run must **never**
be reused against another; doing so silently inverts Donor and Recipient
with no error. Always re-infer. Samples below `min_margin` are reported
`"ambiguous"` and left `NA` rather than guessed, because the separation
is normally near-binary (0.00 vs 0.99) — a low margin means something is
wrong, not merely uncertain.

## Analysis Caveat - 3′/5′ Chemistry Confound

10X chemistry is nested inside disease Stage in the lung cohort: **CLAD
is 100% 3′ and NoALAD/ALI/OP/ALAD\* are 100% 5′**. Only ACR is balanced
(3 vs 3).

This matters because the two chemistries interrogate different variant
space — 3′ reads the last few hundred bp of a transcript, 5′ reads the
TSS end:

| sample | chemistry | SNPs    |
|--------|-----------|---------|
| CLAD1  | 3v3       | 217,970 |
| TBX42  | 5v2       | 62,423  |
| TBX79  | 5v3       | 380,482 |

`CLAD1 ∩ TBX79 = 111,187` — only 51% of CLAD1 and 29% of TBX79. Two 5′
samples by contrast are nested (`TBX42 ∩ TBX79` = 97% of TBX42).

[`findSNPsByGroup()`](https://github.com/cchmc/variantCell/reference/findSNPsByGroup.md)
does **not** yield false positives from this: it requires
`dp >= min_depth` on both the present and the absent side, so uncovered
SNPs are dropped rather than miscalled as absent. The cost is power loss
and selection bias — cross-chemistry contrasts silently test only the
intersection, biased toward short/highly-expressed transcripts, and the
eligible-SNP denominator is not currently surfaced.

Note also that group `dp` is
[`Matrix::rowSums`](https://rdrr.io/pkg/Matrix/man/colSums-methods.html)
across samples, so `min_depth = 10` summed over 11 samples is a weak
per-sample guarantee. Use `findSNPsByGroup_SampleStratified` for
anything cross-chemistry, and lead with ACR for headline results.

## Documentation Build - IMPORTANT

**Do not run `devtools::document()` unless roxygen2 is pinned to 7.x.**

This package documents its R6 methods with roxygen blocks preceding
`variantCell$set("public", "method", ...)` calls. roxygen2 7.3.2 renders
those; **roxygen2 8.0.0 does not, and silently deletes them** — a single
`roxygenize()` removed 15 `man/*.Rd` files covering most of the public
API, and rewrote DESCRIPTION. Recovery is `git checkout -- man/` plus
restoring `RoxygenNote: 7.3.2`. Fixing the R6 doc blocks for roxygen 8
is separate, unstarted work.

## File Structure

    variantCell/
    ├── R/
    │   ├── 00-pkg-rqmts.r          # Package requirements
    │   ├── 01-init.r               # R6 class initialization
    │   ├── 02-process-vireo.R      # Vireo integration functions
    │   ├── 03-snp-database.R       # Core database management
    │   ├── 04-DE-SNP.R            # SNP analysis functions
    │   ├── 05-plots.R             # Visualization functions
    │   ├── 06-utils.R             # Utility functions
    │   ├── 07-donor-inference.R   # Donor/Recipient inference from Vireo + composition
    │   ├── 08-genotype-qc.R       # checkGenotypeConcordance - the three-regime guard
    │   ├── 09-allele-fraction.R   # AD/DP testing (editing, ASE, heteroplasmy)
    │   ├── 10-import-cellsnp.R    # non-germline cellsnp import (no vireo)
    │   └── 11-simulate.R          # simulateVariantCellData - makes examples runnable
    ├── History/                    # gitignored - interactive session logs
    ├── tests/testthat/             # 197 assertions across 6 files
    ├── vignettes/                  # all four executable as of 08/13/26
    │   ├── build_snp_database.Rmd # import path, addSampleData + buildSNPDatabase
    │   ├── donor_recipient.Rmd    # inferDonorType, validated against known truth
    │   ├── cell_level_DE.Rmd      # findDESNPs, plotSNPHeatmap, getCellsForSNPs
    │   └── group_level_DE.Rmd     # findSNPsByGroup, the genotype guard, plotSNPs
    └── CLAUDE.md                   # This documentation

## Usage Workflow

1.  **Infer donor identity**:
    [`inferDonorType()`](https://github.com/cchmc/variantCell/reference/inferDonorType.md)
    on the merged object — re-run this for every new Vireo run, never
    reuse a stored mapping
2.  **Initialize Project**: `project <- variantCell$new()`
3.  **Add Samples**: Use
    [`addSampleData()`](https://github.com/cchmc/variantCell/reference/addSampleData.md)
    for each sample with cellSNP and vireo data, passing
    `donor_type = inf$donor_types[[sample]]`
4.  **Build Database**:
    [`buildSNPDatabase()`](https://github.com/cchmc/variantCell/reference/buildSNPDatabase.md)
    to create unified SNP database
5.  **Set Identity**:
    [`setProjectIdentity()`](https://github.com/cchmc/variantCell/reference/setProjectIdentity.md)
    to define cell groupings
6.  **Check the regime**:
    [`checkGenotypeConcordance()`](https://github.com/cchmc/variantCell/reference/checkGenotypeConcordance.md)
    before any presence/absence contrast — it decides whether the
    question is answerable at all
7.  **Analyze**:
    [`findSNPsByGroup()`](https://github.com/cchmc/variantCell/reference/findSNPsByGroup.md)
    for genotype,
    [`findDESNPs()`](https://github.com/cchmc/variantCell/reference/findDESNPs.md)
    for depth at a *targeted* gene,
    [`computeAlleleFractionIndex()`](https://github.com/cchmc/variantCell/reference/computeAlleleFractionIndex.md)
    for AD/DP
8.  **Visualize**:
    [`plotSNPs()`](https://github.com/cchmc/variantCell/reference/plotSNPs.md)
    or
    [`plotSNPHeatmap()`](https://github.com/cchmc/variantCell/reference/plotSNPHeatmap.md)

Every step runs end to end in the vignettes against
[`simulateVariantCellData()`](https://github.com/cchmc/variantCell/reference/simulateVariantCellData.md),
which needs no external data.

### Note on `addSampleData(data_type = "dataframe")`

It requires `rownames(metadata_df) == prefix_text + barcode`, and
`vireo_path` points at the `donor_ids.tsv` **file**, not its directory.
Derive `prefix_text` by stripping the barcode off the existing cell name
(`substr(cn, 1, nchar(cn) - nchar(barcode))`) rather than rebuilding it
from metadata columns — cell names in a merged object are baked at merge
time and can carry stale labels that no longer match the current
metadata.

## Key Features

- **Multi-format Support**: Works with Seurat, SingleCellExperiment, and
  data frames
- **Population AF Integration**: Compare observed frequencies to
  population databases
- **Memory Efficient**: Intelligent downsampling for large datasets
- **Comprehensive Analysis**: From data import to statistical analysis
  and visualization
- **Transplant Focused**: Specialized tools for donor/recipient analysis

## Dependencies

The package integrates with the broader single-cell ecosystem including
Seurat, Matrix (sparse matrices), data.table (fast I/O), and standard R
statistical functions. \## New Module - Allele-Fraction Testing
(08/12/26)

### `R/09-allele-fraction.R`

The package tested genotype (`findSNPsByGroup`, presence/absence) and
read depth (`findDESNPs`, transcript abundance). **Neither tested
`AD/DP` itself** — the fraction of reads carrying the alternative base.
That is the statistic for anything where the variant lives in the RNA
rather than the DNA: A-to-I editing rate, allele-specific expression,
mtDNA heteroplasmy. The ASE work had to hand-roll it.

Two entry points, because they cost very different amounts of
statistical power:

| function | tests | multiple-testing burden |
|----|----|----|
| [`computeAlleleFractionIndex()`](https://github.com/cchmc/variantCell/reference/computeAlleleFractionIndex.md) | one index per pseudobulk | **none** — a single statistic |
| [`findDEAlleleFraction()`](https://github.com/cchmc/variantCell/reference/findDEAlleleFraction.md) | each site separately | Bonferroni over all sites |

**The index is the one that works at this cohort’s group sizes.** A
per-site scan over ~250,000 editing sites needs roughly 22 paired
samples before its permutation floor (2/2^n) clears Bonferroni (~2e-7).
Six pairs floor at 0.031, three at 0.25.
[`findDEAlleleFraction()`](https://github.com/cchmc/variantCell/reference/findDEAlleleFraction.md)
reports that ceiling and **warns when it cannot be cleared**, rather
than returning a table of hopeful p-values — the same philosophy as
[`checkGenotypeConcordance()`](https://github.com/cchmc/variantCell/reference/checkGenotypeConcordance.md).

### The unit of observation is the sample, never the cell

Cells within a library share a genome, a capture and a sequencing run.
Treating 8,154 of them as independent observations is pseudoreplication
and inflates significance by orders of magnitude — the classic
single-cell DE failure. Both functions aggregate to pseudobulk first via
a sparse indicator-matrix product, then test across pseudobulks. A
regression test asserts that a 10x increase in cell count leaves both
the index and the attainable p unchanged.

### `common_sites`

Defaults TRUE: within each `split_by` unit, both arms are restricted to
sites callable in all of that unit’s arms. Without it an index
difference can be manufactured purely by the two arms using different
positions.

The index is **read-weighted** (`sum(AD)/sum(DP)`) by default rather
than a mean of per-site fractions, which is dominated by low-coverage
sites.

### Tests

`tests/testthat/test-alleleFraction.R` — 26 assertions, deterministic
fixtures, no RNG. Full suite is now 64 assertions across 3 files.

### What to use it for

The strong designs in this cohort are **within-sample**, because every
contrast that has failed here failed on between-patient arithmetic.
Donor vs recipient inside one library is paired, so patient, ancestry,
chemistry and batch all cancel; cell-type contrasts likewise. Pass
`within = "cell_type"` for the stratified form.

## New Module - cellsnp-lite Import (08/12/26)

### `R/10-import-cellsnp.R`

[`addSampleData()`](https://github.com/cchmc/variantCell/reference/addSampleData.md) +
[`buildSNPDatabase()`](https://github.com/cchmc/variantCell/reference/buildSNPDatabase.md)
is the germline workflow: a cellSNP run paired with vireo donor
assignments. A run against a **non-germline region list** — RNA editing
sites, chrM, a candidate panel — has no vireo output, because there is
no genotype to deconvolve.
[`buildCellSNPDatabase()`](https://github.com/cchmc/variantCell/reference/buildCellSNPDatabase.md)
imports those runs directly.

``` r

project <- variantCellFromCellSNP(
  dirs = file.path("editing_pilot/run", samples), sample_ids = samples,
  cell_metadata = seurat_meta, identity = "donor_type")
idx <- project$computeAlleleFractionIndex(within = "cell_type")
```

Sites are **unioned** across samples, not intersected — different
samples cover different subsets (436k vs 470k in the pilot). Matrices
are assembled once from accumulated triplets; cbind-per-sample is
quadratic and has bitten this package before.

### OTH is carried through, and it matters

At an A\>G site the two remaining bases can only be sequencing error, so
`OTH/DP` is an **internal, position-matched error floor** for the
`AD/DP` signal. That comparison is what makes an editing result
interpretable rather than a plausible artifact — on the pilot it gave
7.4x (TBX3) and 9.5x (CLAD-2). `qc_report` reports `alt_fraction`,
`per_base_error_floor` and `signal_to_error`.

### Two bugs found by chasing test warnings

- **Leaked file connections.**
  [`Matrix::readMM()`](https://rdrr.io/pkg/Matrix/man/externalFormats.html)
  does not close a connection handed to it. Three leaked per sample
  against R’s ~128 connection cap would hard-fail a large import long
  before the GC reclaimed them. Pinned by a regression test.
- **Pattern matrices.** `writeMM` emits a `pattern` header whenever
  every stored value is 1, and `readMM` returns that as an `ngTMatrix`
  with **no `x` slot** — reading `@x` errors out. An OTH matrix from a
  shallow run really can be all ones. Same failure family as the
  pattern-matrix alt-fraction defect fixed earlier.

Neither surfaced as a test failure; both showed up only as accumulated
warnings. **Do not dismiss a clean pass with a warning tail.**

### Tests

`tests/testthat/test-importCellSNP.R` — 34 assertions. Suite is now **98
across 4 files**. `R CMD check` Status: OK.

### Build note

`R CMD build` needs pandoc for the vignettes and it is not on the shell
PATH. Prepend RStudio’s bundled copy:
`/usr/lib/rstudio/resources/app/bin/quarto/bin/tools/x86_64`

## API Cleanup - Internals Made Private (08/12/26)

Six methods were public without being API: `calculate_optimal_range`,
`calculate_y_positions`, `create_main_plot`, `create_distribution_plots`
(plotting internals), `process_tsv` (the Vireo TSV reader), and
`getNumericSubset` — which had **no call sites anywhere**, dead public
surface.

All moved to `private`, call sites updated to `private$`,
`man/process_tsv.Rd` removed. Public surface is now **25 members**, all
deliberate.

## Measured: What SNP-Site Depth Actually Captures (08/12/26)

The claim “Seurat dominates findDESNPs” was asserted for months and
finally **measured** on 2026-08-12, over the 52,809 cells shared between
the variantCell project and the merged Seurat object:

|                                       |                                  |
|---------------------------------------|----------------------------------|
| genes in Seurat matrix                | 27,765                           |
| genes with \>=1 covered SNP           | 16,572                           |
| **genes invisible to the SNP view**   | **11,193 (40.3%)**               |
| total reads: gene matrix vs SNP sites | 458.0M vs 234.9M = **51.3%**     |
| capture ratio per gene, median        | **0.202** (75th 0.46, 90th 0.90) |
| SNPs per gene, median / 99th          | **13 / 307**                     |

**The measurement is not thin** — half of all reads land at SNP
positions, and well-covered genes retain 90%. The earlier framing (“a
small fraction of the evidence”) was wrong.

**The domination is a multiple-testing result, not a read-count one.** A
median gene carries 13 SNPs, so it produces 13 correlated tests of one
fact, and correction runs over ~739,000 sites instead of ~20,000 genes —
**37x stricter**. Per test you hold ~1.5% of a gene’s evidence against
37x the penalty, on 60% of the genes.

So the corrected position, now in the roxygen and `man/findDESNPs.Rd`:

- **Not a genome-wide discovery scan.** Ever.
- **Legitimate for a targeted gene** — if you already know which gene
  matters the multiple-testing argument largely evaporates.
- **Legitimate as the coverage companion** to
  [`findSNPsByGroup()`](https://github.com/cchmc/variantCell/reference/findSNPsByGroup.md),
  and as a coverage-comparability check before
  [`computeAlleleFractionIndex()`](https://github.com/cchmc/variantCell/reference/computeAlleleFractionIndex.md).

**Unexploited and genuinely unique to this data:** depth *ratios between
sites within one gene* give 3’UTR length / alternative polyadenylation,
and ~12.5% of sites sit outside annotated genes. Gene counts cannot see
either. No current function extracts them — `findDESNPs` tests one site
against groups, not sites against each other. That is the one real gap
in this area.

## Vignettes Made Executable (08/13/26)

All four vignettes were `eval=FALSE` across **24 of 24 chunks**. They
pasted function *signatures* as if they were code, and two contained
`\dontrun{` — Rd syntax — inside R chunks with unbalanced braces.
Nothing in the repo ran the public API end to end.

Now **71 of 78 chunks execute**. The seven that do not need a
multi-gigabyte reference VCF or the reader’s own Seurat object, and each
says so in the prose.

### `simulateVariantCellData()` (`R/11-simulate.R`)

The enabling piece. Writes a synthetic experiment **as files** —
`cellSNP.base.vcf`, `cellSNP.samples.tsv`, AD/DP/OTH MatrixMarket
matrices, one `donor_ids.tsv` per sample — so vignettes exercise the
real parsing path rather than reaching into the object. 5 samples / 4
patients / 800 cells / ~2,000 sites, ~12 s end to end.

Design points that matter:

- **Sites sit at real exon coordinates** pulled from
  `EnsDb.Hsapiens.v86`, so `annotate_snps()` assigns real gene names and
  `plotSNPs(gene = "CD74")` works. The gene panel avoids symbols with
  alt-contig or readthrough duplicates (`HLA-A`, `B2M`, `PTPRC`,
  `EPCAM`) — those return several records and break
  [`plotSNPs()`](https://github.com/cchmc/variantCell/reference/plotSNPs.md),
  which calls
  [`start()`](https://rdrr.io/r/stats/start.html)/[`end()`](https://rdrr.io/r/stats/start.html)
  on the lookup result.
- **Genotypes are drawn under HWE from AF \> 0.05**, matching the
  common-variant panel constraint; per-cell alt counts are binomial with
  a 1% error rate, so nothing separates perfectly.
- **Which genome Vireo calls `donor0` is randomised per sample**,
  reproducing the real trap. Ground truth is returned in `$truth`.
- Depth carries a per-gene immune/structural shift, which is the signal
  [`findDESNPs()`](https://github.com/cchmc/variantCell/reference/findDESNPs.md)
  exists to detect.
- Denser than real 10x (~56% non-zero vs a few percent) so a small
  example still exercises everything.

Validation against the real cohort, all reproduced by the vignettes:
[`inferDonorType()`](https://github.com/cchmc/variantCell/reference/inferDonorType.md)
5/5 including the flipped sample; concordance 0.499 / 0.999 / 0.689
against the real 0.611 / 0.992 / 0.723;
[`findDESNPs()`](https://github.com/cchmc/variantCell/reference/findDESNPs.md)
recovers the simulated depth ordering at **r = 0.995**.

### Four defects found by running the code

All four were silent. None was a test failure.

**1. `annotate_snps()` attached transcript IDs to the wrong SNPs.**
`annotations$transcript_ids[chunk_indices] <- vapply(split(tx_id, chunk_indices), ...)`
indexed the LHS with one element per *overlap* and the RHS with one per
distinct *SNP*, so R recycled. Measured on a 16-SNP fixture spanning two
chromosomes: **10 of 16 SNPs carried transcripts from the other
chromosome.** Surfaced only as a “number of items to replace is not a
multiple of replacement length” warning on a passing run. `exon_ids` was
aligned but last-wins, discarding overlapping exons; both now collapse
via the [`split()`](https://rdrr.io/r/base/split.html) names.

**2. The genotype guard was inert in the normal case.**
`genotype_contrast_guard()` resolved idents against
`current_project_ident`, but the contrast is defined by whatever
[`aggregateByGroup()`](https://github.com/cchmc/variantCell/reference/aggregateByGroup.md)
grouped on — and aggregating by `condition` while the identity is
`cell_type` is the *common* case. The lookup matched no cells,
[`checkGenotypeConcordance()`](https://github.com/cchmc/variantCell/reference/checkGenotypeConcordance.md)
errored, the `tryCatch` swallowed it, and
`findSNPsByGroup(check_genotype = TRUE)` returned
`$genotype_check = NULL`. **The documented default protection did
nothing.** Now takes `group_col` and `donor_type` from
`aggregated_data$parameters`, and warns rather than returning NULL when
it cannot resolve.

**3.
[`getCellsForSNPs()`](https://github.com/cchmc/variantCell/reference/getCellsForSNPs.md)
never implemented its documented interface.** Every example in the
roxygen and `man/` used `"1:12345"`, but only `CHROM_POS_REF_ALT` and rs
IDs matched — so every documented example returned `character(0)` with a
warning. Added `chromosome:position` (with optional `chr` prefix).
Separately, an identifier resolving to several rows — which any rs ID at
a multi-allelic position does — indexed without `drop = FALSE` and
yielded a *matrix*, so `dp_row > 0`, `length(dp_row)` and the alt
fraction all ran over the whole matrix. Alt counts now pool across
alleles; DP is counted once.

**4.
[`plotSNPs()`](https://github.com/cchmc/variantCell/reference/plotSNPs.md)
clipped its own axis labels.** Group/split labels were anchored near the
right of the left gutter with `hjust = 1`, growing leftward past
`x_min`, where `scale_x_continuous(limits=)` clipped them. “Endothelial”
rendered as **“lial”** and “Recipient” as **“cipient”**. Now anchored at
the left edge with `hjust = 0` and a wider gutter. `05-plots.R` is 1,142
lines and had never been executed by anything.

Also silenced a `promoters()` out-of-bound warning that fired on every
[`buildSNPDatabase()`](https://github.com/cchmc/variantCell/reference/buildSNPDatabase.md)
call — it comes from unplaced scaffolds in the EnsDb gene set, not from
user data.

### Tests

`test-annotateSNPs.R` (new) and `test-getCellsForSNPs.R` (new); guard
cases appended to `test-checkGenotypeConcordance.R`. Suite went **102
-\> 197 assertions across 6 files**. `R CMD check` Status: OK.

### Version 0.1.3

Bumped for this work specifically because defect 1 **changes stored
output**. Any database built at 0.1.2 or earlier has `transcript_ids`
scrambled across SNPs, and `exon_ids` truncated to one arbitrary
overlapping exon. Without the bump, two databases both labelled 0.1.2
would differ with nothing to distinguish them.

`gene_name` was never affected — it was aligned, only last-wins on
overlapping genes — so
[`plotSNPs()`](https://github.com/cchmc/variantCell/reference/plotSNPs.md),
[`findDESNPs()`](https://github.com/cchmc/variantCell/reference/findDESNPs.md)
and everything gene-level is unaffected. **Rebuild only if something
downstream read `transcript_ids`.**

**The lesson is the one already in this file:** two of these four
surfaced only as warnings on a green suite, and two only when a human
looked at a rendered plot. Executable documentation is a test harness.

## Science Direction (as of 08/17/26)

**What the package is:** a multi-genome single-cell analysis layer on
top of cellsnp-lite/Vireo, with a transplant-specific identification
module. It consumes their output and replaces neither.
`07-donor-inference.R` is the only transplant-specific part — 349 of
7,665 lines. Everything else applies to any multi-donor experiment.

**Closed. Do not re-open without new data:**

| line | why |
|----|----|
| differential germline variants for disease state | genotype is fixed; 4 independent failures |
| presence/absence across patients | permutation ceiling, min p 0.057 at 3v4 |
| ASE | censored at the interesting end; needs pure recipient DNA |
| global editing index vs disease | capped, and ISG score beats it 7/7 vs 6/7 |
| passenger leukocyte decay | already published |
| donor-vs-recipient within cell type | 1 of 35 cells has both arms; graft/recipient is confounded with structural/immune |

**Open, ranked:**

1.  **Spatial chimerism** (`/mnt/d-drive/CLAD_spatial`, 189 GB, shipped
    Mar 2026, pipeline complete May 2026). Directly resolves the
    airway-restricted recipient epithelium question — sampling artifact
    vs anastomotic migration. Novel; nobody maps donor/recipient
    genotype spatially in lung transplant.

    **Not blocked** — corrected 08/17/26. The platform is **Curio
    TREKKER**, which is spatial *positioning* on top of an ordinary 10x
    single-nuclei 3′v4 run, not a spot-based assay: it maps 10x cell
    barcodes to bead spatial barcodes. Genotype and transcriptome both
    come from the standard `possorted_genome_bam.bam`, which is present
    for all four samples, and position joins on the cell barcode. So the
    existing pipeline runs unmodified — cellsnp-lite on the BAM, Vireo,
    [`inferDonorType()`](https://github.com/cchmc/variantCell/reference/inferDonorType.md),
    then attach `donor_type` to the TREKKER coordinates. **No per-spot
    deconvolution and no new method.** The earlier “blocked on locating
    the spatial BAMs” note was wrong on both counts.

    What actually limits it is n. Confidently positioned nuclei:
    20129-CTRL 7,695 (21.2%), 10292-CLAD 2,261 (10.8%), 20264-CTRL 750
    (32.0%), KM-CLAD 1,225 (34.3%). **Only the two CLAD samples are
    transplants** — the controls are non-transplant donor lung,
    single-genome, and
    [`checkGenotypeConcordance()`](https://github.com/cchmc/variantCell/reference/checkGenotypeConcordance.md)
    will correctly call them `same_genome`. At 3.5% recipient epithelium
    that is order **10–25 positioned recipient epithelial nuclei per
    CLAD sample**: enough for the geometric question (one contiguous
    airway vs dispersed through parenchyma), which is what actually
    distinguishes the two explanations, but n = 2 patients and it cannot
    carry more than that.

    **Assay caveat:** single-**nuclei** 3′v4 is a third chemistry for
    this project. Reads run 5.6–49.3% intronic and 3.4–15.6% intergenic
    with median UMI 782–3,167, so the variant catalog *expands* outside
    exons while per-cell depth drops. The 0.99 cross-chemistry
    concordance was measured on cytoplasmic 3′ vs 5′ only and does
    **not** transfer across the nuclei/whole-cell boundary unverified.

2.  **Compartment-specific chimerism vs clinical outcome.** One test per
    hypothesis, so the permutation ceiling does not bite. Clinical
    metadata is in the two LungChat .docx files. Check first whether the
    published decay is compartment-resolved.

3.  **Genotype QC as a standalone utility.**
    [`checkGenotypeConcordance()`](https://github.com/cchmc/variantCell/reference/checkGenotypeConcordance.md)
    recovers distinct-genome counts with no patient column, so a
    mislabelled sample cannot hide. Sample swaps are endemic in
    multi-donor single-cell work and there is no standard tool. Probably
    the most broadly citable piece after
    [`inferDonorType()`](https://github.com/cchmc/variantCell/reference/inferDonorType.md).

4.  **Cross-individual doublet detection.** Use Vireo’s calls, do not
    rebuild them — its mixture model is validated and variantCell does
    not do the clustering. The additive step is using genotype doublets
    to estimate the *total* doublet rate (inter-individual doublets are
    2pq of all doublets) and to benchmark scDblFinder against a real
    labelled set rather than simulation.

5.  **Recoding sites** — the one remaining differential idea, small
    enough that per-site testing is affordable under correction.
    `findDEAlleleFraction(sites=)`.

6.  **Within-gene depth ratios (APA / 3’UTR length), and intergenic
    sites.** Surfaced by the 08/12 measurement and *not currently
    reachable by any function* —
    [`findDESNPs()`](https://github.com/cchmc/variantCell/reference/findDESNPs.md)
    tests one site against groups, never sites against each other. A
    median gene carries 13 covered sites, so the ratio between them is a
    within-gene, within-cell measurement that needs no multiple-testing
    correction across genes and no between-patient comparison. ~12.5% of
    sites sit outside annotated genes, which gene counts cannot see at
    all. This is the only idea on the list that is unique to
    variant-level data rather than merely possible with it. Needs a new
    function; nothing to reuse.

**Engineering priority:** the vignette work is done (08/13/26) — 71 of
78 chunks execute, the suite went 102 → 197 assertions, and it
immediately found four silent defects. That closes the systemic reason
defects survived here.

What remains is the **roxygen 8 / R6 documentation debt**. `man/*.Rd`
and `NAMESPACE` are now hand-maintained, so every new export takes two
manual edits and one stray `devtools::document()` deletes 15 man files
(see the build note above). Fixing the R6 doc blocks so roxygen 8
renders them is unstarted, and it is the thing most likely to bite next.
It is still below the science items — it costs correctness of the
*docs*, not of the results.
