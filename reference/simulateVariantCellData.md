# Simulate a cellSNP-lite plus Vireo dataset

Writes a complete synthetic multi-genome experiment to disk in exactly
the layout
[`addSampleData()`](https://github.com/cchmc/variantCell/reference/addSampleData.md)
and
[`inferDonorType()`](https://github.com/cchmc/variantCell/reference/inferDonorType.md)
expect, and returns the paths and cell metadata needed to load it.
Intended for vignettes, examples, reproducible bug reports, and checking
an installation end to end.

## Usage

``` r
simulateVariantCellData(
  path = NULL,
  design = NULL,
  n_cells = 160,
  snps_per_gene = 180,
  genes = NULL,
  seed = 1,
  verbose = TRUE
)
```

## Arguments

- path:

  Directory to write to. Default: a new session temp directory.

- design:

  Data frame with columns `sample`, `patient` and `condition`. Samples
  sharing a `patient` share both genomes. Default: five samples across
  four patients, two of them repeat biopsies of the same patient.

- n_cells:

  Integer. Cells per sample. Default 160.

- snps_per_gene:

  Integer. Sites per gene. Default 180.

- genes:

  Data frame of genes with `gene`, `base_depth` and `immune_log2`
  columns, or NULL for the built-in panel of twelve lung and immune
  genes.

- seed:

  Integer seed. The draw does not disturb the caller's RNG stream.

- verbose:

  Logical. Print progress. Default TRUE.

## Value

A list with `root`, `cellsnp_paths`, `vireo_paths`, `vireo_dir`,
`prefixes`, `metadata` (row names are prefixed cell IDs, ready for
`addSampleData(data_type = "dataframe")`), `design`, `sites`, and
`truth` - the genome behind each Vireo label, plus the per-genome
genotype matrix.

## Details

Each simulated library contains two genetically distinct individuals.
The first supplies mostly structural cells and the second mostly immune
cells, reproducing the near-total lineage segregation seen in real
transplant data - which is what makes
[`inferDonorType()`](https://github.com/cchmc/variantCell/reference/inferDonorType.md)
work.

Which individual Vireo calls `donor0` is randomised per sample, because
that is what Vireo actually does: the labels are arbitrary and are
reassigned on every run. The ground truth is returned in `$truth` so a
vignette can show the consequence rather than assert it.

Genotypes are drawn under Hardy-Weinberg from allele frequencies above
0.05, matching the common-variant panels this package is normally run
against, and per-cell alt counts are binomial with a small
sequencing-error rate so that no site separates perfectly. Read depth
carries a per-gene shift between immune and structural cells, which is
the signal
[`findDESNPs()`](https://github.com/cchmc/variantCell/reference/findDESNPs.md)
is designed to detect.

The output is deliberately denser than a real 10x run - roughly half the
site-by-cell entries are non-zero, against a few percent in practice -
so that a small example still exercises every function.

## See also

[`inferDonorType`](https://github.com/cchmc/variantCell/reference/inferDonorType.md)
for recovering the donor/recipient mapping,
[`checkGenotypeConcordance`](https://github.com/cchmc/variantCell/reference/checkGenotypeConcordance.md)
for verifying the genotype structure.

## Examples

``` r
sim <- simulateVariantCellData(n_cells = 40, snps_per_gene = 10,
                               verbose = FALSE)
list.files(sim$cellsnp_paths[["S1"]])
#> [1] "cellSNP.base.vcf"    "cellSNP.samples.tsv" "cellSNP.tag.AD.mtx" 
#> [4] "cellSNP.tag.DP.mtx"  "cellSNP.tag.OTH.mtx"
head(sim$truth$labels)
#>   sample patient donor0 donor1
#> 1     S1      P1   host  graft
#> 2     S2      P1  graft   host
#> 3     S3      P2  graft   host
#> 4     S4      P3  graft   host
#> 5     S5      P4   host  graft
```
