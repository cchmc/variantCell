# Building the SNP database

``` r

library(variantCell)
```

## What this vignette needs

variantCell sits on top of two upstream tools and replaces neither:

- **cellsnp-lite**, which counts reference and alternative reads per
  cell at a list of variant positions, and
- **Vireo**, which clusters those counts into genotypes and assigns each
  cell to a genome.

Every function here consumes their output. So that this vignette
actually runs,
[`simulateVariantCellData()`](https://github.com/cchmc/variantCell/reference/simulateVariantCellData.md)
writes a synthetic experiment in exactly the layout a real run produces
— `cellSNP.base.vcf`, `cellSNP.samples.tsv`, the AD/DP/OTH MatrixMarket
matrices, and one `donor_ids.tsv` per sample.

``` r

sim <- simulateVariantCellData(verbose = FALSE)

list.files(sim$cellsnp_paths[["S1"]])
#> [1] "cellSNP.base.vcf"    "cellSNP.samples.tsv" "cellSNP.tag.AD.mtx" 
#> [4] "cellSNP.tag.DP.mtx"  "cellSNP.tag.OTH.mtx"
```

Five libraries, two genetically distinct individuals in each. Two of
them (`S1`, `S2`) are repeat biopsies of the same patient, which will
matter later.

``` r

sim$design
#>   sample patient condition
#> 1     S1      P1 Rejection
#> 2     S2      P1 Rejection
#> 3     S3      P2 Rejection
#> 4     S4      P3    Stable
#> 5     S5      P4    Stable
```

The cell metadata is an ordinary data frame whose **row names are the
prefixed cell IDs**. In a real analysis this is your Seurat or
SingleCellExperiment metadata.

``` r

head(sim$metadata, 3)
#>                                  Barcode Biopsy patient condition
#> S1_AAGCCCTATTGGGTAG-1 AAGCCCTATTGGGTAG-1     S1      P1 Rejection
#> S1_ATAACCCACCGTTTGT-1 ATAACCCACCGTTTGT-1     S1      P1 Rejection
#> S1_GGGCGGGTGAACTGCT-1 GGGCGGGTGAACTGCT-1     S1      P1 Rejection
#>                                cell_type compartment
#> S1_AAGCCCTATTGGGTAG-1           Ciliated  Epithelial
#> S1_ATAACCCACCGTTTGT-1           Ciliated  Epithelial
#> S1_GGGCGGGTGAACTGCT-1 Classical Monocyte     Myeloid
```

> The simulated data is deliberately denser than a real 10x run —
> roughly half the site-by-cell entries are non-zero, against a few
> percent in practice — so that a small example still has enough
> coverage to demonstrate every function. Cell counts and site counts
> are correspondingly small.

## Creating a project

``` r

project <- variantCell$new()
```

## Adding samples

[`addSampleData()`](https://github.com/cchmc/variantCell/reference/addSampleData.md)
reads one sample’s cellsnp matrices, matches barcodes against your cell
metadata, applies the Vireo donor assignment, and stores the result.

The `donor_type` argument maps Vireo’s cluster labels onto meaningful
names. **Do not hard-code it.** Vireo’s `donor0` / `donor1` labels are
arbitrary and are reassigned on every run — see the donor/recipient
vignette, which covers this in full. Infer the mapping from the data
instead:

``` r

inf <- inferDonorType(sim$metadata, vireo_dir = sim$vireo_dir, verbose = FALSE)
```

``` r

inf$donor_types[["S1"]]
#>      donor0      donor1 
#> "Recipient"     "Donor"
inf$donor_types[["S3"]]
#>      donor0      donor1 
#>     "Donor" "Recipient"
```

Those two are not the same mapping, which is the whole point.

``` r

for (s in sim$design$sample) {
  project$addSampleData(
    sample_id    = s,
    vireo_path   = sim$vireo_paths[[s]],
    cellsnp_path = sim$cellsnp_paths[[s]],
    cell_data    = sim$metadata[sim$metadata$Biopsy == s, ],
    data_type    = "dataframe",
    prefix_text  = sim$prefixes[[s]],
    donor_type   = inf$donor_types[[s]]
  )
}
```

``` r

names(project$samples)
#> [1] "S1" "S2" "S3" "S4" "S5"
```

A few notes on the arguments:

- `data_type` is `"seurat"`, `"sce"` or `"dataframe"`. For a data frame,
  `vireo_path` points at the `donor_ids.tsv` **file**, not its
  directory, and row names must equal `prefix_text` + barcode.
- `prefix_text` is stripped from your cell names to recover the plain
  barcode that cellsnp-lite wrote. Derive it by removing the barcode
  from an existing cell name rather than rebuilding it from metadata
  columns — cell names in a merged object are fixed at merge time and
  can carry stale labels.
- For a first pass, leave `min_cells = 0` and `min_alt_frac = 0`.
  Filtering early discards sites that downstream plotting and DE need.
- `non_transplant_mode = TRUE` skips Vireo entirely, for a single-genome
  experiment.

## Building the unified database

[`buildSNPDatabase()`](https://github.com/cchmc/variantCell/reference/buildSNPDatabase.md)
takes the union of sites across samples, annotates them against
`EnsDb.Hsapiens.v86`, and assembles one matrix per quantity.

``` r

project$buildSNPDatabase()
```

``` r

db <- project$snp_database
dim(db$dp_matrix)
#> [1] 1994  800
names(db)
#> [1] "ad_matrix"            "dp_matrix"            "dp_matrix_normalized"
#> [4] "cell_metadata"        "snp_info"             "snp_annotations"     
#> [7] "snp_metrics"          "qc_report"
```

Sites are **unioned**, not intersected — different samples cover
different positions, and requiring a site in every sample would discard
most of them.

``` r

head(db$snp_info, 3)
#>      CHROM      POS REF ALT         snp_id rs_id population_AF
#> S1.1     1 27666066   A   T 1_27666066_A_T  <NA>            NA
#> S1.2     1 27666069   T   A 1_27666069_T_A  <NA>            NA
#> S1.3     1 27666071   G   T 1_27666071_G_T  <NA>            NA
```

Annotation lands in `snp_metrics`, which is where the gene name lives:

``` r

table(db$snp_metrics$feature_type)
#> 
#> exonic 
#>   1994
head(sort(table(db$snp_metrics$gene_name), decreasing = TRUE), 5)
#> 
#>  CLDN5    DCN   FCN1   JAK1 PECAM1 
#>    180    180    180    180    180
```

The cell metadata carries everything you supplied plus the columns
variantCell adds — `sample_id`, `donor_id`, `donor_type`:

``` r

colnames(db$cell_metadata)
#>  [1] "cell_id"     "sample_id"   "donor_id"    "donor_type"  "Barcode"    
#>  [6] "Biopsy"      "patient"     "condition"   "cell_type"   "compartment"
table(db$cell_metadata$donor_type, db$cell_metadata$compartment)
#>            
#>             Endothelial Epithelial Lymphoid Myeloid Stromal
#>   Donor              60        207       23      11      56
#>   Recipient           8         19      239     161       6
```

That table is the single most important property of transplant data: the
graft supplies structural cells and the recipient supplies immune cells,
almost completely. It is what makes donor inference reliable — and it is
also why phenotypic donor-versus-recipient comparisons are largely
unavailable, since the two genomes rarely coexist within one cell type.

## Setting an identity

Most analysis functions work off a project-wide identity column:

``` r

project$setProjectIdentity("compartment")
```

``` r

project$getCurrentIdentity()
#> 
#> Current identity: compartment
#> 
#> Cell distribution:
#>   Endothelial: 68 cells
#>   Epithelial: 229 cells
#>   Lymphoid: 265 cells
#>   Myeloid: 175 cells
#>   Stromal: 63 cells
```

## Adding rs IDs and population allele frequencies

If you have a reference VCF — 1000 Genomes, gnomAD —
[`buildSNPDatabase()`](https://github.com/cchmc/variantCell/reference/buildSNPDatabase.md)
can attach rs identifiers and population allele frequencies:

``` r

project$buildSNPDatabase(
  add_rs_ids        = TRUE,
  add_population_AF = TRUE,
  VCF_file_path     = "path/to/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf"
)
```

This chunk is not evaluated because it needs a multi-gigabyte reference
file that cannot ship with the package. Everything else in this vignette
runs.

Note what such a panel implies: if you called variants against a
common-variant reference, **every site in your database is common by
construction**. Nothing rare, and nothing high-penetrance. That is a
property of the panel, not of the method — cellsnp-lite will run against
any region list, which is what
[`buildCellSNPDatabase()`](https://github.com/cchmc/variantCell/reference/buildCellSNPDatabase.md)
exists to import.

## Next

- **Determining donor and recipient cells** — why the Vireo labels must
  be re-inferred every run.
- **Cell-level DE SNP analysis** —
  [`findDESNPs()`](https://github.com/cchmc/variantCell/reference/findDESNPs.md),
  and what it does and does not measure.
- **Group-level DE SNP analysis** —
  [`findSNPsByGroup()`](https://github.com/cchmc/variantCell/reference/findSNPsByGroup.md),
  and the genotype guard that decides whether your contrast is
  answerable at all.
