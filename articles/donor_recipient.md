# Determining Donor and Recipient Cells

``` r

library(variantCell)
```

## The problem

Vireo separates the cells in a multiplexed library into genotype
clusters and labels them `donor0`, `donor1`, and so on. **Those labels
carry no meaning.** They are assigned independently on every run, and
re-processing the same data will happily swap them.

This is not hypothetical. Re-running a 22-sample lung transplant cohort
under a newer CellRanger version swapped the labels on **16 of 22
samples** — the same two cell populations, opposite labels, no error and
no warning.

A donor/recipient mapping recorded against one Vireo run must never be
reused against another. Doing so silently inverts Donor and Recipient,
and every downstream result inverts with it.

[`inferDonorType()`](https://github.com/cchmc/variantCell/reference/inferDonorType.md)
re-derives the mapping from the data every time.

## How it works

The graft supplies the organ’s structural tissue — epithelium,
endothelium, stroma. The recipient supplies infiltrating immune cells.
So for each Vireo cluster, compute the fraction of its cells that are
structural; the higher one is the donor.

``` r

sim <- simulateVariantCellData(verbose = FALSE)
```

The simulator randomises which genome Vireo calls `donor0`, exactly as a
real run does. Here is the ground truth it recorded:

``` r

sim$truth$labels
#>   sample patient donor0 donor1
#> 1     S1      P1   host  graft
#> 2     S2      P1   host  graft
#> 3     S3      P2  graft   host
#> 4     S4      P3  graft   host
#> 5     S5      P4  graft   host
```

Now infer it, using only the cell metadata and the Vireo output:

``` r

inf <- inferDonorType(
  sim$metadata,
  vireo_dir         = sim$vireo_dir,
  sample_col        = "Biopsy",
  barcode_col       = "Barcode",
  structural_col    = "compartment",
  structural_values = c("Epithelial", "Endothelial", "Stromal")
)
#> === Inferring donor/recipient identity ===
#> Samples: 5 | structural = Epithelial, Endothelial, Stromal
#>   S1       donor0=Recipient (struct 0.114, n=88) | donor1=Donor     (struct 0.903, n=72) | margin 0.789
#>   S2       donor0=Recipient (struct 0.071, n=85) | donor1=Donor     (struct 0.861, n=72) | margin 0.791
#>   S3       donor0=Donor     (struct 0.914, n=70) | donor1=Recipient (struct 0.034, n=88) | margin 0.880
#>   S4       donor0=Donor     (struct 0.875, n=72) | donor1=Recipient (struct 0.093, n=86) | margin 0.782
#>   S5       donor0=Donor     (struct 0.972, n=71) | donor1=Recipient (struct 0.070, n=86) | margin 0.902
#> 
#> Called 5/5 samples (0 ambiguous, 0 skipped)
#> Margin: min 0.782, median 0.791
#> donor0 = Recipient in 2 sample(s): S1, S2
#> 
#> Vireo donor labels are arbitrary per run - re-infer after any re-run.
```

The per-sample audit shows the separation it is working from:

``` r

inf$mapping[, c("sample", "n_donor0", "n_donor1",
                "struct_frac_donor0", "struct_frac_donor1",
                "margin", "donor0", "donor1", "call")]
#>   sample n_donor0 n_donor1 struct_frac_donor0 struct_frac_donor1    margin
#> 1     S1       88       72         0.11363636         0.90277778 0.7891414
#> 2     S2       85       72         0.07058824         0.86111111 0.7905229
#> 3     S3       70       88         0.91428571         0.03409091 0.8801948
#> 4     S4       72       86         0.87500000         0.09302326 0.7819767
#> 5     S5       71       86         0.97183099         0.06976744 0.9020635
#>      donor0    donor1 call
#> 1 Recipient     Donor   ok
#> 2 Recipient     Donor   ok
#> 3     Donor Recipient   ok
#> 4     Donor Recipient   ok
#> 5     Donor Recipient   ok
```

And it recovers the truth, including the sample whose labels are
flipped:

``` r

truth <- sim$truth$labels
got   <- inf$mapping[match(truth$sample, inf$mapping$sample), ]
data.frame(
  sample    = truth$sample,
  simulated = ifelse(truth$donor0 == "graft", "Donor", "Recipient"),
  inferred  = got$donor0,
  agree     = (truth$donor0 == "graft") == (got$donor0 == "Donor")
)
#>   sample simulated  inferred agree
#> 1     S1 Recipient Recipient  TRUE
#> 2     S2 Recipient Recipient  TRUE
#> 3     S3     Donor     Donor  TRUE
#> 4     S4     Donor     Donor  TRUE
#> 5     S5     Donor     Donor  TRUE
```

The separation is normally near-binary, which is why a *low* margin is
treated as a problem rather than mild uncertainty. Samples below
`min_margin` are reported `"ambiguous"` and left `NA` rather than
guessed:

``` r

range(inf$mapping$margin)
#> [1] 0.7819767 0.9020635
```

## Using the result

`inf$donor_types` is a per-sample named vector ready to hand to
[`addSampleData()`](https://github.com/cchmc/variantCell/reference/addSampleData.md):

``` r

inf$donor_types[["S1"]]
#>      donor0      donor1 
#> "Recipient"     "Donor"
inf$donor_types[["S3"]]
#>      donor0      donor1 
#>     "Donor" "Recipient"
```

``` r

project <- variantCell$new()
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
project$buildSNPDatabase()
```

With the mapping correct, the composition separates as it should:

``` r

db <- project$snp_database
table(db$cell_metadata$donor_type, db$cell_metadata$compartment)
#>            
#>             Endothelial Epithelial Lymphoid Myeloid Stromal
#>   Donor              60        207       23      11      56
#>   Recipient           8         19      239     161       6
```

## What getting it wrong looks like

Hard-coding `donor0 = "Donor"` for every sample produces this instead:

``` r

naive <- data.frame(
  sample   = inf$mapping$sample,
  inferred = inf$mapping$donor0,
  hardcoded = "Donor",
  result   = ifelse(inf$mapping$donor0 == "Donor", "correct", "INVERTED")
)
naive
#>   sample  inferred hardcoded   result
#> 1     S1 Recipient     Donor INVERTED
#> 2     S2 Recipient     Donor INVERTED
#> 3     S3     Donor     Donor  correct
#> 4     S4     Donor     Donor  correct
#> 5     S5     Donor     Donor  correct
```

``` r

sprintf("%d of %d samples would have donor and recipient swapped",
        sum(naive$result == "INVERTED"), nrow(naive))
#> [1] "2 of 5 samples would have donor and recipient swapped"
```

Sixteen of twenty-two in the real cohort. The failure is silent — the
analysis runs, the plots render, and the biology is backwards for those
samples.

## Confirming the clusters really are two genomes

The composition heuristic assumes Vireo split the library along
genotype. That assumption is checkable directly, without any cell-type
labels:

``` r

cm <- db$cell_metadata
gc_within <- checkGenotypeConcordance(
  db,
  which(cm$sample_id == "S1" & cm$donor_type == "Donor"),
  which(cm$sample_id == "S1" & cm$donor_type == "Recipient"),
  verbose = FALSE
)
c(concordance = round(gc_within$concordance, 3), verdict = gc_within$verdict)
#>         concordance             verdict 
#>             "0.499" "different_genomes"
```

Two genomes in one library agree at only about half the callable sites.
Compare that with the same individual sampled twice — `S1` and `S2` are
repeat biopsies of the same patient:

``` r

gc_same <- checkGenotypeConcordance(
  db,
  which(cm$sample_id == "S1" & cm$donor_type == "Recipient"),
  which(cm$sample_id == "S2" & cm$donor_type == "Recipient"),
  verbose = FALSE
)
c(concordance = round(gc_same$concordance, 3), verdict = gc_same$verdict)
#>   concordance       verdict 
#>       "0.999" "same_genome"
```

Near-perfect agreement, because it is one person. This is also a
sample-swap detector: the distinct-genome count comes from clustering on
concordance, not from any patient column, so a mislabelled sample cannot
hide behind its metadata.

## Other data types

[`inferDonorType()`](https://github.com/cchmc/variantCell/reference/inferDonorType.md)
accepts a Seurat object or a SingleCellExperiment as well as a data
frame, and with `add_metadata = TRUE` returns the input object with
`vireo_donor` and `donor_type` columns attached:

``` r

inf <- inferDonorType(merged_seurat, vireo_dir = "Vireo/cellranger8")
merged_seurat <- inf$object
```

If you only want Vireo’s assignments joined onto your metadata without
the donor/recipient call, the lower-level helpers do that directly:

``` r

res <- project$process_vireo_seurat(seurat_obj, vireo_path, prefix_text)
res <- project$process_vireo_sce(sce_obj, vireo_path, prefix_text)
res <- project$process_vireo_dataframe(meta_df, vireo_path, prefix_text)
```

These two chunks are not evaluated because they need your own object.

## Next

- **Building the SNP database** — the import path in full.
- **Group-level DE SNP analysis** — what a donor/recipient contrast can
  and cannot answer once you have it.
