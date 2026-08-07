# inferDonorType: Infer Donor vs Recipient from Vireo Assignments

Infers which Vireo genetic donor corresponds to the transplanted organ
("Donor") and which corresponds to the transplant recipient
("Recipient"), using the cell-type composition of each genotype cluster.

In a transplant biopsy the graft contributes the structural tissue
(epithelium, endothelium, stroma), while the recipient contributes
infiltrating immune cells. The genotype cluster holding the larger share
of structural cells is therefore the donor.

## Usage

``` r
inferDonorType(
  object,
  vireo_dir,
  sample_col = "Biopsy",
  barcode_col = "Barcode",
  structural_col = "compartment",
  structural_values = c("Epithelial", "Endothelial", "Stromal"),
  vireo_suffix = "_vireo",
  samples = NULL,
  min_cells = 50,
  min_margin = 0.1,
  add_metadata = TRUE,
  verbose = TRUE
)
```

## Arguments

- object:

  A Seurat object, SingleCellExperiment, or a data frame of cell
  metadata.

- vireo_dir:

  Character. Directory holding the per-sample Vireo output folders.

- sample_col:

  Character. Metadata column identifying the sample. Must match the
  Vireo folder names (after `vireo_suffix`). Default: "Biopsy".

- barcode_col:

  Character. Metadata column holding the plain cell barcode as it
  appears in Vireo's `donor_ids.tsv` (e.g. "AAACCCAAGATGCGAC-1"). For a
  merged object this is usually *not*
  [`colnames()`](https://rdrr.io/r/base/colnames.html). Default:
  "Barcode".

- structural_col:

  Character. Metadata column used to identify structural cells. Default:
  "compartment".

- structural_values:

  Character vector. Values of `structural_col` that count as
  graft-derived structural tissue. Default: c("Epithelial",
  "Endothelial", "Stromal").

- vireo_suffix:

  Character. Appended to the sample name to locate its Vireo folder
  (e.g. "\_vireo" gives "TBX1_vireo"). Default: "\_vireo".

- samples:

  Character vector. Restrict to these samples. Default: NULL (all
  samples found in `sample_col`).

- min_cells:

  Integer. Minimum matched cells required to call a sample. Default: 50.

- min_margin:

  Numeric. Minimum separation in structural fraction needed to make a
  call; below this the sample is flagged "ambiguous" and left NA.
  Default: 0.1.

- add_metadata:

  Logical. Attach `vireo_donor` and `donor_type` columns to `object` and
  return it. Default: TRUE.

- verbose:

  Logical. Print progress and a summary table. Default: TRUE.

## Value

A list with:

- mapping:

  Data frame, one row per sample: structural fractions, cell counts,
  `margin`, `call` ("ok"/"ambiguous"/"skipped"), and the inferred
  `donor0`/`donor1` roles.

- donor_types:

  Named list, one entry per confidently called sample, each a named
  vector `c(donor0 = ..., donor1 = ...)` ready to pass straight to
  `addSampleData(donor_type = ...)`.

- assignments:

  Data frame of per-cell assignments (sample, barcode, vireo donor,
  inferred donor_type).

- object:

  The input object with metadata columns added, when
  `add_metadata = TRUE`.

## Details

**Why this function is necessary.** Vireo's `donor0` / `donor1` labels
are assigned arbitrarily and independently on every run. They carry no
biological meaning and are *not* stable across re-runs of the same
sample: re-processing a sample (for example after a CellRanger upgrade)
will frequently swap them. A donor/recipient mapping recorded against
one Vireo run must never be reused against another. Always re-infer.

**Validation.** On a 22-sample lung transplant cohort with an
independently established mapping, this heuristic reproduced the known
assignment for 22/22 samples, with a median separation in structural
fraction of \>0.97.

**Interpreting `margin`.** `margin` is the absolute difference in
structural fraction between the two genotype clusters. Well-separated
samples score near 1.0. Samples below `min_margin` are reported as
`"ambiguous"` and are assigned `NA` rather than guessed, because an
inverted call silently swaps donor and recipient for every downstream
comparison. Always inspect the returned `mapping` table before use.

## Examples

``` r
if (FALSE) { # \dontrun{
inf <- inferDonorType(merged_obj, vireo_dir = "Vireo/cellranger8")

# ALWAYS review before use - especially any row flagged "ambiguous"
inf$mapping

# feed a confidently-called sample straight into addSampleData()
project$addSampleData(
  sample_id   = "TBX1",
  vireo_path  = "Vireo/cellranger8/TBX1_vireo",
  cellsnp_path= "cellSNP/cellranger8/TBX1",
  cell_data   = inf$object,
  donor_type  = inf$donor_types[["TBX1"]],
  prefix_text = "TBX1_"
)
} # }
```
