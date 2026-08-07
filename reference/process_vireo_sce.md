# process_vireo_sce: Integrate Vireo donor assignments with a SingleCellExperiment object

Processes donor genetic identity data from Vireo and integrates it with
a SingleCellExperiment object's metadata. The function adds donor
assignments to the colData, matches cell identifiers between Vireo
output and the SCE object, and generates summary statistics of cell type
distributions per donor.

## Arguments

- sce_obj:

  A SingleCellExperiment object. The single-cell dataset to be annotated
  with donor information.

- vireo_path:

  Character. Path to the Vireo donor_ids.tsv file containing donor
  assignments.

- prefix_text:

  Character. Text to prepend to cell identifiers in the Vireo data to
  match the cell barcodes in the SingleCellExperiment object.

## Value

A list containing:

- sce_object:

  The SingleCellExperiment object with donor assignment metadata added

- donor_data:

  Data frame containing donor assignments for matched cells

- matching_cells:

  Character vector of cell identifiers that matched between SCE and
  Vireo

- summaries:

  List of summary statistics including:

  - donor_summaries: Per-donor cell counts and cell type distributions
    (if available)

  - cells_matched: Total number of cells successfully matched

  - total_sce_cells: Total number of cells in the SingleCellExperiment
    object

  - total_vireo_cells: Total number of cells in the Vireo data

## Details

This function first processes the Vireo TSV file using the `process_tsv`
function, adding the prefix to cell identifiers. It then adds the donor
assignments to the SingleCellExperiment object's colData and generates
summary statistics of cell type distributions for each donor (if
cell_type information is available in colData).

## Note

The function checks for a "cell_type" column in the SingleCellExperiment
object's colData. If present, it will generate cell type distribution
summaries per donor. If not, the cell_types field in donor_summaries
will be NULL.

## Examples

``` r
if (FALSE) { # \dontrun{
# Process a SingleCellExperiment object with Vireo donor assignments
results <- process_vireo_sce(
  sce_obj = my_sce_object,
  vireo_path = "path/to/vireo/donor_ids.tsv",
  prefix_text = "Patient1_Sample3_"
)

# Access the updated SingleCellExperiment object
updated_sce <- results$sce_object

# Check matching statistics
results$summaries$cells_matched
results$summaries$total_sce_cells
} # }
```
