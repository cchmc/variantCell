# process_vireo_seurat: Integrate Vireo donor assignments with a Seurat object

Processes donor genetic identity data from Vireo and integrates it with
a Seurat object's metadata. The function adds donor assignments to each
cell's metadata, matches cell identifiers between Vireo output and the
Seurat object, and generates summary statistics of cell type
distributions per donor.

## Arguments

- seurat_obj:

  A Seurat object. The single-cell dataset to be annotated with donor
  information.

- vireo_path:

  Character. Path to the Vireo donor_ids.tsv file containing donor
  assignments.

- prefix_text:

  Character. Text to prepend to cell identifiers in the Vireo data to
  match the cell barcodes in the Seurat object.

## Value

A list containing:

- seurat_object:

  The Seurat object with donor assignment metadata added

- donor_data:

  Data frame containing donor assignments for matched cells

- matching_cells:

  Character vector of cell identifiers that matched between Seurat and
  Vireo

- summaries:

  List of summary statistics including:

  - donor_summaries: Per-donor cell counts and cell type distributions

  - cells_matched: Total number of cells successfully matched

  - total_seurat_cells: Total number of cells in the Seurat object

  - total_vireo_cells: Total number of cells in the Vireo data

## Details

This function first processes the Vireo TSV file using the `process_tsv`
function, adding the prefix to cell identifiers. It then adds the donor
assignments to the Seurat object's metadata and generates summary
statistics of cell type distributions for each donor. The function
prints verbose diagnostic information during execution to help track the
matching process.

## Note

The function assumes that the Seurat object has a "cell_type" column in
its metadata, which is used to generate the cell type distribution
summaries per donor.

## Examples

``` r
if (FALSE) { # \dontrun{
# Process a Seurat object with Vireo donor assignments
results <- process_vireo_seurat(
  seurat_obj = my_seurat_object,
  vireo_path = "path/to/vireo/donor_ids.tsv",
  prefix_text = "Patient1_Sample3_"
)

# Access the updated Seurat object
updated_seurat <- results$seurat_object

# View donor summaries
results$summaries$donor_summaries

# Check matching statistics
results$summaries$cells_matched
results$summaries$total_seurat_cells
} # }
```
