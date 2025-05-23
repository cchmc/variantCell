#'@title process_tsv: Process a TSV file with prefix addition
#'@name process_tsv
#'
#' @description
#' Reads a TSV file, adds a prefix to a specified column, and sets that column
#' as row names in the resulting data frame. This function is typically used to preprocess
#' sample metadata files and ensure cell identifiers match across different data sources.
#'
#' @param file_path Character. Path to the TSV file to be read.
#' @param prefix_column Character. Name of the column to which the prefix should be added
#'   and which will be used as row names.
#' @param prefix_text Character. Text to prepend to each value in the specified column.
#'
#' @return A data frame where:
#'   - The specified column has been prefixed with the given text
#'   - The modified column has been set as row names
#'   - The column itself has been removed from the data frame
#'
#' @details
#' This function is particularly useful when processing Vireo donor assignment files
#' that need to be integrated with other single-cell data (like Seurat objects), where
#' cell identifiers need to match exactly for proper merging.
#'
#' @examples
#' \dontrun{
#' # Process a donor assignment file for Vireo
#' donor_meta <- process_tsv(
#'   file_path = "path/to/donor_ids.tsv",
#'   prefix_column = "cell",
#'   prefix_text = "Patient1_Sample3_"
#' )
#'
#' # Now donor_meta has the prefixed cell IDs as row names
#' # and can be easily matched with Seurat object cell barcodes
#' }
#'
#' @note
#' The function assumes that the file exists and is formatted as a proper TSV file.
#' No validation is performed to check if the column specified by `prefix_column` exists
#' in the TSV file.
variantCell$set("public","process_tsv", function(file_path, prefix_column, prefix_text) {
  # Read the TSV file
  df <- read.delim(file_path, sep="\t", stringsAsFactors=FALSE)

  # Add prefix to specified column
  df[[prefix_column]] <- paste0(prefix_text, df[[prefix_column]])

  # Set the modified column as row names
  rownames(df) <- df[[prefix_column]]

  # Remove the column since it's now in row names
  df[[prefix_column]] <- NULL

  return(df)
})
#' @title process_vireo_seurat: Integrate Vireo donor assignments with a Seurat object
#' @name process_vireo_seurat
#'
#' @description
#' Processes donor genetic identity data from Vireo and integrates it with a Seurat object's metadata.
#' The function adds donor assignments to each cell's metadata, matches cell identifiers between
#' Vireo output and the Seurat object, and generates summary statistics of cell type distributions
#' per donor.
#'
#' @param seurat_obj A Seurat object. The single-cell dataset to be annotated with donor information.
#' @param vireo_path Character. Path to the Vireo donor_ids.tsv file containing donor assignments.
#' @param prefix_text Character. Text to prepend to cell identifiers in the Vireo data to match
#'   the cell barcodes in the Seurat object.
#'
#' @return A list containing:
#'   \item{seurat_object}{The Seurat object with donor assignment metadata added}
#'   \item{donor_data}{Data frame containing donor assignments for matched cells}
#'   \item{matching_cells}{Character vector of cell identifiers that matched between Seurat and Vireo}
#'   \item{summaries}{List of summary statistics including:
#'     \itemize{
#'       \item donor_summaries: Per-donor cell counts and cell type distributions
#'       \item cells_matched: Total number of cells successfully matched
#'       \item total_seurat_cells: Total number of cells in the Seurat object
#'       \item total_vireo_cells: Total number of cells in the Vireo data
#'     }
#'   }
#'
#' @details
#' This function first processes the Vireo TSV file using the `process_tsv` function,
#' adding the prefix to cell identifiers. It then adds the donor assignments to the Seurat
#' object's metadata and generates summary statistics of cell type distributions for each donor.
#' The function prints verbose diagnostic information during execution to help track the
#' matching process.
#'
#' @note
#' The function assumes that the Seurat object has a "cell_type" column in its metadata,
#' which is used to generate the cell type distribution summaries per donor.
#'
#' @examples
#' \dontrun{
#' # Process a Seurat object with Vireo donor assignments
#' results <- process_vireo_seurat(
#'   seurat_obj = my_seurat_object,
#'   vireo_path = "path/to/vireo/donor_ids.tsv",
#'   prefix_text = "Patient1_Sample3_"
#' )
#'
#' # Access the updated Seurat object
#' updated_seurat <- results$seurat_object
#'
#' # View donor summaries
#' results$summaries$donor_summaries
#'
#' # Check matching statistics
#' results$summaries$cells_matched
#' results$summaries$total_seurat_cells
#' }
variantCell$set("public","process_vireo_seurat", function(seurat_obj, vireo_path, prefix_text) {

  seurat_add_metadata <- function(seurat_obj, new_metadata) {
    # Print initial information
    cat("Seurat object cells:", nrow(seurat_obj@meta.data), "\n")
    cat("Total Vireo cells:", nrow(new_metadata), "\n")

    # Check first few rownames from each
    cat("\nFirst few Seurat rownames:", head(rownames(seurat_obj@meta.data)), "\n")
    cat("First few metadata rownames:", head(rownames(new_metadata)), "\n")

    # Find exact matches
    matching_cells <- intersect(rownames(seurat_obj@meta.data), rownames(new_metadata))
    cat("\nNumber of matching cells:", length(matching_cells), "\n")

    if(length(matching_cells) > 0) {
      # Direct assignment for each column
      for(col in colnames(new_metadata)) {
        seurat_obj@meta.data[[col]] <- NA  # Initialize with NA
        seurat_obj@meta.data[matching_cells, col] <- new_metadata[matching_cells, col]
      }

      # Verify addition
      cat("\nVerification - first few rows of updated metadata:\n")
      print(head(seurat_obj@meta.data[matching_cells, colnames(new_metadata)]))
    }

    return(seurat_obj)
  }


  # Print initial cell identities
  cat("\nInitial cell identities in Seurat object:\n")
  print(table(Idents(seurat_obj)))

  # Process the vireo data using process_tsv helper function
  donor_data <- self$process_tsv(vireo_path, "cell", prefix_text)

  # Add metadata using debug_and_add_metadata helper function
  seurat_obj <- seurat_add_metadata(seurat_obj, donor_data)

  # Get matching cells
  matching_cells <- intersect(rownames(seurat_obj@meta.data), rownames(donor_data))

  # Generate donor summaries with cell identities
  donor_summaries <- lapply(unique(donor_data$donor_id), function(donor) {
    # Get cells for this donor using direct indexing
    donor_cells <- matching_cells[donor_data[matching_cells, "donor_id"] == donor]

    cat(sprintf("\nCell type distribution for %s:\n", donor))
    cell_type_dist <- table(seurat_obj@meta.data[donor_cells, "cell_type"])
    print(cell_type_dist)

    list(
      total_cells = length(donor_cells),
      cell_types = cell_type_dist
    )
  })
  names(donor_summaries) <- unique(donor_data$donor_id)

  return(list(
    seurat_object = seurat_obj,
    donor_data = donor_data[matching_cells, , drop=FALSE],
    matching_cells = matching_cells,
    summaries = list(
      donor_summaries = donor_summaries,
      cells_matched = length(matching_cells),
      total_seurat_cells = ncol(seurat_obj),
      total_vireo_cells = nrow(donor_data)
    )
  ))
})
#' @title process_vireo_sce: Integrate Vireo donor assignments with a SingleCellExperiment object
#' @name process_vireo_sce
#'
#' @description
#' Processes donor genetic identity data from Vireo and integrates it with a SingleCellExperiment
#' object's metadata. The function adds donor assignments to the colData, matches cell identifiers
#' between Vireo output and the SCE object, and generates summary statistics of cell type
#' distributions per donor.
#'
#' @param sce_obj A SingleCellExperiment object. The single-cell dataset to be annotated with
#'   donor information.
#' @param vireo_path Character. Path to the Vireo donor_ids.tsv file containing donor assignments.
#' @param prefix_text Character. Text to prepend to cell identifiers in the Vireo data to match
#'   the cell barcodes in the SingleCellExperiment object.
#'
#' @return A list containing:
#'   \item{sce_object}{The SingleCellExperiment object with donor assignment metadata added}
#'   \item{donor_data}{Data frame containing donor assignments for matched cells}
#'   \item{matching_cells}{Character vector of cell identifiers that matched between SCE and Vireo}
#'   \item{summaries}{List of summary statistics including:
#'     \itemize{
#'       \item donor_summaries: Per-donor cell counts and cell type distributions (if available)
#'       \item cells_matched: Total number of cells successfully matched
#'       \item total_sce_cells: Total number of cells in the SingleCellExperiment object
#'       \item total_vireo_cells: Total number of cells in the Vireo data
#'     }
#'   }
#'
#' @details
#' This function first processes the Vireo TSV file using the `process_tsv` function,
#' adding the prefix to cell identifiers. It then adds the donor assignments to the
#' SingleCellExperiment object's colData and generates summary statistics of cell type
#' distributions for each donor (if cell_type information is available in colData).
#'
#' @note
#' The function checks for a "cell_type" column in the SingleCellExperiment object's colData.
#' If present, it will generate cell type distribution summaries per donor. If not, the
#' cell_types field in donor_summaries will be NULL.
#'
#' @examples
#' \dontrun{
#' # Process a SingleCellExperiment object with Vireo donor assignments
#' results <- process_vireo_sce(
#'   sce_obj = my_sce_object,
#'   vireo_path = "path/to/vireo/donor_ids.tsv",
#'   prefix_text = "Patient1_Sample3_"
#' )
#'
#' # Access the updated SingleCellExperiment object
#' updated_sce <- results$sce_object
#'
#' # Check matching statistics
#' results$summaries$cells_matched
#' results$summaries$total_sce_cells
#' }
variantCell$set("public",  "process_vireo_sce", function(sce_obj, vireo_path, prefix_text) {
  # Read the TSV file
  donor_data <- self$process_tsv(vireo_path, "cell", prefix_text)

  # Find matching cells between SCE and donor data
  matching_cells <- intersect(colnames(sce_obj), rownames(donor_data))

  cat(sprintf("\nFound %d matching cells between SingleCellExperiment and Vireo data",
              length(matching_cells)))

  # Add donor information to colData
  if(length(matching_cells) > 0) {
    colData(sce_obj)$donor_id <- NA
    colData(sce_obj)[matching_cells, "donor_id"] <- donor_data[matching_cells, "donor_id"]
  }

  # Generate donor summaries
  donor_summaries <- lapply(unique(donor_data$donor_id), function(donor) {
    # Get cells for this donor
    donor_cells <- matching_cells[donor_data[matching_cells, "donor_id"] == donor]

    # Get cell type information if available
    cell_type_dist <- NULL
    if("cell_type" %in% colnames(colData(sce_obj))) {
      cell_type_dist <- table(colData(sce_obj)[donor_cells, "cell_type"])
      cat(sprintf("\nCell type distribution for %s:\n", donor))
      print(cell_type_dist)
    }

    list(
      total_cells = length(donor_cells),
      cell_types = cell_type_dist
    )
  })
  names(donor_summaries) <- unique(donor_data$donor_id)

  return(list(
    sce_object = sce_obj,
    donor_data = donor_data[matching_cells, , drop=FALSE],
    matching_cells = matching_cells,
    summaries = list(
      donor_summaries = donor_summaries,
      cells_matched = length(matching_cells),
      total_sce_cells = ncol(sce_obj),
      total_vireo_cells = nrow(donor_data)
    )
  ))
})
#' @title process_vireo_dataframe: Integrate Vireo donor assignments with a metadata data frame
#' @name process_vireo_dataframe
#'
#' @description
#' Processes donor genetic identity data from Vireo and integrates it with an existing
#' metadata data frame. The function matches cell identifiers between the Vireo output
#' and the metadata, and generates basic cell count statistics per donor.
#'
#' @param metadata_df A data frame. The metadata data frame to match with donor information,
#'   with cell identifiers as row names.
#' @param vireo_path Character. Path to the Vireo donor_ids.tsv file containing donor assignments.
#' @param prefix_text Character. Text to prepend to cell identifiers in the Vireo data to match
#'   the cell identifiers in the metadata data frame.
#'
#' @return A list containing:
#'   \item{metadata}{The original metadata data frame, unchanged}
#'   \item{donor_data}{Data frame containing donor assignments for matched cells}
#'   \item{matching_cells}{Character vector of cell identifiers that matched between metadata and Vireo}
#'   \item{summaries}{List of summary statistics including:
#'     \itemize{
#'       \item donor_summaries: Per-donor cell counts
#'       \item cells_matched: Total number of cells successfully matched
#'       \item total_cells: Total number of cells in the metadata data frame
#'       \item total_vireo_cells: Total number of cells in the Vireo data
#'     }
#'   }
#'
#' @details
#' This function provides a simpler alternative to the Seurat and SingleCellExperiment integrations
#' when you only have a data frame of metadata. It first processes the Vireo TSV file using the
#' `process_tsv` function, adding the prefix to cell identifiers, then finds matching cells
#' between the Vireo data and the metadata data frame, and generates basic summary statistics.
#'
#' @note
#' Unlike the Seurat and SingleCellExperiment integration functions, this function does not
#' modify the input metadata data frame. It only returns the matching information.
#'
#' @examples
#' \dontrun{
#' # Process a metadata data frame with Vireo donor assignments
#' results <- process_vireo_dataframe(
#'   metadata_df = cell_metadata,
#'   vireo_path = "path/to/vireo/donor_ids.tsv",
#'   prefix_text = "Patient1_Sample3_"
#' )
#'
#' # Check matching statistics
#' results$summaries$cells_matched
#' results$summaries$total_cells
#'
#' # Access donor assignments for matched cells
#' donor_assignments <- results$donor_data
#' }
variantCell$set("public",  "process_vireo_dataframe", function(metadata_df, vireo_path, prefix_text) {
  # Read the TSV file
  donor_data <- self$process_tsv(vireo_path, "cell", prefix_text)

  # Find matching cells between metadata and donor data
  matching_cells <- intersect(rownames(metadata_df), rownames(donor_data))

  cat(sprintf("\nFound %d matching cells between metadata and Vireo data", length(matching_cells)))

  # Generate summary stats
  donor_summaries <- lapply(unique(donor_data$donor_id), function(donor) {
    # Get cells for this donor
    donor_cells <- matching_cells[donor_data[matching_cells, "donor_id"] == donor]

    list(
      total_cells = length(donor_cells)
      # Add any other summary stats
    )
  })
  names(donor_summaries) <- unique(donor_data$donor_id)

  return(list(
    metadata = metadata_df,
    donor_data = donor_data[matching_cells, , drop=FALSE],
    matching_cells = matching_cells,
    summaries = list(
      donor_summaries = donor_summaries,
      cells_matched = length(matching_cells),
      total_cells = nrow(metadata_df),
      total_vireo_cells = nrow(donor_data)
    )
  ))
})
