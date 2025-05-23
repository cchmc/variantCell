#' @title setProjectIdentity: Set Project-Wide Cell Identity Variable
#' @name setProjectIdentity
#'
#' @description
#' Sets the metadata column to use as the primary cell identity variable for all downstream
#' analyses. This function establishes which cell grouping will be used in functions like
#' findDESNPs() and other SNP analysis methods.
#'
#' @param ident_col Character. Name of the column in cell_metadata to use as cell identity.
#'                 Must be an existing column in the metadata.
#'
#' @return Returns the object invisibly (for method chaining).
#'
#' @details
#' This function validates that the specified column exists in the cell metadata,
#' then sets it as the active identity for all subsequent analyses. It also prints
#' a summary of the unique identities found and the number of cells in each category.
#'
#' The identity column is crucial for many analysis functions as it defines how cells
#' are grouped when comparing SNP expression or presence between populations.
#'
#' Common identity variables include cell type annotations, cluster IDs, condition labels,
#' or any other categorical metadata that meaningfully separates cell populations.
#'
#' @note
#' - This function must be called before using analysis methods that rely on cell identities
#' - The function prints a summary of the identities found, which can be useful for verification
#' - The previously set identity (if any) is overwritten by this function
#'
#' @examples
#'
#' \dontrun{
#'
#' # Set cell type as the active identity
#'
#' project$setProjectIdentity("cell_type")
#'
#' # Use cluster IDs instead
#' project$setProjectIdentity("seurat_clusters")
#'
#' # Set disease status as identity for case-control comparisons
#' project$setProjectIdentity("disease_status")
#'
#' # Method chaining example
#' results <- project$setProjectIdentity("cell_type")$findDESNPs(
#'   ident.1 = "T_cells",
#'   donor_type = "Donor"
#' )
#'}
#' @seealso
#' \code{\link{getCurrentIdentity}} for checking the currently active identity
#' \code{\link{findDESNPs}} and \code{\link{findSNPsByGroup}} which use the set identity
variantCell$set("public",  "setProjectIdentity", function(ident_col) {
  # Check if identity exists in metadata
  if(!ident_col %in% colnames(self$snp_database$cell_metadata)) {
    stop(sprintf("Column '%s' not found in cell metadata", ident_col))
  }

  # Get unique identities
  unique_idents <- unique(self$snp_database$cell_metadata[[ident_col]])

  # Store current project-wide identity
  self$current_project_ident <- ident_col

  # Print summary
  cat(sprintf("\nSetting project-wide identity to: %s", ident_col))
  cat("\nUnique identities found:")
  for(ident in sort(unique_idents)) {
    n_cells <- sum(self$snp_database$cell_metadata[[ident_col]] == ident)
    cat(sprintf("\n  %s: %d cells", ident, n_cells))
  }

  invisible(self)
})
#' @title getCurrentIdentity: Get Current Project-Wide Cell Identity Information
#' @name getCurrentIdentity
#'
#' @description
#' Retrieves and displays information about the currently active cell identity variable.
#' This function returns details about which metadata column is being used for cell grouping
#' and provides a summary of the cell distribution across the different identity categories.
#'
#' @return Invisibly returns a list containing:
#'   \item{identity}{Character. Name of the current identity column.}
#'   \item{distribution}{Table. Distribution of cells across identity categories.}
#'   \item{total_cells}{Integer. Total number of cells in the dataset.}
#'   If no identity is set, returns NULL.
#'
#' @details
#' This function checks whether a project-wide identity has been set using
#' \code{setProjectIdentity()}. If an identity is active, it prints the name of the
#' identity column and displays a summary of how many cells belong to each identity category.
#'
#' The function is useful for:
#' - Verifying which grouping variable is currently active
#' - Checking the cell distribution across groups before analysis
#' - Confirming that identity assignments are as expected
#' - Debugging when analysis results seem unexpected
#'
#' @note
#' - If no identity has been set, the function returns NULL and displays a notification
#' - The function both prints information to the console and returns data invisibly
#' - The returned list can be captured and used programmatically if needed
#'
#' @examples
#'
#' \dontrun{
#'
#' # Check current identity
#' project$getCurrentIdentity()
#'
#' # Capture the return value for programmatic use
#' id_info <- project$getCurrentIdentity()
#' if(!is.null(id_info)) {
#'   # Find the most abundant cell type
#'   most_common <- names(which.max(id_info$distribution))
#'   cat("Most common cell type:", most_common)
#' }
#'
#' # Use in a workflow
#' project$setProjectIdentity("cell_type")
#' project$getCurrentIdentity()  # Verify it worked
#' project$findDESNPs(ident.1 = "T_cells")
#'}
#' @seealso
#' \code{\link{setProjectIdentity}} for setting the active identity
variantCell$set("public",  "getCurrentIdentity",  function() {
  if(is.null(self$current_project_ident)) {
    cat("No project-wide identity currently set\n")
    return(NULL)
  }

  # Get current identity information
  ident_col <- self$current_project_ident
  ident_values <- self$snp_database$cell_metadata[[ident_col]]

  cat(sprintf("\nCurrent identity: %s\n", ident_col))
  cat("\nCell distribution:")

  # Print distribution
  ident_table <- table(ident_values)
  for(ident in names(ident_table)) {
    cat(sprintf("\n  %s: %d cells", ident, ident_table[ident]))
  }

  # Return invisible summary
  invisible(list(
    identity = ident_col,
    distribution = ident_table,
    total_cells = length(ident_values)
  ))
})
#' @title subsetVariantCell: Subset Cells Based on Metadata Values
#' @name subsetVariantCell
#'
#' @description
#' Filters cells in the project based on values in a specified metadata column.
#' This function can either create a new variantCell object with the subset data (default)
#' or modify the current object in-place.
#'
#' @param column Character. Name of the metadata column to filter on.
#' @param values Vector. Values to include or exclude (depending on `invert` parameter).
#' @param invert Logical. If FALSE (default), keeps cells matching the values;
#'              if TRUE, keeps cells NOT matching the values.
#' @param copy Logical. If TRUE (default), returns a new variantCell object with the subset;
#'            if FALSE, modifies the current object in-place.
#'
#' @return
#' If copy=TRUE: Returns a new variantCell object containing only the subset data.
#' If copy=FALSE: Returns the modified object invisibly (for method chaining).
#'
#' @details
#' This function subsets all project data components, including:
#' - Alternative allele (AD) matrix
#' - Depth (DP) matrix
#' - Normalized depth matrix (if available)
#' - Cell metadata
#'
#' The function also handles sample management, removing samples that no longer
#' have any cells after filtering. This ensures data consistency throughout the object.
#'
#' With copy=TRUE (default), the original object remains untouched and a new object
#' with only the selected data is returned, allowing for exploration of subsets
#' without risk of data loss. With copy=FALSE, the subsetting operation permanently
#' modifies the original object, which is more memory-efficient but irreversible.
#'
#' @note
#' - When copy=TRUE, this function performs a deep clone which may use significant memory
#'   for large datasets
#' - In most bioinformatics workflows, it's recommended to keep the original data intact
#'   and work with copies to enable different analysis branches
#' - The function prints a summary of changes for verification
#'
#' @examples
#' \dontrun{
#' # Create a new variantCell object with only T cells (default behavior)
#' t_cell_project <- project$subsetvariantCell("cell_type", c("CD4", "CD8", "Treg"))
#'
#' # Modify the original project in-place (more memory efficient but irreversible)
#' project$subsetvariantCell("patient_id", "Patient1", copy = FALSE)
#'
#' # Create a subset excluding cells from a specific condition
#' no_acr_project <- project$subsetvariantCell("condition", "ACR", invert = TRUE)
#'
#' # Method chaining example (when using copy=FALSE)
#' results <- project$subsetvariantCell("donor_type", "Donor", copy = FALSE)$findDESNPs(
#'   ident.1 = "T_cells",
#'   ident.2 = "B_cells"
#' )
#'}
#' @seealso
#' \code{\link{aggregateByGroup}} for grouping cells after subsetting
variantCell$set("public", "subsetVariantCell", function(column, values, invert = FALSE, copy = TRUE) {
  # Make a deep copy if requested
  if(copy) {
    new_object <- self$clone(deep = TRUE)

    # Input validation
    if(!column %in% colnames(new_object$snp_database$cell_metadata)) {
      stop(sprintf("Column '%s' not found in metadata", column))
    }

    # Create mask for filtering
    if(invert) {
      cells_to_keep <- !new_object$snp_database$cell_metadata[[column]] %in% values
    } else {
      cells_to_keep <- new_object$snp_database$cell_metadata[[column]] %in% values
    }

    n_cells_before <- nrow(new_object$snp_database$cell_metadata)

    # Update matrices and metadata
    new_object$snp_database$ad_matrix <- new_object$snp_database$ad_matrix[, cells_to_keep]
    new_object$snp_database$dp_matrix <- new_object$snp_database$dp_matrix[, cells_to_keep]
    if(!is.null(new_object$snp_database$dp_matrix_normalized)) {
      new_object$snp_database$dp_matrix_normalized <- new_object$snp_database$dp_matrix_normalized[, cells_to_keep]
    }
    new_object$snp_database$cell_metadata <- new_object$snp_database$cell_metadata[cells_to_keep, ]

    # Check remaining samples
    remaining_samples <- unique(new_object$snp_database$cell_metadata$sample_id)
    samples_to_remove <- setdiff(names(new_object$samples), remaining_samples)

    if(length(samples_to_remove) > 0) {
      new_object$samples[samples_to_remove] <- NULL
      if(!is.null(new_object$metadata) && nrow(new_object$metadata) > 0) {
        new_object$metadata <- new_object$metadata[!new_object$metadata$sample_id %in% samples_to_remove, ]
      }
    }

    # Print summary
    cat(sprintf("\nSubset summary (copied object):"))
    cat(sprintf("\nFiltered by %s %s values: %s",
                if(invert) "excluding" else "keeping",
                column,
                paste(values, collapse=", ")))
    cat(sprintf("\nCells before: %d", n_cells_before))
    cat(sprintf("\nCells after: %d", nrow(new_object$snp_database$cell_metadata)))
    cat(sprintf("\nRemaining samples: %d", length(remaining_samples)))

    return(new_object)
  } else {
    # Original in-place implementation
    if(!column %in% colnames(self$snp_database$cell_metadata)) {
      stop(sprintf("Column '%s' not found in metadata", column))
    }

    # Create mask for filtering
    if(invert) {
      cells_to_keep <- !self$snp_database$cell_metadata[[column]] %in% values
    } else {
      cells_to_keep <- self$snp_database$cell_metadata[[column]] %in% values
    }

    n_cells_before <- nrow(self$snp_database$cell_metadata)

    # Update matrices and metadata
    self$snp_database$ad_matrix <- self$snp_database$ad_matrix[, cells_to_keep]
    self$snp_database$dp_matrix <- self$snp_database$dp_matrix[, cells_to_keep]
    if(!is.null(self$snp_database$dp_matrix_normalized)) {
      self$snp_database$dp_matrix_normalized <- self$snp_database$dp_matrix_normalized[, cells_to_keep]
    }
    self$snp_database$cell_metadata <- self$snp_database$cell_metadata[cells_to_keep, ]

    # Check remaining samples
    remaining_samples <- unique(self$snp_database$cell_metadata$sample_id)
    samples_to_remove <- setdiff(names(self$samples), remaining_samples)

    if(length(samples_to_remove) > 0) {
      self$samples[samples_to_remove] <- NULL
      if(!is.null(self$metadata) && nrow(self$metadata) > 0) {
        self$metadata <- self$metadata[!self$metadata$sample_id %in% samples_to_remove, ]
      }
    }

    # Print summary
    cat(sprintf("\nSubset summary (in-place):"))
    cat(sprintf("\nFiltered by %s %s values: %s",
                if(invert) "excluding" else "keeping",
                column,
                paste(values, collapse=", ")))
    cat(sprintf("\nCells before: %d", n_cells_before))
    cat(sprintf("\nCells after: %d", nrow(self$snp_database$cell_metadata)))
    cat(sprintf("\nRemaining samples: %d", length(remaining_samples)))

    invisible(self)
  }
})
#' @title downsampleVariant: Downsample Cells by Group to Balance Cell Numbers
#' @name downsampleVariant
#'
#' @description
#' Reduces the number of cells in the variantCell object by downsampling each group to a maximum
#' number of cells. This function is useful for balancing cell numbers across groups, reducing
#' computational burden, and mitigating the effects of groups with very different cell counts on
#' downstream analyses.
#'
#' @param max_cells Integer. Maximum number of cells to keep from each group. Groups with fewer
#'                 cells than this threshold will retain all their cells. Default: 1000.
#' @param group_by Character, optional. Metadata column to use for grouping cells. If NULL,
#'                uses the current project identity set by setProjectIdentity(). Default: NULL.
#' @param seed Integer. Random seed for reproducible downsampling. Default: 42.
#'
#' @return Returns the modified object invisibly (for method chaining).
#'
#' @details
#' The function performs downsampling by:
#' 1. Grouping cells based on the specified metadata column
#' 2. For each group, if cell count exceeds max_cells, randomly selecting max_cells cells to keep
#' 3. Updating all matrices and metadata to include only the selected cells
#' 4. Maintaining consistency across all data structures in the object
#'
#' This operation modifies the object in-place, permanently removing cells that aren't selected.
#' It's particularly useful when working with imbalanced datasets, where some cell types or
#' conditions have many more cells than others, which could bias analytical results.
#'
#' The function automatically handles updates to all relevant data structures, including:
#' - Alternative allele (AD) matrix
#' - Depth (DP) matrix
#' - Normalized depth matrix (if available)
#' - Cell metadata
#' - Sample-level information
#'
#' @note
#' - This function modifies the object in-place (no copy is created)
#' - Downsampling is performed randomly for each group
#' - The seed parameter ensures reproducibility of random sampling
#' - Groups with fewer cells than max_cells will keep all their cells
#' - If after downsampling a sample has no remaining cells, it will be removed from the object
#' - A detailed summary of the downsampling is printed to the console
#'
#' @examples
#' \dontrun{
#' # Basic usage - downsample to 500 cells per cell type
#' project$setProjectIdentity("cell_type")
#' project$downsampleVariant(max_cells = 500)
#'
#' # Downsample by a different grouping variable
#' project$downsampleVariant(
#'   max_cells = 200,
#'   group_by = "condition",
#'   seed = 123  # Use different seed for different random selection
#' )
#'
#' # Use with method chaining
#' results <- project$downsampleVariant(max_cells = 300)$findDESNPs(
#'   ident.1 = "T_cells",
#'   ident.2 = "B_cells"
#' )
#' }
#'
#' @seealso
#' \code{\link{setProjectIdentity}} for setting the grouping identity
#' \code{\link{subsetvariantCell}} for other filtering operations
variantCell$set("public", "downsampleVariant", function(max_cells = 1000,
                                                        group_by = NULL,
                                                        seed = 42) {
  # Input validation
  if(!is.numeric(max_cells) || max_cells <= 0) {
    stop("max_cells must be a positive number")
  }

  # Get metadata
  meta <- self$snp_database$cell_metadata
  n_cells_before <- nrow(meta)

  # Determine grouping column
  group_col <- group_by

  # If no group_by provided, use current project identity
  if(is.null(group_col)) {
    if(is.null(self$current_project_ident)) {
      stop("No grouping column specified and no current identity set. Either provide a group_by parameter or use setProjectIdentity() first.")
    }
    group_col <- self$current_project_ident
    cat(sprintf("\nNo grouping column specified, using current identity: %s\n", group_col))
  }

  # Verify column exists
  if(!group_col %in% colnames(meta)) {
    stop(sprintf("Column '%s' not found in metadata", group_col))
  }

  # Set seed for reproducibility
  set.seed(seed)

  # Identify cells to keep
  cells_to_keep <- character(0)

  # Get unique groups
  unique_groups <- unique(meta[[group_col]])
  cat(sprintf("\nDownsampling by '%s' with %d unique groups\n", group_col, length(unique_groups)))

  # Loop through groups and sample cells
  for(group in unique_groups) {
    group_cells <- rownames(meta)[meta[[group_col]] == group]
    n_group_cells <- length(group_cells)

    if(n_group_cells <= max_cells) {
      # Keep all cells for this group
      cat(sprintf("Group '%s': keeping all %d cells\n", group, n_group_cells))
      cells_to_keep <- c(cells_to_keep, group_cells)
    } else {
      # Randomly sample cells for this group
      sampled_cells <- sample(group_cells, max_cells)
      cat(sprintf("Group '%s': downsampled from %d to %d cells\n",
                  group, n_group_cells, max_cells))
      cells_to_keep <- c(cells_to_keep, sampled_cells)
    }
  }

  # Create a logical vector for subset indexing
  cell_mask <- rownames(meta) %in% cells_to_keep

  # Update matrices and metadata
  self$snp_database$ad_matrix <- self$snp_database$ad_matrix[, cell_mask]
  self$snp_database$dp_matrix <- self$snp_database$dp_matrix[, cell_mask]
  if(!is.null(self$snp_database$dp_matrix_normalized)) {
    self$snp_database$dp_matrix_normalized <- self$snp_database$dp_matrix_normalized[, cell_mask]
  }
  self$snp_database$cell_metadata <- self$snp_database$cell_metadata[cell_mask, ]

  # Check remaining samples and update as needed
  remaining_samples <- unique(self$snp_database$cell_metadata$sample_id)
  samples_to_remove <- setdiff(names(self$samples), remaining_samples)

  if(length(samples_to_remove) > 0) {
    self$samples[samples_to_remove] <- NULL
    cat(sprintf("\nRemoved %d samples with no remaining cells\n", length(samples_to_remove)))
  }

  # Update sample information
  for(sample_id in names(self$samples)) {
    # Update cell counts and filter metadata
    sample_cells <- rownames(meta)[meta$sample_id == sample_id & cell_mask]
    if(length(sample_cells) > 0) {
      self$samples[[sample_id]]$metadata <- self$samples[[sample_id]]$metadata[
        rownames(self$samples[[sample_id]]$metadata) %in% sample_cells, ]
    }
  }

  # Print summary
  n_cells_after <- nrow(self$snp_database$cell_metadata)
  cat(sprintf("\nDownsampling summary:"))
  cat(sprintf("\n - Cells before: %d", n_cells_before))
  cat(sprintf("\n - Cells after: %d", n_cells_after))
  cat(sprintf("\n - Reduction: %.2f%%", (1 - n_cells_after/n_cells_before) * 100))
  cat(sprintf("\n - Remaining samples: %d", length(remaining_samples)))

  # Print distribution after downsampling
  cat(sprintf("\n\nDistribution after downsampling by '%s':", group_col))
  group_table <- table(self$snp_database$cell_metadata[[group_col]])
  for(group in names(group_table)) {
    cat(sprintf("\n  %s: %d cells", group, group_table[group]))
  }

  invisible(self)
})
