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
#' @param drop_empty_snps Logical. If TRUE (default), SNPs left with no coverage in the
#'            retained cells are removed, along with their rows in snp_info,
#'            snp_annotations and snp_metrics. Subsets otherwise inherit the parent's full
#'            SNP dimension - subsetting a 32-sample database to 2 samples keeps all
#'            738,596 rows, ~81% of them empty, and every downstream call pays for them.
#'            Set FALSE to keep the SNP dimension aligned with the parent object.
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
variantCell$set("public", "subsetVariantCell", function(column, values, invert = FALSE, copy = TRUE,
                                                        drop_empty_snps = TRUE) {
  # Drop SNP rows that have no coverage left after the cell filter, keeping every
  # SNP-indexed table aligned. Without this a subset keeps the full SNP dimension
  # of the parent database - subsetting 32 samples down to 2 still carried all
  # 738,596 rows, ~81% of them empty, and every downstream call paid for them.
  prune_snps <- function(db) {
    keep <- Matrix::rowSums(db$dp_matrix) > 0
    if(all(keep)) return(db)
    db$ad_matrix <- db$ad_matrix[keep, , drop = FALSE]
    db$dp_matrix <- db$dp_matrix[keep, , drop = FALSE]
    if(!is.null(db$dp_matrix_normalized)) {
      db$dp_matrix_normalized <- db$dp_matrix_normalized[keep, , drop = FALSE]
    }
    # These are all indexed by the same SNP order as the matrices; any one of
    # them left unpruned would silently misalign every result.
    for(tbl in c("snp_info", "snp_annotations", "snp_metrics")) {
      if(!is.null(db[[tbl]]) && nrow(db[[tbl]]) == length(keep)) {
        db[[tbl]] <- db[[tbl]][keep, , drop = FALSE]
      }
    }
    attr(db, "n_snps_dropped") <- sum(!keep)
    db
  }

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

    n_snps_before <- nrow(new_object$snp_database$dp_matrix)
    if(drop_empty_snps) {
      new_object$snp_database <- prune_snps(new_object$snp_database)
    }

    # Print summary
    cat(sprintf("\nSubset summary (copied object):"))
    cat(sprintf("\nFiltered by %s %s values: %s",
                if(invert) "excluding" else "keeping",
                column,
                paste(values, collapse=", ")))
    cat(sprintf("\nCells before: %d", n_cells_before))
    cat(sprintf("\nCells after: %d", nrow(new_object$snp_database$cell_metadata)))
    if(drop_empty_snps) {
      cat(sprintf("\nSNPs: %d -> %d (dropped %d with no coverage)",
                  n_snps_before, nrow(new_object$snp_database$dp_matrix),
                  n_snps_before - nrow(new_object$snp_database$dp_matrix)))
    }
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

    n_snps_before <- nrow(self$snp_database$dp_matrix)
    if(drop_empty_snps) {
      self$snp_database <- prune_snps(self$snp_database)
    }

    # Print summary
    cat(sprintf("\nSubset summary (in-place):"))
    cat(sprintf("\nFiltered by %s %s values: %s",
                if(invert) "excluding" else "keeping",
                column,
                paste(values, collapse=", ")))
    cat(sprintf("\nCells before: %d", n_cells_before))
    cat(sprintf("\nCells after: %d", nrow(self$snp_database$cell_metadata)))
    if(drop_empty_snps) {
      cat(sprintf("\nSNPs: %d -> %d (dropped %d with no coverage)",
                  n_snps_before, nrow(self$snp_database$dp_matrix),
                  n_snps_before - nrow(self$snp_database$dp_matrix)))
    }
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
#' @title getCellsForSNPs: Extract Cell IDs Based on SNP Expression Criteria
#' @name getCellsForSNPs
#'
#' @description
#' Identifies cells that express specific SNPs at defined thresholds. This function
#' extracts cell IDs based on alternative allele fraction and depth criteria for
#' specified SNPs. Useful for annotating cells in Seurat or other single-cell objects
#' based on genetic variant expression patterns.
#'
#' @param snp_ids Character vector. SNP identifiers to query. Can be chromosome:position
#'   format (e.g., "1:12345") or rs IDs if available in the database.
#' @param min_alt_frac Numeric. Minimum alternative allele fraction threshold.
#'   Cells must have alt_frac >= this value to be included. Use 0 to include
#'   cells not expressing the alternative allele. Default: 0.2.
#' @param max_alt_frac Numeric. Maximum alternative allele fraction threshold.
#'   Cells must have alt_frac <= this value to be included. Use 1.0 to include
#'   all expressing cells. Set to 0 to get only reference allele cells. Default: 1.0.
#' @param min_dp Integer. Minimum depth (read coverage) required for the SNP
#'   in each cell. Cells with DP < min_dp are excluded. Default: 5.
#' @param sample_ids Character vector, optional. Restrict analysis to specific samples.
#'   If NULL, uses all samples in the database. Default: NULL.
#'
#' @return A named list where each element corresponds to a queried SNP:
#'   \item{snp_id}{Character vector of cell IDs meeting the criteria for this SNP}
#'   If a SNP is not found, returns an empty character vector with a warning.
#'   The list also includes an attribute "summary" with per-SNP statistics.
#'
#' @details
#' This function searches through the SNP database to find cells expressing
#' specified SNPs within defined thresholds. It's particularly useful for:
#' - Identifying cells carrying specific mutations
#' - Finding cells not expressing certain alleles (min_alt_frac = 0, max_alt_frac = 0)
#' - Annotating cell subsets based on genetic variants for downstream analysis
#' - Quality control based on read depth requirements
#'
#' The function handles both chromosome:position identifiers and rs IDs if the
#' database was built with rs ID annotation. SNPs are matched using exact string
#' matching on the SNP identifier.
#'
#' @note
#' - Cells must meet ALL criteria (alt_frac range AND min_dp) to be included
#' - The function only searches SNPs present in the current database
#' - For cells not expressing the alt allele, use min_alt_frac=0, max_alt_frac=0
#' - Empty results may indicate SNP not found or no cells meeting criteria
#'
#' @examples
#' \dontrun{
#' # Find cells expressing specific SNPs above 20% alt fraction
#' expressing_cells <- project$getCellsForSNPs(
#'   snp_ids = c("1:12345", "2:67890"),
#'   min_alt_frac = 0.2,
#'   min_dp = 5
#' )
#'
#' # Find cells NOT expressing the alternative allele (reference only)
#' ref_cells <- project$getCellsForSNPs(
#'   snp_ids = "1:12345",
#'   min_alt_frac = 0,
#'   max_alt_frac = 0,
#'   min_dp = 5
#' )
#'
#' # Find cells with high expression of alt allele (>80%)
#' high_alt_cells <- project$getCellsForSNPs(
#'   snp_ids = "rs123456",
#'   min_alt_frac = 0.8,
#'   min_dp = 10
#' )
#'
#' # Use results to annotate Seurat object
#' seurat_obj$expressing_snp1 <- ifelse(
#'   colnames(seurat_obj) %in% expressing_cells$`1:12345`,
#'   "expressing", "not_expressing"
#' )
#' }
#'
#' @seealso
#' \code{\link{findDESNPs}} for differential SNP analysis
#' \code{\link{plotSNPs}} for SNP visualization
variantCell$set("public", "getCellsForSNPs", function(snp_ids,
                                                      min_alt_frac = 0.2,
                                                      max_alt_frac = 1.0,
                                                      min_dp = 5,
                                                      sample_ids = NULL) {
  
  # Input validation
  if(is.null(self$snp_database) || is.null(self$snp_database$snp_info)) {
    stop("SNP database not built. Use buildSNPDatabase() first.")
  }
  
  if(!is.character(snp_ids) || length(snp_ids) == 0) {
    stop("snp_ids must be a non-empty character vector")
  }
  
  if(!is.numeric(min_alt_frac) || !is.numeric(max_alt_frac) || 
     !is.numeric(min_dp) || length(min_alt_frac) != 1 || 
     length(max_alt_frac) != 1 || length(min_dp) != 1) {
    stop("min_alt_frac, max_alt_frac, and min_dp must be single numeric values")
  }
  
  if(min_alt_frac < 0 || min_alt_frac > 1 || max_alt_frac < 0 || max_alt_frac > 1) {
    stop("alt_frac thresholds must be between 0 and 1")
  }
  
  if(min_alt_frac > max_alt_frac) {
    stop("min_alt_frac cannot be greater than max_alt_frac")
  }
  
  if(min_dp < 0) {
    stop("min_dp must be >= 0")
  }
  
  # Filter samples if specified
  if(!is.null(sample_ids)) {
    if(!all(sample_ids %in% names(self$samples))) {
      missing_samples <- setdiff(sample_ids, names(self$samples))
      stop(sprintf("Sample(s) not found: %s", paste(missing_samples, collapse=", ")))
    }
    available_cells <- self$snp_database$cell_metadata$sample_id %in% sample_ids
    cell_subset <- which(available_cells)
  } else {
    cell_subset <- 1:nrow(self$snp_database$cell_metadata)
  }
  
  # Initialize result list and summary
  result_list <- list()
  summary_stats <- data.frame(
    snp_id = character(0),
    snp_found = logical(0),
    total_cells_with_data = integer(0),
    cells_meeting_criteria = integer(0),
    mean_alt_frac = numeric(0),
    mean_dp = numeric(0),
    stringsAsFactors = FALSE
  )
  
  # Process each SNP
  for(snp_id in snp_ids) {
    cat(sprintf("Processing SNP: %s\n", snp_id))
    
    # Find SNP in database - check both exact match and if rs IDs exist
    snp_idx <- NULL
    if(snp_id %in% self$snp_database$snp_info$snp_id) {
      snp_idx <- which(self$snp_database$snp_info$snp_id == snp_id)
    } else if("rs_id" %in% colnames(self$snp_database$snp_info) &&
              snp_id %in% self$snp_database$snp_info$rs_id) {
      snp_idx <- which(self$snp_database$snp_info$rs_id == snp_id)
    }
    
    if(is.null(snp_idx) || length(snp_idx) == 0) {
      warning(sprintf("SNP '%s' not found in database", snp_id))
      result_list[[snp_id]] <- character(0)
      
      # Add to summary
      summary_stats <- rbind(summary_stats, data.frame(
        snp_id = snp_id,
        snp_found = FALSE,
        total_cells_with_data = 0,
        cells_meeting_criteria = 0,
        mean_alt_frac = NA,
        mean_dp = NA,
        stringsAsFactors = FALSE
      ))
      next
    }
    
    # Get AD and DP data for this SNP across specified cells
    ad_row <- self$snp_database$ad_matrix[snp_idx, cell_subset]
    dp_row <- self$snp_database$dp_matrix[snp_idx, cell_subset]
    
    # Calculate alt fractions
    alt_fractions <- rep(0, length(dp_row))
    valid_dp <- dp_row > 0
    alt_fractions[valid_dp] <- ad_row[valid_dp] / dp_row[valid_dp]
    
    # Apply filters
    depth_filter <- dp_row >= min_dp
    alt_frac_filter <- alt_fractions >= min_alt_frac & alt_fractions <= max_alt_frac
    combined_filter <- depth_filter & alt_frac_filter
    
    # Get qualifying cell IDs
    qualifying_cell_indices <- cell_subset[combined_filter]
    qualifying_cell_ids <- self$snp_database$cell_metadata$cell_id[qualifying_cell_indices]
    result_list[[snp_id]] <- qualifying_cell_ids
    
    # Calculate summary statistics for cells with data
    cells_with_data <- cell_subset[valid_dp]
    
    summary_stats <- rbind(summary_stats, data.frame(
      snp_id = snp_id,
      snp_found = TRUE,
      total_cells_with_data = length(cells_with_data),
      cells_meeting_criteria = length(qualifying_cell_ids),
      mean_alt_frac = ifelse(length(cells_with_data) > 0, 
                            mean(alt_fractions[valid_dp]), NA),
      mean_dp = ifelse(length(cells_with_data) > 0, 
                      mean(dp_row[valid_dp]), NA),
      stringsAsFactors = FALSE
    ))
    
    cat(sprintf("  Found: %d cells meeting criteria (%.1f%% alt_frac range, %d+ DP)\n",
                length(qualifying_cell_ids), min_alt_frac*100, min_dp))
  }
  
  # Add summary as attribute
  attr(result_list, "summary") <- summary_stats
  
  # Print overall summary
  cat(sprintf("\nSummary for %d SNP(s):\n", length(snp_ids)))
  for(i in 1:nrow(summary_stats)) {
    s <- summary_stats[i, ]
    if(s$snp_found) {
      cat(sprintf("  %s: %d/%d cells (%.1f%%)\n", 
                  s$snp_id, s$cells_meeting_criteria, s$total_cells_with_data,
                  ifelse(s$total_cells_with_data > 0, 
                        100 * s$cells_meeting_criteria / s$total_cells_with_data, 0)))
    } else {
      cat(sprintf("  %s: NOT FOUND\n", s$snp_id))
    }
  }
  
  return(result_list)
})
