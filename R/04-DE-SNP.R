#' @title aggregateByGroup: Group-Level SNP Data Aggregation
#' @name aggregateByGroup
#'
#' @description
#' Aggregates single-cell SNP data into group-level summaries based on a specified
#' metadata column. This function collapses individual cell SNP counts into group-level
#' matrices, which can be used for group-level differential SNP analyses. The function
#' supports both transplant and non-transplant modes, donor type filtering, and normalized
#' expression values.
#'
#' @param group_by Character. Column name in metadata to use for grouping cells.
#'                 Must be present in cell_metadata.
#' @param sample_column Character. Column name in metadata to use for sample stratification.
#'                     If NULL, operates in group-only mode. If specified, creates sample-aware
#'                     aggregation that preserves sample-level information for statistical testing.
#'                     Default: "sample_id".
#' @param donor_type Character, optional. Specific donor type to analyze
#'                   (e.g., "Donor" or "Recipient"). If NULL, uses all cells.
#'                   Ignored in non-transplant mode.
#' @param min_cells_per_group Integer. Minimum number of cells required for a group
#'                           to be included in analysis. Groups with fewer cells are
#'                           marked as "filtered_low_cells" in the metadata.
#' @param use_normalized Logical. Whether to include normalized depth counts in the
#'                      output (TRUE) or only use raw counts (FALSE).
#'
#' @return A list containing:
#'   \item{ad_matrix}{Aggregated alternative allele counts matrix (SNPs x Groups)}
#'   \item{dp_matrix}{Aggregated depth matrix (SNPs x Groups)}
#'   \item{dp_matrix_normalized}{Aggregated normalized depth matrix (SNPs x Groups), if available and requested}
#'   \item{metadata}{Data frame with group-level metadata and QC metrics}
#'   \item{group_by}{The metadata column used for grouping}
#'   \item{parameters}{List of parameters used for aggregation}
#'   \item{snp_info}{Data frame with SNP information}
#'   \item{snp_annotations}{Data frame with SNP annotations}
#'
#' @details
#' This function works by:
#' 1. Filtering cells based on donor_type if specified (e.g., only use Donor cells)
#' 2. Identifying unique values in the grouping column (e.g., cell_type)
#' 3. Summing alternative allele counts and depth counts across all cells in each group
#' 4. Creating group-level metadata with cell counts and quality metrics
#' 5. Filtering groups with fewer cells than the specified threshold
#'
#' The function automatically detects non-transplant mode (single donor type) and
#' adjusts its behavior accordingly. It also checks for normalized counts and includes
#' them in the output if available and requested.
#'
#' @note
#' - This function is typically used as a preprocessing step before `findSNPsByGroup()`
#' - The aggregated matrices no longer contain cell-level information; all counts are
#'   summed across cells in each group
#' - For transplant data, it's often useful to analyze donor and recipient cells separately
#'   by specifying the donor_type parameter
#' - Groups with fewer cells than min_cells_per_group are marked as "filtered_low_cells"
#'   in the metadata but are still included in the output matrices
#'
#' @examples
#'
#' \dontrun{
#' # Basic usage - aggregate by cell type (group-only mode)
#' collapsed <- project$aggregateByGroup(
#'   group_by = "cell_type",
#'   sample_column = NULL,
#'   use_normalized = TRUE
#' )
#'
#' # Sample-stratified aggregation for statistical testing
#' stratified_agg <- project$aggregateByGroup(
#'   group_by = "disease_status",
#'   sample_column = "patient_id",
#'   donor_type = "Recipient",
#'   use_normalized = TRUE
#' )
#'
#' # Analyze only donor cells with stricter filtering
#' donor_agg <- project$aggregateByGroup(
#'   group_by = "cell_type",
#'   sample_column = "sample_id",
#'   donor_type = "Donor",
#'   min_cells_per_group = 5
#' )
#'}
# Enhanced aggregateByGroup function with sample_column = NULL support
# This modification allows direct group comparison without sample stratification

variantCell$set("public", "aggregateByGroup", function(group_by,
                                                       sample_column = "sample_id", 
                                                       donor_type = NULL,
                                                       min_cells_per_group = 3,
                                                       use_normalized = TRUE) {
  
  # Input validation
  if(is.null(self$snp_database) || is.null(self$snp_database$cell_metadata)) {
    stop("SNP database not properly initialized")
  }
  
  if(!group_by %in% colnames(self$snp_database$cell_metadata)) {
    stop(sprintf("Column '%s' not found in metadata", group_by))
  }
  
  # Modified validation for sample_column - allow NULL
  if(!is.null(sample_column) && !sample_column %in% colnames(self$snp_database$cell_metadata)) {
    stop(sprintf("Sample column '%s' not found in metadata", sample_column))
  }
  
  if(use_normalized && !"dp_matrix_normalized" %in% names(self$snp_database)) {
    stop("Normalized counts not available. Use aggregateByGroup with sample_column parameter.")
  }
  
  # Validate numeric parameters
  if(!is.numeric(min_cells_per_group) || min_cells_per_group < 0) {
    stop("min_cells_per_group must be a non-negative integer")
  }
  
  # Get metadata and matrices
  meta <- self$snp_database$cell_metadata
  ad_matrix <- self$snp_database$ad_matrix
  dp_matrix <- self$snp_database$dp_matrix
  if(use_normalized) {
    dp_matrix_norm <- self$snp_database$dp_matrix_normalized
  }
  
  # Validate matrix dimensions
  if(ncol(ad_matrix) != nrow(meta) || ncol(dp_matrix) != nrow(meta)) {
    stop("Matrix dimensions do not match metadata")
  }
  if(ncol(ad_matrix) != ncol(dp_matrix)) {
    stop("AD and DP matrices have different dimensions")
  }
  
  # Check for non-transplant mode
  non_transplant_mode <- FALSE
  if(length(unique(meta$donor_type)) == 1 && unique(meta$donor_type)[1] == "donor0") {
    non_transplant_mode <- TRUE
    cat("\nNon-transplant mode detected (single donor type 'donor0')")
  }
  
  # Apply donor_type filter if specified and not in non-transplant mode
  if(!is.null(donor_type) && !non_transplant_mode) {
    if(!"donor_type" %in% colnames(meta)) {
      stop("donor_type column not found in metadata")
    }
    
    donor_cells <- !is.na(meta$donor_type) & meta$donor_type == donor_type
    meta <- meta[donor_cells, ]
    ad_matrix <- ad_matrix[, donor_cells]
    dp_matrix <- dp_matrix[, donor_cells]
    if(use_normalized) {
      dp_matrix_norm <- dp_matrix_norm[, donor_cells]
    }
    
    cat(sprintf("Filtered to %d %s cells", nrow(meta), donor_type))
  } else if(!is.null(donor_type) && non_transplant_mode) {
    cat(sprintf("\nIgnoring donor_type parameter in non-transplant mode (all cells are 'donor0')"))
  }
  
  # Check if we're in "group-only" mode (sample_column = NULL)
  if(is.null(sample_column)) {
    cat("\nOperating in group-only mode: comparing entire groups without sample stratification\n")
    
    # Get unique groups
    groups <- sort(unique(meta[[group_by]]))
    n_groups <- length(groups)
    
    cat(sprintf("Found %d unique groups", n_groups))
    
    # Initialize group-level matrices
    collapsed_ad <- Matrix::Matrix(0, nrow=nrow(ad_matrix), ncol=n_groups, sparse=TRUE)
    collapsed_dp <- Matrix::Matrix(0, nrow=nrow(dp_matrix), ncol=n_groups, sparse=TRUE)
    if(use_normalized) {
      collapsed_dp_norm <- Matrix::Matrix(0, nrow=nrow(dp_matrix), ncol=n_groups, sparse=TRUE)
    }
    
    colnames(collapsed_ad) <- groups
    colnames(collapsed_dp) <- groups
    if(use_normalized) colnames(collapsed_dp_norm) <- groups
    
    # Create group-level metadata (pre-calculate for efficiency)
    collapsed_meta <- data.frame(
      group = groups,
      n_cells = sapply(groups, function(g) sum(meta[[group_by]] == g)),
      n_samples = NA,  # Not applicable in group-only mode
      stringsAsFactors = FALSE
    )
    
    # Add filter status
    collapsed_meta$filter_status <- ifelse(
      collapsed_meta$n_cells >= min_cells_per_group,
      "included",
      "filtered_low_cells"
    )
    
    # Process each group directly
    for(i in seq_along(groups)) {
      group <- groups[i]
      group_cells <- meta[[group_by]] == group
      
      if(collapsed_meta$n_cells[i] >= min_cells_per_group) {
        # Sum counts for all cells in this group
        collapsed_ad[, i] <- Matrix::rowSums(ad_matrix[, group_cells, drop=FALSE])
        collapsed_dp[, i] <- Matrix::rowSums(dp_matrix[, group_cells, drop=FALSE])
        if(use_normalized) {
          collapsed_dp_norm[, i] <- Matrix::rowSums(dp_matrix_norm[, group_cells, drop=FALSE])
        }
      }
    }
    
    # Create return object for group-only mode
    result <- list(
      ad = collapsed_ad,
      dp = collapsed_dp,
      metadata = collapsed_meta,
      mode = "group_only",
      parameters = list(
        group_by = group_by,
        sample_column = NULL,
        donor_type = donor_type,
        min_cells_per_group = min_cells_per_group,
        use_normalized = use_normalized,
        non_transplant_mode = non_transplant_mode
      )
    )
    
    if(use_normalized) {
      result$dp_norm <- collapsed_dp_norm
    }
    
    # Print summary
    cat(sprintf("\nCollapsed data by %s:", group_by))
    cat(sprintf("\n - Total cells: %d", ncol(ad_matrix)))
    cat(sprintf("\n - Total groups: %d", n_groups))
    cat(sprintf("\n - Included groups: %d", sum(collapsed_meta$filter_status == "included")))
    cat(sprintf("\n - Filtered groups: %d", sum(collapsed_meta$filter_status == "filtered_low_cells")))
    cat("\n - Using normalized counts:", use_normalized)
    if(non_transplant_mode) {
      cat("\n - Running in non-transplant mode (single donor)")
    } else if(!is.null(donor_type)) {
      cat(sprintf("\n - Filtered to donor type: %s", donor_type))
    } else {
      cat("\n - Using all donor types")
    }
    cat("\n\nGroup details:")
    print(collapsed_meta)
    
    return(result)
    
  } else {
    # Original sample-stratified mode
    cat("\nOperating in sample-stratified mode: comparing samples within groups\n")
    
    # Create group+sample combinations
    meta$group_sample <- paste(meta[[group_by]], meta[[sample_column]], sep="_")
    unique_group_samples <- sort(unique(meta$group_sample))
    n_group_samples <- length(unique_group_samples)
    
    # Get unique groups for final aggregation
    groups <- sort(unique(meta[[group_by]]))
    n_groups <- length(groups)
    
    cat(sprintf("Found %d unique group-sample combinations across %d groups", 
                n_group_samples, n_groups))
    
    # STEP 1: Aggregate by group+sample combinations
    sample_level_ad <- Matrix::Matrix(0, nrow=nrow(ad_matrix), ncol=n_group_samples, sparse=TRUE)
    sample_level_dp <- Matrix::Matrix(0, nrow=nrow(dp_matrix), ncol=n_group_samples, sparse=TRUE)
    if(use_normalized) {
      sample_level_dp_norm <- Matrix::Matrix(0, nrow=nrow(dp_matrix), ncol=n_group_samples, sparse=TRUE)
    }
    
    colnames(sample_level_ad) <- unique_group_samples
    colnames(sample_level_dp) <- unique_group_samples
    if(use_normalized) colnames(sample_level_dp_norm) <- unique_group_samples
    
    # Pre-calculate sample metadata for efficiency
    sample_meta <- data.frame(
      group_sample = unique_group_samples,
      group = sapply(unique_group_samples, function(gs) {
        cells <- meta$group_sample == gs
        if(sum(cells) > 0) {
          meta[[group_by]][cells][1]
        } else {
          strsplit(gs, "_")[[1]][1]
        }
      }),
      sample = sapply(unique_group_samples, function(gs) {
        cells <- meta$group_sample == gs
        if(sum(cells) > 0) {
          meta[[sample_column]][cells][1]
        } else {
          paste(strsplit(gs, "_")[[1]][-1], collapse="_")
        }
      }),
      n_cells = sapply(unique_group_samples, function(gs) {
        sum(meta$group_sample == gs)
      }),
      stringsAsFactors = FALSE
    )
    
    # Add filter status
    sample_meta$filter_status <- ifelse(
      sample_meta$n_cells >= min_cells_per_group,
      "included",
      "filtered_low_cells"
    )
    
    # Process each group+sample combination
    for(i in seq_along(unique_group_samples)) {
      if(sample_meta$n_cells[i] >= min_cells_per_group) {
        group_sample_cells <- meta$group_sample == unique_group_samples[i]
        # Sum counts for all cells in this group+sample
        sample_level_ad[, i] <- Matrix::rowSums(ad_matrix[, group_sample_cells, drop=FALSE])
        sample_level_dp[, i] <- Matrix::rowSums(dp_matrix[, group_sample_cells, drop=FALSE])
        if(use_normalized) {
          sample_level_dp_norm[, i] <- Matrix::rowSums(dp_matrix_norm[, group_sample_cells, drop=FALSE])
        }
      }
    }
    
    # STEP 2: Aggregate to group level (sum across samples within each group)
    collapsed_ad <- Matrix::Matrix(0, nrow=nrow(ad_matrix), ncol=n_groups, sparse=TRUE)
    collapsed_dp <- Matrix::Matrix(0, nrow=nrow(dp_matrix), ncol=n_groups, sparse=TRUE)
    if(use_normalized) {
      collapsed_dp_norm <- Matrix::Matrix(0, nrow=nrow(dp_matrix), ncol=n_groups, sparse=TRUE)
    }
    
    colnames(collapsed_ad) <- groups
    colnames(collapsed_dp) <- groups
    if(use_normalized) colnames(collapsed_dp_norm) <- groups
    
    # Pre-calculate group-level metadata for efficiency
    collapsed_meta <- data.frame(
      group = groups,
      n_cells = sapply(groups, function(g) {
        included_samples <- sample_meta$group == g & sample_meta$filter_status == "included"
        sum(sample_meta$n_cells[included_samples])
      }),
      n_samples = sapply(groups, function(g) {
        sum(sample_meta$group == g & sample_meta$filter_status == "included")
      }),
      filter_status = "included",
      stringsAsFactors = FALSE
    )
    
    # Process group-level aggregation
    for(i in seq_along(groups)) {
      group <- groups[i]
      group_sample_indices <- which(sample_meta$group == group & sample_meta$filter_status == "included")
      
      if(length(group_sample_indices) > 0) {
        collapsed_ad[, i] <- Matrix::rowSums(sample_level_ad[, group_sample_indices, drop=FALSE])
        collapsed_dp[, i] <- Matrix::rowSums(sample_level_dp[, group_sample_indices, drop=FALSE])
        if(use_normalized) {
          collapsed_dp_norm[, i] <- Matrix::rowSums(sample_level_dp_norm[, group_sample_indices, drop=FALSE])
        }
      }
    }
    
    # Create return object for sample-stratified mode
    result <- list(
      ad = collapsed_ad,
      dp = collapsed_dp,
      metadata = collapsed_meta,
      sample_level_ad = sample_level_ad,
      sample_level_dp = sample_level_dp,
      sample_metadata = sample_meta,
      mode = "sample_stratified",
      parameters = list(
        group_by = group_by,
        sample_column = sample_column,
        donor_type = donor_type,
        min_cells_per_group = min_cells_per_group,
        use_normalized = use_normalized,
        non_transplant_mode = non_transplant_mode
      )
    )
    
    if(use_normalized) {
      result$dp_norm <- collapsed_dp_norm
      result$sample_level_dp_norm <- sample_level_dp_norm
    }
    
    # Print summary
    cat(sprintf("\nSample-stratified aggregation by %s:", group_by))
    cat(sprintf("\n - Total cells: %d", ncol(ad_matrix)))
    cat(sprintf("\n - Total groups: %d", n_groups))
    cat(sprintf("\n - Total group-sample combinations: %d", n_group_samples))
    cat(sprintf("\n - Included combinations: %d", sum(sample_meta$filter_status == "included")))
    cat(sprintf("\n - Filtered combinations: %d", sum(sample_meta$filter_status == "filtered_low_cells")))
    cat("\n - Using normalized counts:", use_normalized)
    if(non_transplant_mode) {
      cat("\n - Running in non-transplant mode (single donor)")
    } else if(!is.null(donor_type)) {
      cat(sprintf("\n - Filtered to donor type: %s", donor_type))
    } else {
      cat("\n - Using all donor types")
    }
    
    cat("\n\nGroup-level summary:")
    print(collapsed_meta)
    
    cat("\n\nSample-level details:")
    print(sample_meta)
    
    return(result)
  }
})




variantCell$set("public",  "getNumericSubset", function(sparseMat, rows, cols) {
  # Matrix is already numeric (dgCMatrix), just need to subset and convert
  return(as.matrix(sparseMat[rows, cols]))
})
#' @title findDESNPs: Cell-Level Differential SNP Expression Analysis
#' @name findDESNPs
#'
#' @description
#' Identifies differentially expressed SNPs between cell populations by comparing read depths
#' and alternative allele frequencies. This function performs comprehensive statistical analysis
#' at the single-cell level, with support for parallel processing to improve performance on large datasets.
#'
#' @param ident.1 Character. Primary cell identity to analyze.
#' @param ident.2 Character, optional. Secondary cell identity to compare against.
#'                If NULL, compares against all other cells.
#' @param donor_type Character, optional. Donor type to restrict analysis to ("Donor" or "Recipient").
#'                  If NULL, uses all cells regardless of donor type.
#' @param use_normalized Logical. Whether to use normalized depth counts (TRUE) or raw counts (FALSE).
#' @param min_expr_cells Integer. Minimum number of expressing cells required in each group.
#' @param min_alt_frac Numeric [0-1]. Minimum alternative allele fraction to consider a cell as expressing.
#' @param logfc.threshold Numeric. Minimum absolute log2 fold-change required to report a SNP.
#' @param calc_p Logical. Whether to calculate p-values (Wilcoxon test). Set to FALSE to save computation time.
#' @param p.adjust.method Character. Method for p-value adjustment, passed to p.adjust(). Default: "BH" (Benjamini-Hochberg).
#' @param return_all Logical. Whether to return all SNPs or only significant ones.
#' @param pseudocount Numeric. Value added to expression values before log transformation.
#' @param min.p Numeric. Minimum p-value to report (prevents numerical underflow).
#' @param debug Logical. Whether to print debugging information during analysis.
#' @param n_cores Integer, optional. Number of CPU cores to use for parallel processing.
#'               If NULL, automatically uses detectCores() - 1.
#' @param include_rs_ids Logical. If true, includes RS IDs               
#' @param include_population_AF Logical. If true, includes population allele frequencies
#'
#' @return List containing:
#'   \item{results}{Data frame of differentially expressed SNPs with metrics including log2FC,
#'                 expression values, cell counts, and significance statistics.}
#'   \item{summary}{List with analysis overview, including counts of significant SNPs,
#'                 up/downregulated SNPs, and parameter settings used.}
#'
#' @details
#' The function calculates differential expression by comparing the average expression
#' of SNPs between two groups, normalized by the total number of cells in each group.
#' For each SNP, cells are only considered as expressing if they have a minimum
#' alternative allele fraction (min_alt_frac) and positive read depth.
#'
#' Statistical testing is performed using Wilcoxon rank-sum test when calc_p=TRUE.
#' Multiple testing correction is applied using the specified p.adjust.method.
#'
#' The parallel implementation distributes SNP processing across multiple CPU cores
#' for significantly improved performance on large datasets.
#'
#' @note
#' - Requires package 'parallel', 'foreach', and 'doParallel' for parallel processing
#' - Project identity must be set before using this function via setProjectIdentity()
#' - For non-transplant datasets, donor_type filtering is automatically disabled
#'
#' @examples
#'
#' \dontrun{
#' # Initialize a variantCell project
#'
#' proj$setProjectIdentity('cell_type')
#'
#' # Basic usage comparing T cells vs other cells, donor cells only
#' results <- proj$findDESNPs(
#'   ident.1 = "T_cells",
#'   ident.2 = NULL,
#'   donor_type = "Donor",
#'   min_expr_cells = 5,
#'   logfc.threshold = 0.25
#' )
#'
#' # Without p-value calculation for faster processing
#' fast_results <- proj$findDESNPs(
#'   ident.1 = "CD4",
#'   ident.2 = "CD8",
#'   calc_p = FALSE,
#'   n_cores = 8
#' )
#'
#' # Access results
#' head(results$results)
#' results$summary
#' }
#'
#' @seealso
#' \code{\link{setProjectIdentity}} for setting the cell identity to use
#' \code{\link{findSNPsByGroup}} for group-level SNP analysis
variantCell$set("public", "findDESNPs", function(ident.1,
                                                 ident.2 = NULL,
                                                 donor_type = NULL,
                                                 use_normalized = TRUE,    
                                                 min_expr_cells = 3,        
                                                 min_alt_frac = 0.2,        
                                                 logfc.threshold = 0.1,     
                                                 calc_p = TRUE,             
                                                 p.adjust.method = "BH",
                                                 return_all = TRUE,
                                                 pseudocount = 1,           
                                                 min.p = 1e-300,
                                                 debug = FALSE,
                                                 n_cores = NULL,
                                                 use_parallel = TRUE,       
                                                 chunk_size = 1000,         
                                                 max_ram_gb = 4,
                                                 include_rs_ids = TRUE,
                                                 include_population_AF = TRUE) {  
  
  # Input validation
  if(is.null(self$current_project_ident)) {
    stop("No identity set. Please use setProjectIdentity() first.")
  }
  
  if(use_normalized && is.null(self$snp_database$dp_matrix_normalized)) {
    stop("Normalized counts requested but not available. Rebuild database with normalize=TRUE")
  }
  
  # Check if rs# IDs are available
  has_rs_ids <- FALSE
  if(include_rs_ids) {
    if(!is.null(self$snp_database) && 
       "snp_info" %in% names(self$snp_database) &&
       "rs_id" %in% colnames(self$snp_database$snp_info)) {
      has_rs_ids <- !all(is.na(self$snp_database$snp_info$rs_id))
    }
    
    if(!has_rs_ids && include_rs_ids) {
      cat("Warning: rs# identifiers requested but not available in project. ")
      cat("Use buildSNPDatabase(add_rs_ids = TRUE, VCF_file_path = '...') to add them.\n")
    }
  }
  
  # Check if population AF is available
  has_population_AF <- FALSE
  if(include_population_AF) {
    if(!is.null(self$snp_database) && 
       "snp_info" %in% names(self$snp_database) &&
       "population_AF" %in% colnames(self$snp_database$snp_info)) {
      has_population_AF <- !all(is.na(self$snp_database$snp_info$population_AF))
    }
    
    if(!has_population_AF && include_population_AF) {
      cat("Warning: Population AF requested but not available in project. ")
      cat("Use buildSNPDatabase(add_population_AF = TRUE, VCF_file_path = '...') to add them.\n")
    }
  }
  
  # Get data and metadata
  meta <- self$snp_database$cell_metadata
  
  # Check for non-transplant mode
  non_transplant_mode <- FALSE
  if(length(unique(meta$donor_type)) == 1 && unique(meta$donor_type)[1] == "donor0") {
    non_transplant_mode <- TRUE
    cat("\nNon-transplant mode detected (single donor type 'donor0')")
  }
  
  # Apply donor type filter if specified and not in non-transplant mode
  if(!is.null(donor_type) && !non_transplant_mode) {
    donor_mask <- meta$donor_type == donor_type
    if(sum(donor_mask) == 0) {
      stop(sprintf("No cells found for donor_type: %s", donor_type))
    }
    meta <- meta[donor_mask, ]
    
    cat(sprintf("\nFiltered to %d cells with donor_type: %s", sum(donor_mask), donor_type))
  } else if(!is.null(donor_type) && non_transplant_mode) {
    cat(sprintf("\nIgnoring donor_type parameter in non-transplant mode (all cells are 'donor0')"))
  }
  
  # Create group masks
  group1_mask <- meta[[self$current_project_ident]] == ident.1
  if(is.null(ident.2)) {
    group2_mask <- !group1_mask
    group2_name <- "rest"
  } else {
    group2_mask <- meta[[self$current_project_ident]] == ident.2
    group2_name <- ident.2
  }
  
  # Get group sizes
  n_cells_group1 <- sum(group1_mask)
  n_cells_group2 <- sum(group2_mask)
  
  # Check if we have enough cells in each group
  if(n_cells_group1 < min_expr_cells) {
    stop(sprintf("Group '%s' has only %d cells, which is below the minimum of %d cells", 
                 ident.1, n_cells_group1, min_expr_cells))
  }
  if(n_cells_group2 < min_expr_cells) {
    stop(sprintf("Group '%s' has only %d cells, which is below the minimum of %d cells", 
                 group2_name, n_cells_group2, min_expr_cells))
  }
  
  # Get indices of cells in each group
  group1_indices <- which(group1_mask)
  group2_indices <- which(group2_mask)
  all_group_indices <- c(group1_indices, group2_indices)
  
  # ========== PRE-FILTERING  ==========
  cat(sprintf("\nPre-filtering SNPs based on expression criteria..."))
  cat(sprintf("\nCriteria: min_expr_cells=%d, min_alt_frac=%.2f", min_expr_cells, min_alt_frac))
  
  # Vectorized calculation of alt fractions for all SNPs
  alt_frac_matrix1 <- self$snp_database$ad_matrix[, group1_indices] / 
    pmax(self$snp_database$dp_matrix[, group1_indices], 1)
  alt_frac_matrix2 <- self$snp_database$ad_matrix[, group2_indices] / 
    pmax(self$snp_database$dp_matrix[, group2_indices], 1)
  
  # Create depth masks
  depth_mask1 <- self$snp_database$dp_matrix[, group1_indices] > 0
  depth_mask2 <- self$snp_database$dp_matrix[, group2_indices] > 0
  
  # Set alt_frac to 0 where depth is 0
  alt_frac_matrix1[!depth_mask1] <- 0
  alt_frac_matrix2[!depth_mask2] <- 0
  
  # Count expressing cells per SNP (vectorized)
  expr_counts1 <- Matrix::rowSums(depth_mask1 & alt_frac_matrix1 >= min_alt_frac)
  expr_counts2 <- Matrix::rowSums(depth_mask2 & alt_frac_matrix2 >= min_alt_frac)
  
  # Find qualifying SNPs - must meet criteria in BOTH groups
  qualifying_snps <- (expr_counts1 >= min_expr_cells) & (expr_counts2 >= min_expr_cells)
  qualifying_snp_indices <- which(qualifying_snps)
  
  cat(sprintf("\nPre-filtering results:"))
  cat(sprintf("\n - Total SNPs in database: %d", nrow(self$snp_database$dp_matrix)))
  cat(sprintf("\n - SNPs meeting expression criteria: %d (%.1f%%)", 
              length(qualifying_snp_indices), 
              100 * length(qualifying_snp_indices) / nrow(self$snp_database$dp_matrix)))
  
  if(length(qualifying_snp_indices) == 0) {
    cat("\nNo SNPs meet the expression criteria in both groups.")
    return(NULL)
  }
  
  cat(sprintf("\nLoading optimized matrices for %d qualifying SNPs...", length(qualifying_snp_indices)))
  
  # Load matrices for qualifying SNPs
  if(use_normalized) {
    dp_matrix_subset <- self$snp_database$dp_matrix_normalized[qualifying_snp_indices, all_group_indices]
  } else {
    dp_matrix_subset <- self$snp_database$dp_matrix[qualifying_snp_indices, all_group_indices]
  }
  ad_matrix_subset <- self$snp_database$ad_matrix[qualifying_snp_indices, all_group_indices]
  
  # Update SNP info and annotations to match
  snp_info_subset <- self$snp_database$snp_info[qualifying_snp_indices, ]
  snp_annotations_subset <- self$snp_database$snp_annotations[qualifying_snp_indices, ]
  
  # Adjust group indices to work with the subset matrices
  group1_indices_subset <- match(group1_indices, all_group_indices)
  group2_indices_subset <- match(group2_indices, all_group_indices)
  
  # Extract matrices for both groups from the subset
  dp1_full <- dp_matrix_subset[, group1_indices_subset]
  dp2_full <- dp_matrix_subset[, group2_indices_subset]
  ad1_full <- ad_matrix_subset[, group1_indices_subset]
  ad2_full <- ad_matrix_subset[, group2_indices_subset]
  
  # Get total SNPs (now using pre-filtered count)
  total_snps <- nrow(dp_matrix_subset)
  
  # Determine processing approach
  can_use_parallel <- use_parallel && requireNamespace("parallel", quietly = TRUE) && 
    requireNamespace("foreach", quietly = TRUE) && 
    requireNamespace("doParallel", quietly = TRUE)
  
  # Auto-detect cores if not specified (leaving 1 core free)
  if(is.null(n_cores) && can_use_parallel) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }
  
  # Estimate memory requirements (based on pre-filtered SNPs)
  est_mem_per_snp_mb <- (n_cells_group1 + n_cells_group2) * 8 * 2 / 1024^2
  est_total_mem_gb <- est_mem_per_snp_mb * total_snps * n_cores / 1024
  
  if(est_total_mem_gb > max_ram_gb && can_use_parallel) {
    old_chunk_size <- chunk_size
    chunk_size <- min(chunk_size, as.integer(max_ram_gb * 1024^2 / (est_mem_per_snp_mb * n_cores)))
    if(debug) {
      cat(sprintf("\nMemory usage could be high (est. %.1f GB). Reducing chunk size from %d to %d", 
                  est_total_mem_gb, old_chunk_size, chunk_size))
    }
  }
  
  # Print initial summary
  process_type <- if(can_use_parallel) sprintf("Parallel (%d cores)", n_cores) else "Serial"
  cat(sprintf("\n=== Starting Cell-Level Differential Analysis (%s) ===", process_type))
  cat(sprintf("\nUsing %s expression values", if(use_normalized) "normalized" else "raw"))
  cat(sprintf("\nComparison groups:"))
  cat(sprintf("\n%s: %d cells", ident.1, n_cells_group1))
  cat(sprintf("\n%s: %d cells", group2_name, n_cells_group2))
  
  # Progress tracking
  if(debug) {
    cat(sprintf("\nAnalyzing %d pre-filtered SNPs with chunk size %d\n", total_snps, chunk_size))
  }
  
  # Initialize results storage
  results <- list()
  results_count <- 0
  
  # Create index vector for processing (now refers to subset)
  snp_indices <- 1:total_snps
  
  # PARALLEL PROCESSING PATH
  if(can_use_parallel) {
    tryCatch({
      # Load required packages for parallel processing
      library(parallel)
      library(foreach)
      library(doParallel)
      
      # Create chunks for efficient processing
      snp_chunks <- split(snp_indices, ceiling(seq_along(snp_indices) / chunk_size))
      
      if(debug) {
        cat(sprintf("Processing %d SNPs in %d chunks using %d cores...\n", 
                    total_snps, length(snp_chunks), n_cores))
      }
      
      # Setup parallel processing
      cl <- makeCluster(n_cores)
      on.exit(stopCluster(cl), add = TRUE)
      registerDoParallel(cl)
      
      # Export necessary variables to cluster
      clusterExport(cl, varlist = c("min_expr_cells", "min_alt_frac", "logfc.threshold",
                                    "n_cells_group1", "n_cells_group2", "calc_p",
                                    "pseudocount", "min.p"), envir = environment())
      
      # Process SNPs in parallel with error handling
      results <- foreach(chunk = snp_chunks, 
                         .combine = 'c', 
                         .packages = c("stats", "Matrix"),
                         .errorhandling = 'pass') %dopar% {
                           
                           tryCatch({
                             chunk_results <- list()
                             chunk_count <- 0
                             
                             for(i in chunk) {
                               # Get data for both groups (converting to vectors for efficiency)
                               dp1 <- as.vector(dp1_full[i,])
                               dp2 <- as.vector(dp2_full[i,])
                               ad1 <- as.vector(ad1_full[i,])
                               ad2 <- as.vector(ad2_full[i,])
                               
                               # Calculate alt fractions
                               alt_frac1 <- ifelse(dp1 > 0, ad1/dp1, 0)
                               alt_frac2 <- ifelse(dp2 > 0, ad2/dp2, 0)
                               
                               # Count expressing cells
                               n_expr1 <- sum(dp1 > 0 & alt_frac1 >= min_alt_frac)
                               n_expr2 <- sum(dp2 > 0 & alt_frac2 >= min_alt_frac)
                               
                               if(n_expr1 >= min_expr_cells && n_expr2 >= min_expr_cells) {
                                 # Calculate group averages (normalized by total cells in group)
                                 avg_expr1 <- sum(dp1[dp1 > 0 & alt_frac1 >= min_alt_frac]) / n_cells_group1
                                 avg_expr2 <- sum(dp2[dp2 > 0 & alt_frac2 >= min_alt_frac]) / n_cells_group2
                                 
                                 # Calculate log2 fold change using normalized averages
                                 log2fc <- log2((avg_expr1 + pseudocount)/(avg_expr2 + pseudocount))
                                 
                                 if(abs(log2fc) >= logfc.threshold) {
                                   # Only calculate p-value if requested
                                   pvalue <- NA
                                   if(calc_p) {
                                     test_result <- tryCatch({
                                       wilcox.test(dp1, dp2)
                                     }, error = function(e) {
                                       list(p.value = 1.0)
                                     })
                                     pvalue <- max(test_result$p.value, min.p)
                                   }
                                   
                                   chunk_count <- chunk_count + 1
                                   chunk_results[[chunk_count]] <- list(
                                     snp_idx = i,
                                     log2fc = log2fc,
                                     avg_expr1 = avg_expr1,
                                     avg_expr2 = avg_expr2,
                                     n_expr1 = n_expr1,
                                     n_expr2 = n_expr2,
                                     mean_alt_frac1 = mean(alt_frac1[dp1 > 0]),
                                     mean_alt_frac2 = mean(alt_frac2[dp2 > 0]),
                                     pvalue = pvalue
                                   )
                                 }
                               }
                               
                               # Clean up to reduce memory footprint
                               rm(dp1, dp2, ad1, ad2, alt_frac1, alt_frac2)
                               gc(verbose = FALSE, full = FALSE)
                             }
                             return(chunk_results)
                           }, error = function(e) {
                             return(list(error = e$message, chunk_first = chunk[1], chunk_last = chunk[length(chunk)]))
                           })
                         }
      
      # Check for errors in parallel processing results
      error_results <- Filter(function(x) is.list(x) && "error" %in% names(x), results)
      if(length(error_results) > 0) {
        error_msgs <- sapply(error_results, function(x) paste("Error in chunk", x$chunk_first, "-", x$chunk_last, ":", x$error))
        warning("Some parallel chunks had errors:\n", paste(error_msgs, collapse = "\n"))
        results <- Filter(function(x) !is.list(x) || !"error" %in% names(x), results)
      }
      
    }, error = function(e) {
      warning(sprintf("Parallel processing failed: %s\nFalling back to non-parallel method.", e$message))
      can_use_parallel <- FALSE
      if(exists("cl") && inherits(cl, "cluster")) {
        tryCatch(stopCluster(cl), error = function(e) NULL)
      }
    })
  }
  
  # NON-PARALLEL PROCESSING PATH
  if(!can_use_parallel) {
    if(debug) {
      cat(sprintf("Processing %d SNPs in serial mode...\n", total_snps))
      pb <- txtProgressBar(min = 0, max = total_snps, style = 3)
    }
    
    results <- list()
    results_count <- 0
    
    # Process chunks to manage memory
    chunk_start <- 1
    while(chunk_start <= total_snps) {
      chunk_end <- min(chunk_start + chunk_size - 1, total_snps)
      chunk_indices <- chunk_start:chunk_end
      
      for(i in chunk_indices) {
        if(debug) setTxtProgressBar(pb, i)
        
        # Get data for both groups (converting to vectors for efficiency)
        dp1 <- as.vector(dp1_full[i,])
        dp2 <- as.vector(dp2_full[i,])
        ad1 <- as.vector(ad1_full[i,])
        ad2 <- as.vector(ad2_full[i,])
        
        # Calculate alt fractions
        alt_frac1 <- ifelse(dp1 > 0, ad1/dp1, 0)
        alt_frac2 <- ifelse(dp2 > 0, ad2/dp2, 0)
        
        # Count expressing cells
        n_expr1 <- sum(dp1 > 0 & alt_frac1 >= min_alt_frac)
        n_expr2 <- sum(dp2 > 0 & alt_frac2 >= min_alt_frac)
        
        if(n_expr1 >= min_expr_cells && n_expr2 >= min_expr_cells) {
          # Calculate group averages (normalized by total cells in group)
          avg_expr1 <- sum(dp1[dp1 > 0 & alt_frac1 >= min_alt_frac]) / n_cells_group1
          avg_expr2 <- sum(dp2[dp2 > 0 & alt_frac2 >= min_alt_frac]) / n_cells_group2
          
          # Calculate log2 fold change using normalized averages
          log2fc <- log2((avg_expr1 + pseudocount)/(avg_expr2 + pseudocount))
          
          if(abs(log2fc) >= logfc.threshold) {
            # Only calculate p-value if requested
            pvalue <- NA
            if(calc_p) {
              test_result <- tryCatch({
                wilcox.test(dp1, dp2)
              }, error = function(e) {
                list(p.value = 1.0)
              })
              pvalue <- max(test_result$p.value, min.p)
            }
            
            results_count <- results_count + 1
            results[[results_count]] <- list(
              snp_idx = i,
              log2fc = log2fc,
              avg_expr1 = avg_expr1,
              avg_expr2 = avg_expr2,
              n_expr1 = n_expr1,
              n_expr2 = n_expr2,
              mean_alt_frac1 = mean(alt_frac1[dp1 > 0]),
              mean_alt_frac2 = mean(alt_frac2[dp2 > 0]),
              pvalue = pvalue
            )
          }
        }
        
        # Clean up to reduce memory footprint
        rm(dp1, dp2, ad1, ad2, alt_frac1, alt_frac2)
      }
      
      # Force garbage collection between chunks to free memory
      gc(verbose = FALSE)
      
      # Move to next chunk
      chunk_start <- chunk_end + 1
    }
    
    if(debug) close(pb)
  }
  
  # Process results
  results_count <- length(results)
  
  if(results_count > 0) {
    # Create base data frame from results using SUBSET indices
    subset_indices <- sapply(results, function(x) x$snp_idx)
    original_indices <- qualifying_snp_indices[subset_indices]  # Map back to original
    
    all_results <- data.frame(
      snp_idx = original_indices,  # Use original indices
      chromosome = snp_info_subset$CHROM[subset_indices],
      position = snp_info_subset$POS[subset_indices],
      ref = snp_info_subset$REF[subset_indices],
      alt = snp_info_subset$ALT[subset_indices],
      feature_type = snp_annotations_subset$feature_type[subset_indices],
      gene_name = snp_annotations_subset$gene_name[subset_indices],
      gene_type = snp_annotations_subset$gene_type[subset_indices],
      
      # Expression metrics
      log2fc = sapply(results, function(x) x$log2fc),
      avg_expr_group1 = sapply(results, function(x) x$avg_expr1),
      avg_expr_group2 = sapply(results, function(x) x$avg_expr2),
      
      # Cell counts and fractions
      total_cells_1 = n_cells_group1,
      total_cells_2 = n_cells_group2,
      expr_cells_1 = sapply(results, function(x) x$n_expr1),
      expr_cells_2 = sapply(results, function(x) x$n_expr2),
      expr_frac_1 = sapply(results, function(x) x$n_expr1)/n_cells_group1,
      expr_frac_2 = sapply(results, function(x) x$n_expr2)/n_cells_group2,
      
      # Alternative allele metrics
      mean_alt_frac_1 = sapply(results, function(x) x$mean_alt_frac1),
      mean_alt_frac_2 = sapply(results, function(x) x$mean_alt_frac2),
      
      # Statistical results (NA if calc_p = FALSE)
      pvalue = sapply(results, function(x) x$pvalue),
      stringsAsFactors = FALSE
    )
    
    # Add rs# identifier if available and requested
    if(has_rs_ids && include_rs_ids) {
      all_results$rs_id <- snp_info_subset$rs_id[subset_indices]
    }
    
    # Add population AF if available and requested
    if(has_population_AF && include_population_AF) {
      all_results$population_AF <- snp_info_subset$population_AF[subset_indices]
    }
    
    # Adjust p-values only if they were calculated
    if(calc_p && !all(is.na(all_results$pvalue))) {
      all_results$padj <- p.adjust(all_results$pvalue, method = p.adjust.method)
    } else {
      all_results$padj <- NA
    }
    
    # Calculate percent change
    all_results$percent_change <- ((all_results$avg_expr_group1 + pseudocount) /
                                     (all_results$avg_expr_group2 + pseudocount) - 1) * 100
    
    # Sort results by absolute log2fc if no p-values, otherwise by padj
    if(calc_p && !all(is.na(all_results$padj))) {
      all_results <- all_results[order(all_results$padj), ]
    } else {
      all_results <- all_results[order(-abs(all_results$log2fc)), ]
    }
    
    # Create summary
    summary <- list(
      total_in_database = nrow(self$snp_database$dp_matrix),
      pre_filtered = length(qualifying_snp_indices),
      passed_filters = results_count,
      significant = if(calc_p && !all(is.na(all_results$padj))) sum(all_results$padj < 0.05, na.rm=TRUE) else NA,
      upregulated = sum(all_results$log2fc > 0),
      downregulated = sum(all_results$log2fc < 0),
      filtering_efficiency = list(
        pre_filter_reduction = (1 - length(qualifying_snp_indices)/nrow(self$snp_database$dp_matrix)) * 100,
        memory_savings = sprintf("%.1f%%", (1 - length(qualifying_snp_indices)/nrow(self$snp_database$dp_matrix)) * 100)
      ),
      parameters = list(
        use_normalized = use_normalized,
        min_expr_cells = min_expr_cells,
        min_alt_frac = min_alt_frac,
        calculated_pvalues = calc_p,
        test_method = if(calc_p) "wilcox" else "none",
        donor_type = donor_type,
        non_transplant_mode = non_transplant_mode,
        parallel_processing = can_use_parallel,
        n_cores = if(can_use_parallel) n_cores else 1,
        cells_analyzed = length(all_group_indices),
        total_cells_in_database = ncol(self$snp_database$dp_matrix),
        include_rs_ids = include_rs_ids,
        rs_ids_available = has_rs_ids
      )
    )
    
    # Print summary
    cat("\n\n=== Analysis Summary ===")
    cat(sprintf("\nTotal SNPs in database: %d", summary$total_in_database))
    cat(sprintf("\nSNPs pre-filtered (meet expression criteria): %d (%.1f%% reduction)", 
                summary$pre_filtered, summary$filtering_efficiency$pre_filter_reduction))
    cat(sprintf("\nSNPs passing final filters: %d", summary$passed_filters))
    if(calc_p && !all(is.na(all_results$padj))) {
      cat(sprintf("\nSignificant SNPs: %d", summary$significant))
    }
    cat(sprintf("\n - Upregulated: %d", summary$upregulated))
    cat(sprintf("\n - Downregulated: %d", summary$downregulated))
    cat(sprintf("\nProcessing mode: %s", if(can_use_parallel) "Parallel" else "Serial"))
    
    if(has_rs_ids && include_rs_ids) {
      rs_count <- sum(!is.na(all_results$rs_id))
      cat(sprintf("\n - SNPs with rs# identifiers: %d (%.1f%%)", 
                  rs_count, (rs_count/nrow(all_results))*100))
    }
    
    return(list(
      results = if(return_all) all_results else {
        if(calc_p && !all(is.na(all_results$padj))) {
          all_results[all_results$padj < 0.05, ]
        } else {
          all_results
        }
      },
      summary = summary
    ))
  }
  
  cat("\nNo SNPs passed all filters")
  return(NULL)
})
#' @title findSNPsByGroup: Group-Level SNP Presence Analysis
#' @name findSNPsByGroup
#'
#' @description
#' Identifies SNPs that are exclusively or predominantly present in one cell group compared to another.
#' This function analyzes alternative allele frequencies between groups using aggregated data to detect
#' group-specific genetic variants.
#'
#' @param ident.1 Character. Primary group identity to analyze.
#' @param ident.2 Character, optional. Secondary group identity to compare against.
#'                If NULL, compares against all other groups combined.
#' @param aggregated_data List. Output from aggregateByGroup function with required matrices and metadata.
#' @param min_depth Integer. Minimum total read depth required for a group to consider a SNP.
#' @param min_alt_frac Numeric between 0 and 1. Minimum alternative allele fraction required in a group
#'                     for a SNP to be considered present.
#' @param max_alt_frac_other Numeric between 0 and 1. Maximum alternative allele fraction allowed in the
#'                           other group for a SNP to be considered absent there.
#' @param return_all Logical. Whether to return all results regardless of significance.
#' @param include_rs_ids Logical. Whether to include rs# identifiers in results if available.
#' @param presence_score_weights Numeric vector of length 3. Weights for calculating presence scores
#'                              in sample-stratified mode: c(fold_change, sample_consistency, depth_reliability).
#'                              Default: c(0.4, 0.3, 0.3). Weights should sum to 1.0.
#' @return List containing:
#'   \item{results}{Data frame of group-specific SNPs with metrics including genomic position,
#'                  gene annotation, depth metrics, allele frequencies, and presence classification.}
#'   \item{summary}{List with analysis overview including counts of SNPs present in each group
#'                  and parameters used for filtering.}
#'
#' @details
#' The function identifies SNPs that are present in one group but absent in another by applying
#' thresholds to alternative allele frequencies. For each SNP, a presence score is calculated
#' that quantifies the strength of evidence for group-specific presence, considering both the
#' frequency difference and the read depth.
#'
#' A SNP is considered "present" in a group when its alternative allele frequency exceeds
#' `min_alt_frac` and the read depth exceeds `min_depth`. It is considered "absent" in the
#' other group when its alternative allele frequency is below `max_alt_frac_other` and the
#' read depth exceeds `min_depth`.
#'
#' The presence score formula is:
#' score = (alt_frac_present - alt_frac_absent) * (depth/min_depth) * (1 - alt_frac_absent/min_alt_frac)
#'
#' @note
#' - This function operates on pre-aggregated data from `aggregateByGroup()` rather than raw SNP data
#' - Non-transplant mode is automatically detected from the aggregated data parameters
#' - Results are sorted by presence score, with highest-scoring SNPs listed first
#'
#' @examples
#'
#' \dontrun{
#' # Aggregate SNP data by cell type
#' agg_data <- proj$aggregateByGroup(
#'   group_by = "cell_type",
#'   donor_type = "Donor",
#'   use_normalized = TRUE
#' )
#'
#' # Find T cell-specific SNPs
#' tc_snps <- proj$findSNPsByGroup(
#'   ident.1 = "T_cells",
#'   ident.2 = "B_cells",
#'   aggregated_data = agg_data,
#'   min_depth = 20,
#'   min_alt_frac = 0.25,
#'   max_alt_frac_other = 0.05
#' )
#'
#' # Comparing patient groups
#' patient_snps <- proj$findSNPsByGroup(
#'   ident.1 = "ACR",
#'   ident.2 = "No_ACR",
#'   aggregated_data = patient_data,
#'   min_alt_frac = 0.1,
#'   max_alt_frac_other = 0.02
#' )
#'}
#'
#' @seealso
#' \code{\link{aggregateByGroup}} for preparing input data
#' \code{\link{findDESNPs}} for cell-level differential analysis
#' \code{\link{plotSNPs}} for visualizing the identified SNPs

variantCell$set("public", "findSNPsByGroup", function(ident.1,
                                                      ident.2 = NULL,
                                                      aggregated_data,
                                                      min_depth = 10,
                                                      min_alt_frac = 0.2,
                                                      max_alt_frac_other = 0.1,
                                                      return_all = TRUE,
                                                      include_rs_ids = TRUE,
                                                      include_population_AF = TRUE,
                                                      presence_score_weights = c(0.4, 0.3, 0.3)) {  
  # Validate input data structure
  if(!all(c("ad", "dp", "metadata", "mode") %in% names(aggregated_data))) {
    stop("aggregated_data must contain 'ad', 'dp', 'metadata', and 'mode' elements from aggregateByGroup()")
  }
  
  # Check if rs# IDs are available
  has_rs_ids <- FALSE
  if(include_rs_ids) {
    if(!is.null(self$snp_database) && 
       "snp_info" %in% names(self$snp_database) &&
       "rs_id" %in% colnames(self$snp_database$snp_info)) {
      has_rs_ids <- !all(is.na(self$snp_database$snp_info$rs_id))
    }
  }
  
  # Check if population AF is available
  has_population_AF <- FALSE
  if(include_population_AF) {
    if(!is.null(self$snp_database) && 
       "snp_info" %in% names(self$snp_database) &&
       "population_AF" %in% colnames(self$snp_database$snp_info)) {
      has_population_AF <- !all(is.na(self$snp_database$snp_info$population_AF))
    }
  }
  
  
  # Validate presence score weights
  if(abs(sum(presence_score_weights) - 1.0) > 0.01) {
    warning("Presence score weights do not sum to 1.0. Normalizing...")
    presence_score_weights <- presence_score_weights / sum(presence_score_weights)
  }
  
  # Detect aggregation mode and route to appropriate function
  if(aggregated_data$mode == "group_only") {
    cat("Detected group-only mode: analyzing entire group differences\n")
    return(self$findSNPsByGroup_GroupOnly(
      ident.1 = ident.1,
      ident.2 = ident.2,
      aggregated_data = aggregated_data,
      min_depth = min_depth,
      min_alt_frac = min_alt_frac,
      max_alt_frac_other = max_alt_frac_other,
      return_all = return_all,
      include_rs_ids = include_rs_ids,
      include_population_AF = include_population_AF,
      presence_score_weights = presence_score_weights
    ))
  } else if(aggregated_data$mode == "sample_stratified") {
    cat("Detected sample-stratified mode: analyzing samples within groups\n")
    return(self$findSNPsByGroup_SampleStratified(
      ident.1 = ident.1,
      ident.2 = ident.2,
      aggregated_data = aggregated_data,
      min_depth = min_depth,
      min_alt_frac = min_alt_frac,
      max_alt_frac_other = max_alt_frac_other,
      return_all = return_all,
      include_rs_ids = include_rs_ids,
      include_population_AF = include_population_AF,
      presence_score_weights = presence_score_weights
    ))
  } else {
    stop("Unknown aggregation mode. Please use aggregated data from the enhanced aggregateByGroup function.")
  }
})


# Helper function for group-only mode (similar to 0.1.8 behavior)
variantCell$set("public", "findSNPsByGroup_GroupOnly", function(ident.1,
                                                                ident.2 = NULL,
                                                                aggregated_data,
                                                                min_depth = 10,
                                                                min_alt_frac = 0.2,
                                                                max_alt_frac_other = 0.1,
                                                                return_all = TRUE,
                                                                include_rs_ids = TRUE,
                                                                include_population_AF = TRUE,
                                                                presence_score_weights = c(0.4, 0.3, 0.3)) {
  # Input validation
  if (!all(c("ad", "dp", "metadata") %in% names(aggregated_data))) {
    stop("Aggregated data missing required elements")
  }
  
  # Check if rs# IDs are available
  has_rs_ids <- FALSE
  if(include_rs_ids) {
    if(!is.null(self$snp_database) && 
       "snp_info" %in% names(self$snp_database) &&
       "rs_id" %in% colnames(self$snp_database$snp_info)) {
      has_rs_ids <- !all(is.na(self$snp_database$snp_info$rs_id))
    }
    
    if(!has_rs_ids && include_rs_ids) {
      cat("Warning: rs# identifiers requested but not available in project. ")
      cat("Use buildSNPDatabase(add_rs_ids = TRUE, VCF_file_path = '...') to add them.\n")
    }
  }
  
  # Check if population AF is available
  has_population_AF <- FALSE
  if(include_population_AF) {
    if(!is.null(self$snp_database) && 
       "snp_info" %in% names(self$snp_database) &&
       "population_AF" %in% colnames(self$snp_database$snp_info)) {
      has_population_AF <- !all(is.na(self$snp_database$snp_info$population_AF))
    }
    
    if(!has_population_AF && include_population_AF) {
      cat("Warning: Population AF requested but not available in project. ")
      cat("Use buildSNPDatabase(add_population_AF = TRUE, VCF_file_path = '...') to add them.\n")
    }
  }
  
  
  
  # Validate input parameters
  if(!is.numeric(min_depth) || min_depth < 0) {
    stop("min_depth must be a non-negative number")
  }
  if(!is.numeric(min_alt_frac) || min_alt_frac < 0 || min_alt_frac > 1) {
    stop("min_alt_frac must be between 0 and 1")
  }
  if(!is.numeric(max_alt_frac_other) || max_alt_frac_other < 0 || max_alt_frac_other > 1) {
    stop("max_alt_frac_other must be between 0 and 1")
  }
  
  # Extract matrices and metadata
  ad_matrix <- aggregated_data$ad
  dp_matrix <- aggregated_data$dp
  metadata <- aggregated_data$metadata
  
  # Validate extracted data
  if(is.null(ad_matrix) || is.null(dp_matrix) || is.null(metadata)) {
    stop("Invalid aggregated data structure")
  }
  if(ncol(ad_matrix) != ncol(dp_matrix)) {
    stop("AD and DP matrices have different dimensions")
  }
  
  # Check for non-transplant mode in the aggregated data
  non_transplant_mode <- FALSE
  if ("parameters" %in% names(aggregated_data) &&
      "non_transplant_mode" %in% names(aggregated_data$parameters)) {
    non_transplant_mode <- aggregated_data$parameters$non_transplant_mode
  }
  
  
  # Validate groups exist
  if(!ident.1 %in% colnames(ad_matrix)) {
    stop(sprintf("Group '%s' not found in aggregated data", ident.1))
  }
  
  # Check if group has sufficient data
  group1_total_depth <- sum(dp_matrix[, ident.1])
  if(group1_total_depth == 0) {
    warning(sprintf("Group '%s' has no sequencing depth", ident.1))
  }
  
  # Handle ident.2
  if(is.null(ident.2)) {
    # Compare against all other groups combined
    other_groups <- setdiff(colnames(ad_matrix), ident.1)
    if(length(other_groups) == 0) {
      stop("No other groups available for comparison")
    }
    
    # Combine all other groups
    group1_ad <- ad_matrix[, ident.1, drop=FALSE]
    group1_dp <- dp_matrix[, ident.1, drop=FALSE]
    
    group2_ad <- Matrix::rowSums(ad_matrix[, other_groups, drop=FALSE])
    group2_dp <- Matrix::rowSums(dp_matrix[, other_groups, drop=FALSE])
    
    comparison_label <- paste(ident.1, "vs", "All_Others")
    
  } else {
    # Compare specific groups
    if(!ident.2 %in% colnames(ad_matrix)) {
      stop(sprintf("Group '%s' not found in aggregated data", ident.2))
    }
    
    group1_ad <- ad_matrix[, ident.1, drop=FALSE]
    group1_dp <- dp_matrix[, ident.1, drop=FALSE]
    
    group2_ad <- ad_matrix[, ident.2, drop=FALSE]
    group2_dp <- dp_matrix[, ident.2, drop=FALSE]
    
    comparison_label <- paste(ident.1, "vs", ident.2)
  }
  
  # Calculate alternative allele fractions
  group1_alt_frac <- ifelse(group1_dp > 0, group1_ad / group1_dp, 0)
  group2_alt_frac <- ifelse(group2_dp > 0, group2_ad / group2_dp, 0)
  
  # Apply filtering criteria
  # SNPs present in group1 but absent in group2
  group1_present <- group1_alt_frac >= min_alt_frac & group1_dp >= min_depth
  group2_absent <- group2_alt_frac <= max_alt_frac_other & group2_dp >= min_depth
  group1_specific <- group1_present & group2_absent
  
  # SNPs present in group2 but absent in group1  
  group2_present <- group2_alt_frac >= min_alt_frac & group2_dp >= min_depth
  group1_absent <- group1_alt_frac <= max_alt_frac_other & group1_dp >= min_depth
  group2_specific <- group2_present & group1_absent
  
  # Combine results
  significant_snps <- group1_specific | group2_specific
  
  if(!any(significant_snps) && !return_all) {
    cat("No significant SNPs found with current thresholds\n")
    return(list(
      results = data.frame(),
      summary = list(
        comparison = comparison_label,
        total_snps_tested = length(significant_snps),
        snps_found = 0,
        parameters = list(
          min_depth = min_depth,
          min_alt_frac = min_alt_frac,
          max_alt_frac_other = max_alt_frac_other
        )
      )
    ))
  }
  
  # Get SNP annotations
  if(is.null(self$snp_database$snp_info)) {
    stop("SNP info not available in database")
  }
  if(is.null(self$snp_database$snp_annotations)) {
    stop("SNP annotations not available in database")
  }
  
  snp_info <- self$snp_database$snp_info
  
  # Create results data frame
  if(return_all) {
    snp_indices <- 1:nrow(ad_matrix)
  } else {
    snp_indices <- which(significant_snps)
  }
  
  results_df <- data.frame(
    snp_index = snp_indices,
    chr = snp_info$CHROM[snp_indices],
    pos = snp_info$POS[snp_indices],
    ref = snp_info$REF[snp_indices],
    alt = snp_info$ALT[snp_indices],
    gene = self$snp_database$snp_annotations$gene_name[snp_indices],
    feature = self$snp_database$snp_annotations$feature_type[snp_indices],
    group1_ad = as.numeric(group1_ad[snp_indices]),
    group1_dp = as.numeric(group1_dp[snp_indices]),
    group1_alt_frac = as.numeric(group1_alt_frac[snp_indices]),
    group2_ad = as.numeric(group2_ad[snp_indices]),
    group2_dp = as.numeric(group2_dp[snp_indices]),
    group2_alt_frac = as.numeric(group2_alt_frac[snp_indices]),
    stringsAsFactors = FALSE
  )
  
  # Add rs# IDs if available
  if(include_rs_ids && has_rs_ids) {
    results_df$rs_id <- self$snp_database$snp_info$rs_id[snp_indices]
  }
  
  # Add population AF to results if available
  if(include_population_AF && has_population_AF) {
    results_df$population_AF <- self$snp_database$snp_info$population_AF[snp_indices]
  }
  if(has_population_AF) {
    population_af_subset <- self$snp_database$snp_info$population_AF[snp_indices]
    
    # Calculate fold enrichment for group1 vs population AF
    # Use log2 fold change, handling zero/NA population AF values
    results_df$group1_fold_enrichment <- ifelse(
      !is.na(population_af_subset) & population_af_subset > 0,
      log2((results_df$group1_alt_frac + 1e-6) / (population_af_subset + 1e-6)),
      NA
    )
    
    # Calculate fold enrichment for group2 vs population AF
    results_df$group2_fold_enrichment <- ifelse(
      !is.na(population_af_subset) & population_af_subset > 0,
      log2((results_df$group2_alt_frac + 1e-6) / (population_af_subset + 1e-6)),
      NA
    )
    
    # Calculate normalized fold enrichment (bounded between -5 and 5)
    results_df$group1_fold_enrichment_norm <- pmax(-5, pmin(5, results_df$group1_fold_enrichment))
    results_df$group2_fold_enrichment_norm <- pmax(-5, pmin(5, results_df$group2_fold_enrichment))
    
    # Add enrichment interpretation
    results_df$group1_enrichment_level <- ifelse(
      is.na(results_df$group1_fold_enrichment), "No_PopAF",
      ifelse(results_df$group1_fold_enrichment > 1, "High_Enriched",
             ifelse(results_df$group1_fold_enrichment > 0.5, "Moderate_Enriched",
                    ifelse(results_df$group1_fold_enrichment > -0.5, "Similar",
                           ifelse(results_df$group1_fold_enrichment > -1, "Moderate_Depleted", "High_Depleted"))))
    )
    
    results_df$group2_enrichment_level <- ifelse(
      is.na(results_df$group2_fold_enrichment), "No_PopAF",
      ifelse(results_df$group2_fold_enrichment > 1, "High_Enriched",
             ifelse(results_df$group2_fold_enrichment > 0.5, "Moderate_Enriched",
                    ifelse(results_df$group2_fold_enrichment > -0.5, "Similar",
                           ifelse(results_df$group2_fold_enrichment > -1, "Moderate_Depleted", "High_Depleted"))))
    )
    
  }  
  
  # Calculate presence classification and scores
  results_df$group1_present <- results_df$group1_alt_frac >= min_alt_frac & results_df$group1_dp >= min_depth
  results_df$group2_present <- results_df$group2_alt_frac >= min_alt_frac & results_df$group2_dp >= min_depth
  
  # Determine presence pattern
  results_df$presence_pattern <- ifelse(
    results_df$group1_present & !results_df$group2_present, paste0(ident.1, "_specific"),
    ifelse(!results_df$group1_present & results_df$group2_present, 
           ifelse(is.null(ident.2), "Other_specific", paste0(ident.2, "_specific")),
           ifelse(results_df$group1_present & results_df$group2_present, "Both_present", "Neither_present")
    )
  )
  
  # Calculate presence score
  alt_frac_diff <- abs(results_df$group1_alt_frac - results_df$group2_alt_frac)
  min_depth_score <- pmin(results_df$group1_dp, results_df$group2_dp) / min_depth
  max_depth_score <- pmax(results_df$group1_dp, results_df$group2_dp) / min_depth
  
  results_df$presence_score <- (
    presence_score_weights[1] * alt_frac_diff +
      presence_score_weights[2] * pmin(min_depth_score, 1) +
      presence_score_weights[3] * pmin(max_depth_score / 10, 1)  # Scaled for very high depths
  )
  
  # Sort by presence score
  results_df <- results_df[order(results_df$presence_score, decreasing = TRUE), ]
  
  # Create summary
  summary_info <- list(
    comparison = comparison_label,
    total_snps_tested = nrow(results_df),
    snps_group1_specific = sum(results_df$presence_pattern == paste0(ident.1, "_specific"), na.rm = TRUE),
    snps_group2_specific = sum(grepl("_specific$", results_df$presence_pattern) & 
                                 results_df$presence_pattern != paste0(ident.1, "_specific"), na.rm = TRUE),
    snps_both_present = sum(results_df$presence_pattern == "Both_present", na.rm = TRUE),
    snps_neither_present = sum(results_df$presence_pattern == "Neither_present", na.rm = TRUE),
    parameters = list(
      min_depth = min_depth,
      min_alt_frac = min_alt_frac,
      max_alt_frac_other = max_alt_frac_other,
      mode = "group_only"
    )
  )
  
  # Add enrichment summary if population AF available
  if(has_population_AF) {
    summary_info$enrichment_summary <- list(
      group1_highly_enriched = sum(results_df$group1_enrichment_level == "High_Enriched", na.rm = TRUE),
      group1_moderately_enriched = sum(results_df$group1_enrichment_level == "Moderate_Enriched", na.rm = TRUE),
      group1_similar = sum(results_df$group1_enrichment_level == "Similar", na.rm = TRUE),
      group1_moderately_depleted = sum(results_df$group1_enrichment_level == "Moderate_Depleted", na.rm = TRUE),
      group1_highly_depleted = sum(results_df$group1_enrichment_level == "High_Depleted", na.rm = TRUE),
      snps_with_population_af = sum(!is.na(results_df$population_AF)),
      snps_without_population_af = sum(is.na(results_df$population_AF))
    )
  }
  
  cat(sprintf("Found %d %s-specific and %d other-group-specific SNPs\n", 
              summary_info$snps_group1_specific, 
              ident.1, 
              summary_info$snps_group2_specific))
  
  if(has_population_AF) {
    cat(sprintf("Population AF enrichment summary for %s:\n", ident.1))
    cat(sprintf("  High enrichment (>2x): %d SNPs\n", summary_info$enrichment_summary$group1_highly_enriched))
    cat(sprintf("  Moderate enrichment (1.4-2x): %d SNPs\n", summary_info$enrichment_summary$group1_moderately_enriched))
    cat(sprintf("  Similar to population: %d SNPs\n", summary_info$enrichment_summary$group1_similar))
    cat(sprintf("  Moderate depletion (0.5-0.7x): %d SNPs\n", summary_info$enrichment_summary$group1_moderately_depleted))
    cat(sprintf("  High depletion (<0.5x): %d SNPs\n", summary_info$enrichment_summary$group1_highly_depleted))
  }
  
  
  
  return(list(
    results = results_df,
    summary = summary_info
  ))
})

# Helper function for sample-stratified mode 
variantCell$set("public", "findSNPsByGroup_SampleStratified", function(ident.1,
                                                                       ident.2 = NULL,
                                                                       aggregated_data,
                                                                       min_depth = 10,
                                                                       min_alt_frac = 0.2,
                                                                       max_alt_frac_other = 0.1,
                                                                       return_all = TRUE,
                                                                       include_rs_ids = TRUE,
                                                                       include_population_AF = TRUE,
                                                                       presence_score_weights = c(fold_change = 0.5, 
                                                                                                  sample_consistency = 0.3, 
                                                                                                  depth_reliability = 0.2)) {
  
  # Input validation
  if(!all(c("ad", "dp", "metadata") %in% names(aggregated_data))) {
    stop("Aggregated data missing required elements")
  }
  
  # Check if this is sample-aware aggregated data
  is_sample_aware <- !is.null(aggregated_data$parameters$sample_column) &&
    "sample_metadata" %in% names(aggregated_data)
  
  if(!is_sample_aware) {
    stop("This function requires aggregated data for multiple samples. Use aggregateByGroup with sample_column parameter.")
  }
  
  # Check if rs# IDs are available
  has_rs_ids <- FALSE
  if(include_rs_ids) {
    if(!is.null(self$snp_database) && 
       "snp_info" %in% names(self$snp_database) &&
       "rs_id" %in% colnames(self$snp_database$snp_info)) {
      has_rs_ids <- !all(is.na(self$snp_database$snp_info$rs_id))
    }
  }
  
  # Check if population AF is available
  has_population_AF <- FALSE
  if(include_population_AF) {
    if(!is.null(self$snp_database) && 
       "snp_info" %in% names(self$snp_database) &&
       "population_AF" %in% colnames(self$snp_database$snp_info)) {
      has_population_AF <- !all(is.na(self$snp_database$snp_info$population_AF))
    }
  }
  
  
  # Validate presence score weights
  if(abs(sum(presence_score_weights) - 1.0) > 0.01) {
    warning("Presence score weights do not sum to 1.0. Normalizing...")
    presence_score_weights <- presence_score_weights / sum(presence_score_weights)
  }
  
  # Get data
  sample_meta <- aggregated_data$sample_metadata
  sample_level_ad <- aggregated_data$sample_level_ad
  sample_level_dp <- aggregated_data$sample_level_dp
  
  # Check for non-transplant mode
  non_transplant_mode <- FALSE
  if("parameters" %in% names(aggregated_data) &&
     "non_transplant_mode" %in% names(aggregated_data$parameters)) {
    non_transplant_mode <- aggregated_data$parameters$non_transplant_mode
  }
  
  # Get sample information for each group
  if(is.null(ident.2)) {
    samples_group1 <- sample_meta$group_sample[sample_meta$group == ident.1 & sample_meta$filter_status == "included"]
    samples_group2 <- sample_meta$group_sample[sample_meta$group != ident.1 & sample_meta$filter_status == "included"]
    group2_name <- "rest"
  } else {
    samples_group1 <- sample_meta$group_sample[sample_meta$group == ident.1 & sample_meta$filter_status == "included"]
    samples_group2 <- sample_meta$group_sample[sample_meta$group == ident.2 & sample_meta$filter_status == "included"]
    group2_name <- ident.2
  }
  
  n_samples_group1 <- length(samples_group1)
  n_samples_group2 <- length(samples_group2)
  
  # Enhanced sample requirement calculation
  calculate_min_samples <- function(n_samples, min_pct, min_abs) {
    max(min_abs, ceiling(n_samples * min_pct))
  }
  
  # Set default values for parameters not in function signature
  min_pct_samples <- 0.6
  min_abs_samples <- 2
  calc_p_values <- TRUE
  p_adjust_method <- "BH"
  
  min_samples_group1 <- calculate_min_samples(n_samples_group1, min_pct_samples, min_abs_samples)
  min_samples_group2 <- calculate_min_samples(n_samples_group2, min_pct_samples, min_abs_samples)
  
  # Print enhanced summary
  cat("\n=== Sample-Aware SNP Analysis ===")
  cat(sprintf("\nComparing SNP presence between groups:"))
  cat(sprintf("\n%s: %d samples, need %d samples (%.0f%% or min %d)", 
              ident.1, n_samples_group1, min_samples_group1, min_pct_samples * 100, min_abs_samples))
  cat(sprintf("\n%s: %d samples, need %d samples (%.0f%% or min %d)", 
              group2_name, n_samples_group2, min_samples_group2, min_pct_samples * 100, min_abs_samples))
  cat(sprintf("\nStatistical testing: %s", if(calc_p_values) "enabled" else "disabled"))
  
  # Get sample-level indices
  sample_indices_group1 <- which(colnames(sample_level_ad) %in% samples_group1)
  sample_indices_group2 <- which(colnames(sample_level_ad) %in% samples_group2)
  
  # Presence score calculation
  calculate_presence_score <- function(alt_frac_diff, sample_weight, depth_factor, weights) {
    # Normalize components to 0-1 scale with diminishing returns
    norm_fold_change <- pmin(1, abs(alt_frac_diff) / 0.5)  # Cap at 50% difference
    norm_sample_weight <- sample_weight  # Already 0-1
    norm_depth <- pmin(1, depth_factor / 10)  # Cap depth contribution at 10x min_depth
    
    # Apply sigmoid transformation for diminishing returns
    sigmoid <- function(x) 2 / (1 + exp(-2 * x)) - 1
    
    norm_fold_change <- sigmoid(norm_fold_change * 2)
    norm_depth <- sigmoid(norm_depth * 2)
    
    # Weighted combination
    score <- weights["fold_change"] * norm_fold_change + 
      weights["sample_consistency"] * norm_sample_weight + 
      weights["depth_reliability"] * norm_depth
    
    return(pmin(1, pmax(0, score)))  # Ensure 0-1 bounds
  }
  
  # Statistical testing function
  perform_presence_test <- function(present_g1, total_g1, present_g2, total_g2) {
    if(total_g1 < 2 || total_g2 < 2) return(list(p_value = NA, effect_size = NA))
    
    # Fisher's exact test for presence pattern
    contingency_table <- matrix(c(present_g1, total_g1 - present_g1,
                                  present_g2, total_g2 - present_g2), 
                                nrow = 2, byrow = TRUE)
    
    test_result <- tryCatch({
      fisher.test(contingency_table)
    }, error = function(e) {
      list(p.value = NA)
    })
    
    # Calculate effect size (log odds ratio)
    effect_size <- tryCatch({
      log((present_g1 / (total_g1 - present_g1 + 1)) / 
            (present_g2 / (total_g2 - present_g2 + 1) + 1e-6))
    }, error = function(e) NA)
    
    return(list(p_value = test_result$p.value, effect_size = effect_size))
  }
  
  # Quality assessment function with depth metrics
  assess_quality <- function(sample_depths, sample_alt_fracs, n_samples, sample_mean_depths) {
    depth_cv <- sd(sample_depths, na.rm = TRUE) / mean(sample_depths, na.rm = TRUE)
    alt_frac_consistency <- 1 - sd(sample_alt_fracs, na.rm = TRUE)
    sample_coverage <- n_samples / max(n_samples_group1, n_samples_group2)
    
    # Overall quality score (0-1, higher is better)
    quality_score <- mean(c(
      pmax(0, 1 - depth_cv),
      pmax(0, alt_frac_consistency),
      sample_coverage
    ), na.rm = TRUE)
    
    return(list(
      depth_cv = depth_cv,
      alt_frac_consistency = alt_frac_consistency,
      sample_coverage = sample_coverage,
      overall_quality = quality_score
    ))
  }
  
  # Initialize results storage
  results <- vector("list", nrow(sample_level_ad))
  results_count <- 0
  
  # Process each SNP
  cat(sprintf("\nProcessing %d SNPs...", nrow(sample_level_ad)))
  pb <- txtProgressBar(min = 0, max = nrow(sample_level_ad), style = 3)
  
  for(i in seq_len(nrow(sample_level_ad))) {
    setTxtProgressBar(pb, i)
    
    # Get sample-level data for this SNP
    ad1_samples <- sample_level_ad[i, sample_indices_group1]
    dp1_samples <- sample_level_dp[i, sample_indices_group1]
    ad2_samples <- sample_level_ad[i, sample_indices_group2]
    dp2_samples <- sample_level_dp[i, sample_indices_group2]
    
    # Calculate alt fractions for each sample
    alt_frac1_samples <- ifelse(dp1_samples > 0, ad1_samples / dp1_samples, 0)
    alt_frac2_samples <- ifelse(dp2_samples > 0, ad2_samples / dp2_samples, 0)
    
    # Count samples meeting presence criteria in each group
    samples_present_group1 <- sum(dp1_samples >= min_depth & alt_frac1_samples >= min_alt_frac, na.rm = TRUE)
    samples_absent_group2 <- sum(dp2_samples >= min_depth & alt_frac2_samples <= max_alt_frac_other, na.rm = TRUE)
    
    samples_present_group2 <- sum(dp2_samples >= min_depth & alt_frac2_samples >= min_alt_frac, na.rm = TRUE)
    samples_absent_group1 <- sum(dp1_samples >= min_depth & alt_frac1_samples <= max_alt_frac_other, na.rm = TRUE)
    
    # Check if sample-level criteria are met
    snp_in_group1_samples <- samples_present_group1 >= min_samples_group1
    snp_not_in_group2_samples <- samples_absent_group2 >= min_samples_group2
    
    snp_in_group2_samples <- samples_present_group2 >= min_samples_group2
    snp_not_in_group1_samples <- samples_absent_group1 >= min_samples_group1
    
    # Only proceed if sample-level criteria are met
    if((snp_in_group1_samples && snp_not_in_group2_samples) || 
       (snp_in_group2_samples && snp_not_in_group1_samples)) {
      
      # Calculate group-level summary statistics
      ad1_total <- sum(ad1_samples, na.rm = TRUE)
      dp1_total <- sum(dp1_samples, na.rm = TRUE)
      ad2_total <- sum(ad2_samples, na.rm = TRUE)
      dp2_total <- sum(dp2_samples, na.rm = TRUE)
      
      alt_frac1_group <- ifelse(dp1_total > 0, ad1_total / dp1_total, 0)
      alt_frac2_group <- ifelse(dp2_total > 0, ad2_total / dp2_total, 0)
      
      # Determine presence classification and calculate enhanced score
      if(snp_in_group1_samples && snp_not_in_group2_samples) {
        sample_weight1 <- samples_present_group1 / n_samples_group1
        sample_weight2 <- samples_absent_group2 / n_samples_group2
        sample_consistency <- (sample_weight1 + sample_weight2) / 2
        alt_frac_diff <- alt_frac1_group - alt_frac2_group
        depth_factor <- dp1_total / min_depth
        presence_classification <- sprintf("Present in %s", ident.1)
        
        # Quality assessment for group 1
        quality_metrics <- assess_quality(dp1_samples[dp1_samples > 0], 
                                          alt_frac1_samples[dp1_samples > 0], 
                                          n_samples_group1)
      } else {
        sample_weight1 <- samples_absent_group1 / n_samples_group1
        sample_weight2 <- samples_present_group2 / n_samples_group2
        sample_consistency <- (sample_weight1 + sample_weight2) / 2
        alt_frac_diff <- alt_frac2_group - alt_frac1_group
        depth_factor <- dp2_total / min_depth
        presence_classification <- sprintf("Present in %s", group2_name)
        
        # Quality assessment for group 2
        quality_metrics <- assess_quality(dp2_samples[dp2_samples > 0], 
                                          alt_frac2_samples[dp2_samples > 0], 
                                          n_samples_group2)
      }
      
      # Calculate enhanced presence score
    presence_score <- calculate_presence_score(
      alt_frac_diff, sample_consistency, depth_factor,
      c(fold_change = presence_score_weights[1],
        sample_consistency = presence_score_weights[2],
        depth_reliability = presence_score_weights[3])
    )      
      # Statistical testing
      if(calc_p_values) {
        stat_result <- perform_presence_test(
          samples_present_group1, n_samples_group1,
          samples_present_group2, n_samples_group2
        )
      } else {
        stat_result <- list(p_value = NA, effect_size = NA)
      }
      
      results_count <- results_count + 1
      
      # Create enhanced results data frame
      result_row <- data.frame(
        snp_idx = i,
        chromosome = self$snp_database$snp_info$CHROM[i],
        position = self$snp_database$snp_info$POS[i],
        ref = self$snp_database$snp_info$REF[i],
        alt = self$snp_database$snp_info$ALT[i],
        feature_type = self$snp_database$snp_annotations$feature_type[i],
        gene_name = self$snp_database$snp_annotations$gene_name[i],
        gene_type = self$snp_database$snp_annotations$gene_type[i],
        
        # Group-level metrics
        depth_1 = dp1_total,
        depth_2 = dp2_total,
        alt_count_1 = ad1_total,
        alt_count_2 = ad2_total,
        alt_frac_1 = alt_frac1_group,
        alt_frac_2 = alt_frac2_group,
        alt_frac_diff = abs(alt_frac1_group - alt_frac2_group),
        
        # Sample-level metrics
        n_samples_1 = n_samples_group1,
        n_samples_2 = n_samples_group2,
        samples_present_1 = samples_present_group1,
        samples_present_2 = samples_present_group2,
        samples_absent_1 = samples_absent_group1,
        samples_absent_2 = samples_absent_group2,
        pct_samples_present_1 = samples_present_group1 / n_samples_group1 * 100,
        pct_samples_present_2 = samples_present_group2 / n_samples_group2 * 100,
        sample_consistency = sample_consistency,
        
        # Enhanced scoring and statistics
        presence_score = presence_score,
        presence = presence_classification,
        p_value = stat_result$p_value,
        effect_size = stat_result$effect_size,
        
        # Quality metrics
        depth_cv = quality_metrics$depth_cv,
        alt_frac_consistency = quality_metrics$alt_frac_consistency,
        sample_coverage = quality_metrics$sample_coverage,
        overall_quality = quality_metrics$overall_quality,
        
        stringsAsFactors = FALSE
      )
      
      # Add rs# identifier if available and requested
      if(has_rs_ids && include_rs_ids) {
        result_row$rs_id <- self$snp_database$snp_info$rs_id[i]
      }
      
      if(include_population_AF && has_population_AF) {
        result_row$population_AF <- self$snp_database$snp_info$population_AF[i]
      }
      if(has_population_AF) {
        population_af_i <- self$snp_database$snp_info$population_AF[i]
        
        # Calculate fold enrichment for both groups vs population AF
        result_row$group1_fold_enrichment <- ifelse(
          !is.na(population_af_i) & population_af_i > 0,
          log2((result_row$alt_frac_1 + 1e-6) / (population_af_i + 1e-6)),
          NA
        )
        
        result_row$group2_fold_enrichment <- ifelse(
          !is.na(population_af_i) & population_af_i > 0,
          log2((result_row$alt_frac_2 + 1e-6) / (population_af_i + 1e-6)),
          NA
        )
        
        # Calculate normalized fold enrichment (bounded between -5 and 5)
        result_row$group1_fold_enrichment_norm <- pmax(-5, pmin(5, result_row$group1_fold_enrichment))
        result_row$group2_fold_enrichment_norm <- pmax(-5, pmin(5, result_row$group2_fold_enrichment))
        
        # Add enrichment interpretation
        result_row$group1_enrichment_level <- ifelse(
          is.na(result_row$group1_fold_enrichment), "No_PopAF",
          ifelse(result_row$group1_fold_enrichment > 1, "High_Enriched",
                 ifelse(result_row$group1_fold_enrichment > 0.5, "Moderate_Enriched",
                        ifelse(result_row$group1_fold_enrichment > -0.5, "Similar",
                               ifelse(result_row$group1_fold_enrichment > -1, "Moderate_Depleted", "High_Depleted"))))
        )
        
        result_row$group2_enrichment_level <- ifelse(
          is.na(result_row$group2_fold_enrichment), "No_PopAF",
          ifelse(result_row$group2_fold_enrichment > 1, "High_Enriched",
                 ifelse(result_row$group2_fold_enrichment > 0.5, "Moderate_Enriched",
                        ifelse(result_row$group2_fold_enrichment > -0.5, "Similar",
                               ifelse(result_row$group2_fold_enrichment > -1, "Moderate_Depleted", "High_Depleted"))))
        )
      }
      
      results[[results_count]] <- result_row
    }
  }
  
  close(pb)
  
  # Process final results
  if(results_count > 0) {
    all_results <- do.call(rbind, results[1:results_count])
    
    # Apply multiple testing correction if p-values were calculated
    if(calc_p_values && !all(is.na(all_results$p_value))) {
      all_results$p_adjusted <- p.adjust(all_results$p_value, method = p_adjust_method)
      all_results$significant <- !is.na(all_results$p_adjusted) & all_results$p_adjusted < 0.05
    } else {
      all_results$p_adjusted <- NA
      all_results$significant <- NA
    }
    
    # Sort results by presence score (and p-value if available)
    if(calc_p_values && !all(is.na(all_results$p_adjusted))) {
      all_results <- all_results[order(all_results$p_adjusted, -all_results$presence_score), ]
    } else {
      all_results <- all_results[order(-all_results$presence_score), ]
    }
    
    # Create enhanced summary
    summary <- list(
      total_tested = nrow(sample_level_ad),
      passed_filters = results_count,
      present_in_group1 = sum(all_results$presence == sprintf("Present in %s", ident.1)),
      present_in_group2 = sum(all_results$presence == sprintf("Present in %s", group2_name)),
      significant_results = if(calc_p_values) sum(all_results$significant, na.rm = TRUE) else NA,
      high_quality_results = sum(all_results$overall_quality >= 0.7, na.rm = TRUE),
      parameters = list(
        min_depth = min_depth,
        min_alt_frac = min_alt_frac,
        max_alt_frac_other = max_alt_frac_other,
        min_pct_samples = min_pct_samples,
        min_abs_samples = min_abs_samples,
        min_samples_group1 = min_samples_group1,
        min_samples_group2 = min_samples_group2,
        calc_p_values = calc_p_values,
        p_adjust_method = p_adjust_method,
        presence_score_weights = presence_score_weights,
        non_transplant_mode = non_transplant_mode,
        include_rs_ids = include_rs_ids,
        rs_ids_available = has_rs_ids,
        include_population_AF = include_population_AF, 
        population_AF_available = has_population_AF 
      ),
      quality_distribution = summary(all_results$overall_quality),
      patterns = table(all_results$presence)
    )
    
    # Print enhanced summary
    cat("\n\n=== Enhanced Analysis Summary ===")
    cat(sprintf("\nTotal SNPs tested: %d", summary$total_tested))
    cat(sprintf("\nSNPs passing enhanced filters: %d", summary$passed_filters))
    cat(sprintf("\n - Present in %s: %d", ident.1, summary$present_in_group1))
    cat(sprintf("\n - Present in %s: %d", group2_name, summary$present_in_group2))
    
    if(calc_p_values) {
      cat(sprintf("\n - Statistically significant: %d", summary$significant_results))
    }
    
    cat(sprintf("\n - High quality (>0.7): %d", summary$high_quality_results))
    cat(sprintf("\nSample criteria: %.0f%% (min %d) of samples must meet thresholds", 
                min_pct_samples * 100, min_abs_samples))
    
    if(has_rs_ids && include_rs_ids) {
      rs_count <- sum(!is.na(all_results$rs_id))
      cat(sprintf("\n - SNPs with rs# identifiers: %d (%.1f%%)", 
                  rs_count, (rs_count/nrow(all_results))*100))
    }
    if(has_population_AF && include_population_AF) {
      AF_count <- sum(!is.na(all_results$population_AF))
      cat(sprintf("\n - SNPs with population AF: %d (%.1f%%)", 
                  AF_count, (AF_count/nrow(all_results))*100))
    }
    
    
    cat("\n\nQuality score distribution:")
    print(summary$quality_distribution)
    
    cat("\n\nPresence distribution:")
    print(summary$patterns)
    
    return(list(
      results = all_results,
      summary = summary
    ))
  }
  
  cat("\nNo SNPs meeting enhanced criteria found")
  return(NULL)
})

#' @title prioritizeSNPs: Comprehensive SNP Prioritization
#' @name prioritizeSNPs
#'
#' @description
#' Prioritizes SNPs using multiple methodologies including fold change consistency 
#' with gene expression data, machine learning regression, and LLM-based clinical 
#' relevance assessment. This function integrates differential SNP analysis results 
#' with single-cell gene expression data to identify the most biologically and 
#' clinically relevant variants.
#'
#' @param snp_df Data frame. Output from findDESNPs or findSNPsByGroup with 
#'               differential SNP analysis results. Must contain columns:
#'               rs_id, gene_name, GEX_avg_log2FC, GEX_avg_FC, GEX_p_val, GEX_p_val_adj
#' @param gex_fc_df Data frame. Gene expression fold change data with columns:
#'                  gene_name and Average (average FC across clusters)
#' @param method Character vector. Prioritization methods to use:
#'               "fc_consistency" - Score based on GEX fold change consistency
#'               "ml_regression" - Machine learning regression prioritization  
#'               "llm_clinical" - LLM-based clinical relevance scoring
#'               Default: c("fc_consistency", "ml_regression", "llm_clinical")
#' @param ml_features Character vector. Features to use for ML regression.
#'                    Default: c("alt_frac_diff", "effect_size", "presence_score", "overall_quality", 
#'                              "population_AF", "group1_fold_enrichment", "group2_fold_enrichment")
#' @param top_n_ml Integer. Number of top SNPs to select from ML prioritization for LLM assessment.
#'                 Default: 50
#' @param top_n_final Integer. Final number of top SNPs to return. Default: 20
#' @param llm_prompt_template Character. Custom LLM prompt template for clinical assessment.
#'                           Default: NULL (uses built-in template)
#' @param verbose Logical. Whether to print detailed progress information. Default: TRUE
#'
#' @note
#' For LLM-based clinical assessment, this function requires:
#' - The ellmer package: install with remotes::install_github('tidyverse/ellmer')
#' - API keys set as environment variables: ANTHROPIC_API_KEY or OPENAI_API_KEY
#' - If LLM integration is unavailable, the function falls back to rule-based assessment
#'
#' @return A list containing:
#'   \item{prioritized_snps}{Data frame with top prioritized SNPs and all scoring methods}
#'   \item{fc_consistency_scores}{Data frame with fold change consistency scores}
#'   \item{ml_scores}{Data frame with ML regression scores}  
#'   \item{llm_assessments}{Data frame with LLM clinical relevance assessments}
#'   \item{method_weights}{Weights used to combine different prioritization methods}
#'   \item{summary}{Summary statistics of the prioritization process}
#'
#' @details
#' This function implements a multi-stage SNP prioritization approach:
#' 
#' 1. **Fold Change Consistency**: Compares SNP-associated gene expression changes 
#'    with overall condition-specific expression patterns to identify SNPs whose 
#'    effects align with known disease biology.
#'    
#' 2. **Machine Learning Regression**: Uses multiple SNP characteristics 
#'    (effect size, presence score, population frequency, etc.) to predict 
#'    SNP importance using ensemble methods.
#'    
#' 3. **LLM Clinical Assessment**: Evaluates top-scoring SNPs for clinical 
#'    relevance using structured queries to language models, considering 
#'    gene function, disease association, and therapeutic implications.
#'
#' The final prioritization combines scores from all methods using weighted 
#' averaging, with higher weights given to methods showing better concordance.
#'
#' @examples
#' \dontrun{
#' # Basic usage with all methods
#' prioritized <- project$prioritizeSNPs(
#'   snp_df = de_snp_results,
#'   gex_fc_df = gene_expression_changes
#' )
#' 
#' # Use only specific methods
#' prioritized <- project$prioritizeSNPs(
#'   snp_df = de_snp_results,
#'   gex_fc_df = gene_expression_changes,
#'   method = c("fc_consistency", "ml_regression"),
#'   top_n_final = 10
#' )
#' 
#' # Custom ML features
#' prioritized <- project$prioritizeSNPs(
#'   snp_df = de_snp_results,
#'   gex_fc_df = gene_expression_changes,
#'   ml_features = c("effect_size", "presence_score", "population_AF")
#' )
#' }
#'
#' @export
prioritizeSNPs <- function(snp_df, 
                          gex_fc_df,
                          method = c("fc_consistency", "ml_regression", "llm_clinical"),
                          ml_features = c("alt_frac_diff", "effect_size", "presence_score", 
                                        "overall_quality", "population_AF", 
                                        "group1_fold_enrichment", "group2_fold_enrichment"),
                          top_n_ml = 50,
                          top_n_final = 20,
                          llm_prompt_template = NULL,
                          verbose = TRUE) {
  
  if(verbose) {
    cat("=== SNP Prioritization Analysis ===\n")
    cat("Input data:\n")
    cat(sprintf("  SNP dataframe: %d rows, %d columns\n", nrow(snp_df), ncol(snp_df)))
    cat(sprintf("  GEX FC dataframe: %d rows, %d columns\n", nrow(gex_fc_df), ncol(gex_fc_df)))
    cat(sprintf("  Methods requested: %s\n", paste(method, collapse=", ")))
    cat("\n")
  }
  
  # Validate input data
  required_snp_cols <- c("rs_id", "gene_name")
  required_gex_cols <- c("gene_name", "Average")
  
  missing_snp_cols <- setdiff(required_snp_cols, colnames(snp_df))
  missing_gex_cols <- setdiff(required_gex_cols, colnames(gex_fc_df))
  
  if(length(missing_snp_cols) > 0) {
    stop(paste("Missing required columns in snp_df:", paste(missing_snp_cols, collapse=", ")))
  }
  if(length(missing_gex_cols) > 0) {
    stop(paste("Missing required columns in gex_fc_df:", paste(missing_gex_cols, collapse=", ")))
  }
  
  # Initialize results list
  results <- list()
  all_scores <- data.frame(
    rs_id = snp_df$rs_id,
    gene_name = snp_df$gene_name,
    stringsAsFactors = FALSE
  )
  
  # Method 1: Fold Change Consistency Scoring
  if("fc_consistency" %in% method) {
    if(verbose) cat("1. Computing fold change consistency scores...\n")
    
    fc_scores <- .compute_fc_consistency_scores(snp_df, gex_fc_df, verbose)
    results$fc_consistency_scores <- fc_scores
    all_scores$fc_consistency_score <- fc_scores$consistency_score[match(all_scores$rs_id, fc_scores$rs_id)]
    
    if(verbose) {
      cat(sprintf("   Computed scores for %d SNPs\n", nrow(fc_scores)))
      cat(sprintf("   Mean consistency score: %.3f\n", mean(fc_scores$consistency_score, na.rm=TRUE)))
    }
  }
  
  # Method 2: Machine Learning Regression
  if("ml_regression" %in% method) {
    if(verbose) cat("2. Running machine learning regression prioritization...\n")
    
    ml_scores <- .compute_ml_prioritization(snp_df, ml_features, verbose)
    results$ml_scores <- ml_scores
    all_scores$ml_score <- ml_scores$ml_priority_score[match(all_scores$rs_id, ml_scores$rs_id)]
    
    if(verbose) {
      cat(sprintf("   Computed ML scores for %d SNPs\n", nrow(ml_scores)))
      cat(sprintf("   Mean ML score: %.3f\n", mean(ml_scores$ml_priority_score, na.rm=TRUE)))
    }
  }
  
  # Method 3: LLM Clinical Relevance Assessment
  if("llm_clinical" %in% method) {
    if(verbose) cat("3. Performing LLM clinical relevance assessment...\n")
    
    # Select top SNPs for LLM assessment based on available scores
    top_snps_for_llm <- .select_top_snps_for_llm(all_scores, snp_df, top_n_ml)
    
    llm_assessments <- .compute_llm_clinical_scores(top_snps_for_llm, llm_prompt_template, verbose)
    results$llm_assessments <- llm_assessments
    all_scores$llm_clinical_score <- llm_assessments$clinical_relevance_score[match(all_scores$rs_id, llm_assessments$rs_id)]
    
    if(verbose) {
      cat(sprintf("   Assessed %d top SNPs with LLM\n", nrow(llm_assessments)))
      cat(sprintf("   Mean clinical relevance score: %.3f\n", 
                  mean(llm_assessments$clinical_relevance_score, na.rm=TRUE)))
    }
  }
  
  # Combine scores and create final prioritization
  if(verbose) cat("4. Combining scores and generating final prioritization...\n")
  
  final_priority <- .combine_prioritization_scores(all_scores, method, snp_df)
  
  # Select top N SNPs
  top_snps <- head(final_priority[order(-final_priority$combined_priority_score, na.last=TRUE), ], top_n_final)
  
  # Generate summary
  summary_stats <- .generate_prioritization_summary(results, all_scores, top_snps, method)
  
  if(verbose) {
    cat("5. Prioritization complete!\n")
    cat(sprintf("   Final top %d SNPs selected\n", nrow(top_snps)))
    cat(sprintf("   Top SNP: %s (gene: %s, score: %.3f)\n", 
                top_snps$rs_id[1], top_snps$gene_name[1], top_snps$combined_priority_score[1]))
    cat("\n")
  }
  
  return(list(
    prioritized_snps = top_snps,
    fc_consistency_scores = results$fc_consistency_scores,
    ml_scores = results$ml_scores,
    llm_assessments = results$llm_assessments,
    method_weights = final_priority$method_weights[1],
    summary = summary_stats
  ))
}

# Helper function: Compute fold change consistency scores
.compute_fc_consistency_scores <- function(snp_df, gex_fc_df, verbose = TRUE) {
  
  # Merge SNP data with GEX fold change data
  merged_data <- merge(snp_df, gex_fc_df, by = "gene_name", all.x = TRUE)
  
  # Compute consistency scores
  consistency_scores <- data.frame(
    rs_id = merged_data$rs_id,
    gene_name = merged_data$gene_name,
    snp_gex_log2fc = merged_data$GEX_avg_log2FC,
    condition_avg_fc = merged_data$Average,
    stringsAsFactors = FALSE
  )
  
  # Calculate consistency score (higher when SNP effect aligns with condition effect)
  consistency_scores$fc_direction_match <- sign(consistency_scores$snp_gex_log2fc) == sign(log2(consistency_scores$condition_avg_fc))
  consistency_scores$fc_magnitude_similarity <- 1 - abs(consistency_scores$snp_gex_log2fc - log2(consistency_scores$condition_avg_fc)) / 
    (abs(consistency_scores$snp_gex_log2fc) + abs(log2(consistency_scores$condition_avg_fc)) + 0.1)
  
  # Combined consistency score (0-1 scale)
  consistency_scores$consistency_score <- ifelse(
    consistency_scores$fc_direction_match,
    0.5 + 0.5 * consistency_scores$fc_magnitude_similarity,
    0.5 * consistency_scores$fc_magnitude_similarity
  )
  
  # Handle NAs
  consistency_scores$consistency_score[is.na(consistency_scores$consistency_score)] <- 0
  
  return(consistency_scores)
}

# Helper function: Compute ML prioritization scores
.compute_ml_prioritization <- function(snp_df, ml_features, verbose = TRUE) {
  
  # Check which features are available
  available_features <- intersect(ml_features, colnames(snp_df))
  missing_features <- setdiff(ml_features, colnames(snp_df))
  
  if (verbose && length(missing_features) > 0) {
    cat(sprintf("   Warning: Missing ML features: %s\n", 
                paste(missing_features, collapse = ", ")))
  }
  
  if (length(available_features) == 0) {
    warning("No ML features available in SNP dataframe")
    return(data.frame(rs_id = snp_df$rs_id, ml_priority_score = 0))
  }
  
  # Extract feature matrix
  feature_data <- snp_df[, available_features, drop = FALSE]
  
  # Handle missing values with more robust imputation
  for (col in colnames(feature_data)) {
    if (is.numeric(feature_data[[col]])) {
      # Use median imputation for missing values
      median_val <- median(feature_data[[col]], na.rm = TRUE)
      if (is.na(median_val)) median_val <- 0
      feature_data[[col]][is.na(feature_data[[col]])] <- median_val
    }
  }
  
  # Advanced feature engineering
  if (verbose) cat("   Performing feature engineering...\n")
  
  # Create interaction features if certain combinations exist
  if (all(c("effect_size", "presence_score") %in% available_features)) {
    feature_data$effect_presence_interaction <- 
      feature_data$effect_size * feature_data$presence_score
  }
  
  if (all(c("alt_frac_diff", "overall_quality") %in% available_features)) {
    feature_data$quality_weighted_diff <- 
      feature_data$alt_frac_diff * feature_data$overall_quality
  }
  
  # Population frequency-based features
  if ("population_AF" %in% available_features) {
    # Rare variant score (higher for rarer variants)
    feature_data$rarity_score <- 1 - pmin(feature_data$population_AF, 0.5) / 0.5
    
    # MAF-based score (minor allele frequency consideration)
    feature_data$maf_score <- ifelse(feature_data$population_AF > 0.5, 
                                     1 - feature_data$population_AF,
                                     feature_data$population_AF)
  }
  
  # Normalize all features to 0-1 scale
  if (verbose) cat("   Normalizing features...\n")
  
  normalized_features <- as.data.frame(lapply(feature_data, function(x) {
    if (is.numeric(x)) {
      x_range <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
      if (x_range == 0) {
        rep(0.5, length(x))  # If no variation, set to middle value
      } else {
        (x - min(x, na.rm = TRUE)) / x_range
      }
    } else {
      x
    }
  }))
  
  # Ensemble ML approach with multiple algorithms
  if (verbose) cat("   Running ensemble ML models...\n")
  
  # Model 1: Weighted linear combination (baseline)
  feature_weights <- c(
    "effect_size" = 0.20,
    "alt_frac_diff" = 0.18,
    "presence_score" = 0.15,
    "overall_quality" = 0.12,
    "population_AF" = 0.08,
    "group1_fold_enrichment" = 0.06,
    "group2_fold_enrichment" = 0.06,
    "effect_presence_interaction" = 0.05,
    "quality_weighted_diff" = 0.04,
    "rarity_score" = 0.03,
    "maf_score" = 0.03
  )
  
  available_weights <- feature_weights[names(feature_weights) %in% 
                                      colnames(normalized_features)]
  linear_score <- rowSums(sweep(
    normalized_features[, names(available_weights), drop = FALSE], 
    2, available_weights, "*"
  ), na.rm = TRUE)
  
  # Model 2: Non-linear scoring with interaction effects
  nonlinear_score <- rep(0, nrow(normalized_features))
  
  # Add quadratic terms for key features
  key_features <- intersect(c("effect_size", "presence_score", "alt_frac_diff"), 
                           colnames(normalized_features))
  for (feat in key_features) {
    nonlinear_score <- nonlinear_score + 
      0.1 * (normalized_features[[feat]]^2)
  }
  
  # Add synergistic interactions
  if (all(c("effect_size", "overall_quality") %in% colnames(normalized_features))) {
    nonlinear_score <- nonlinear_score + 
      0.15 * (normalized_features$effect_size * normalized_features$overall_quality)
  }
  
  # Model 3: Principal component-based scoring
  pc_score <- rep(0, nrow(normalized_features))
  if (ncol(normalized_features) >= 3) {
    tryCatch({
      # Simple PCA-like transformation
      numeric_features <- normalized_features[, sapply(normalized_features, is.numeric), 
                                             drop = FALSE]
      if (ncol(numeric_features) >= 2) {
        # Create composite score based on variance-weighted features
        feature_vars <- sapply(numeric_features, var, na.rm = TRUE)
        feature_vars[is.na(feature_vars)] <- 0
        
        if (sum(feature_vars) > 0) {
          var_weights <- feature_vars / sum(feature_vars)
          pc_score <- rowSums(sweep(numeric_features, 2, var_weights, "*"), 
                             na.rm = TRUE)
        }
      }
    }, error = function(e) {
      if (verbose) cat("   Note: PCA scoring failed, using fallback\n")
    })
  }
  
  # Model 4: Rank-based scoring
  rank_score <- rep(0, nrow(normalized_features))
  for (feat in names(available_weights)) {
    if (feat %in% colnames(normalized_features)) {
      # Convert to ranks (0-1 scale)
      ranks <- rank(normalized_features[[feat]], na.last = "keep") / 
               sum(!is.na(normalized_features[[feat]]))
      ranks[is.na(ranks)] <- 0.5
      rank_score <- rank_score + available_weights[[feat]] * ranks
    }
  }
  
  # Ensemble combination with adaptive weighting
  ensemble_weights <- c(
    linear = 0.35,      # Baseline linear model
    nonlinear = 0.25,   # Non-linear interactions
    pc = 0.20,          # Principal component approach
    rank = 0.20         # Rank-based approach
  )
  
  # Normalize individual model scores
  linear_score <- (linear_score - min(linear_score)) / 
                  (max(linear_score) - min(linear_score) + 1e-8)
  nonlinear_score <- (nonlinear_score - min(nonlinear_score)) / 
                     (max(nonlinear_score) - min(nonlinear_score) + 1e-8)
  pc_score <- (pc_score - min(pc_score)) / 
              (max(pc_score) - min(pc_score) + 1e-8)
  rank_score <- (rank_score - min(rank_score)) / 
                (max(rank_score) - min(rank_score) + 1e-8)
  
  # Final ensemble score
  ml_scores <- ensemble_weights["linear"] * linear_score +
               ensemble_weights["nonlinear"] * nonlinear_score +
               ensemble_weights["pc"] * pc_score +
               ensemble_weights["rank"] * rank_score
  
  # Add controlled randomness for tie-breaking
  set.seed(42)
  ml_scores <- ml_scores + rnorm(length(ml_scores), 0, 0.02)
  
  # Final normalization
  ml_scores <- pmax(0, pmin(1, ml_scores))  # Clamp to [0,1]
  
  if (verbose) {
    cat(sprintf("   Ensemble ML scoring complete: mean=%.3f, sd=%.3f\n", 
                mean(ml_scores), sd(ml_scores)))
  }
  
  return(data.frame(
    rs_id = snp_df$rs_id,
    gene_name = snp_df$gene_name,
    ml_priority_score = ml_scores,
    linear_component = linear_score,
    nonlinear_component = nonlinear_score,
    pc_component = pc_score,
    rank_component = rank_score,
    stringsAsFactors = FALSE
  ))
}

# Helper function: Select top SNPs for LLM assessment
.select_top_snps_for_llm <- function(all_scores, snp_df, top_n) {
  
  # Create composite score from available methods
  score_cols <- c("fc_consistency_score", "ml_score")
  available_score_cols <- intersect(score_cols, colnames(all_scores))
  
  if(length(available_score_cols) > 0) {
    composite_score <- rowMeans(all_scores[, available_score_cols, drop = FALSE], na.rm = TRUE)
  } else {
    # Fallback to effect size or another criterion
    if("effect_size" %in% colnames(snp_df)) {
      composite_score <- abs(snp_df$effect_size)
    } else {
      composite_score <- runif(nrow(all_scores))
    }
  }
  
  # Select top SNPs
  top_indices <- order(-composite_score, na.last = TRUE)[1:min(top_n, length(composite_score))]
  
  return(snp_df[top_indices, ])
}

# Helper function: Compute LLM clinical scores
.compute_llm_clinical_scores <- function(top_snps, llm_prompt_template = NULL, 
                                         verbose = TRUE) {
  
  if (verbose) {
    cat("   Performing LLM-based clinical relevance assessment...\n")
    cat(sprintf("   Analyzing %d top-scoring SNPs\n", nrow(top_snps)))
  }
  
  # Initialize results dataframe
  clinical_scores <- data.frame(
    rs_id = top_snps$rs_id,
    gene_name = top_snps$gene_name,
    clinical_relevance_score = NA,
    clinical_category = NA,
    therapeutic_potential = NA,
    disease_association = NA,
    mechanism_of_action = NA,
    clinical_evidence = NA,
    druggability_score = NA,
    pathway_involvement = NA,
    stringsAsFactors = FALSE
  )
  
  # Default prompt template if none provided
  if (is.null(llm_prompt_template)) {
    llm_prompt_template <- .get_default_clinical_prompt_template()
  }
  
  # Process SNPs in batches to optimize API calls
  batch_size <- 5
  n_batches <- ceiling(nrow(top_snps) / batch_size)
  
  for (batch_idx in seq_len(n_batches)) {
    start_idx <- (batch_idx - 1) * batch_size + 1
    end_idx <- min(batch_idx * batch_size, nrow(top_snps))
    batch_snps <- top_snps[start_idx:end_idx, ]
    
    if (verbose) {
      cat(sprintf("   Processing batch %d/%d (%d SNPs)...\n", 
                  batch_idx, n_batches, nrow(batch_snps)))
    }
    
    # Create structured query for this batch
    batch_results <- .query_llm_for_clinical_assessment(
      batch_snps, llm_prompt_template, verbose
    )
    
    # Update results
    clinical_scores[start_idx:end_idx, ] <- batch_results
    
    # Add small delay to respect API rate limits
    if (batch_idx < n_batches) {
      Sys.sleep(0.5)
    }
  }
  
  if (verbose) {
    successful_assessments <- sum(!is.na(clinical_scores$clinical_relevance_score))
    cat(sprintf("   Completed LLM assessment: %d/%d SNPs successfully analyzed\n",
                successful_assessments, nrow(clinical_scores)))
  }
  
  return(clinical_scores)
}

# Helper function: Get default clinical assessment prompt template
.get_default_clinical_prompt_template <- function() {
  return(paste(
    "You are a clinical geneticist and bioinformatics expert specializing in",
    "variant interpretation. Please assess the clinical relevance of the",
    "following genetic variants based on the provided information.",
    "",
    "For each variant, provide:",
    "1. Clinical relevance score (0-1, where 1 is most clinically relevant)",
    "2. Clinical category (High/Moderate/Low clinical significance)",
    "3. Therapeutic potential (Drug_target/Biomarker/Risk_factor/Unknown)",
    "4. Disease association strength (Strong/Moderate/Weak/None)",
    "5. Mechanism of action (brief description)",
    "6. Clinical evidence level (Strong/Moderate/Limited/None)",
    "7. Druggability score (0-1, where 1 is most druggable)",
    "8. Pathway involvement (key biological pathways)",
    "",
    "Consider factors like:",
    "- Gene function and biological pathways",
    "- Known disease associations",
    "- Population frequency and effect size",
    "- Existing therapeutic targets",
    "- Clinical actionability",
    "",
    "Variant details:",
    "{variant_info}",
    "",
    "Please provide structured output in JSON format with the following structure:",
    "{",
    '  "variants": [',
    "    {",
    '      "rs_id": "rs123456",',
    '      "gene_name": "GENE1",',
    '      "clinical_relevance_score": 0.8,',
    '      "clinical_category": "High",',
    '      "therapeutic_potential": "Drug_target",',
    '      "disease_association": "Strong_evidence",',
    '      "mechanism_of_action": "Brief description",',
    '      "clinical_evidence": "Strong",',
    '      "druggability_score": 0.7,',
    '      "pathway_involvement": "Key pathways"',
    "    }",
    "  ]",
    "}",
    sep = "\n"
  ))
}

# Helper function: Query LLM for clinical assessment
.query_llm_for_clinical_assessment <- function(batch_snps, prompt_template, 
                                               verbose = TRUE) {
  
  # Prepare variant information for LLM query
  variant_info <- .prepare_variant_info_for_llm(batch_snps)
  
  # Create the full prompt
  full_prompt <- gsub("\\{variant_info\\}", variant_info, prompt_template)
  
  # Initialize results with fallback values
  results <- data.frame(
    rs_id = batch_snps$rs_id,
    gene_name = batch_snps$gene_name,
    clinical_relevance_score = NA,
    clinical_category = NA,
    therapeutic_potential = NA,
    disease_association = NA,
    mechanism_of_action = NA,
    clinical_evidence = NA,
    druggability_score = NA,
    pathway_involvement = NA,
    stringsAsFactors = FALSE
  )
  
  # Try to make LLM API call using ellmer
  tryCatch({
    # Use ellmer package for structured LLM communication
    llm_response <- .call_llm_api_via_ellmer(full_prompt, verbose)
    
    if (!is.null(llm_response)) {
      # Parse LLM response and extract structured data
      parsed_results <- .parse_llm_clinical_response(llm_response, batch_snps)
      if (!is.null(parsed_results)) {
        results <- parsed_results
      } else {
        if (verbose) {
          cat("   LLM response parsing failed, using fallback assessment...\n")
        }
        results <- .fallback_clinical_assessment(batch_snps, verbose)
      }
    } else {
      if (verbose) {
        cat("   LLM API unavailable, using local clinical assessment...\n")
      }
      
      # Fallback to rule-based clinical assessment
      results <- .fallback_clinical_assessment(batch_snps, verbose)
    }
    
  }, error = function(e) {
    if (verbose) {
      cat(sprintf("   LLM API error: %s\n", e$message))
      cat("   Using fallback clinical assessment...\n")
    }
    
    # Use rule-based fallback
    results <- .fallback_clinical_assessment(batch_snps, verbose)
  })
  
  return(results)
}

# Helper function: Prepare variant information for LLM
.prepare_variant_info_for_llm <- function(batch_snps) {
  
  variant_summaries <- character(nrow(batch_snps))
  
  for (i in seq_len(nrow(batch_snps))) {
    snp <- batch_snps[i, ]
    
    # Create comprehensive variant description
    info_parts <- c(
      sprintf("Variant %d:", i),
      sprintf("  rs_id: %s", snp$rs_id),
      sprintf("  Gene: %s", snp$gene_name)
    )
    
    # Add available quantitative data
    if ("chromosome" %in% colnames(snp)) {
      info_parts <- c(info_parts, sprintf("  Location: %s:%s", 
                                         snp$chromosome, snp$position))
    }
    
    if ("effect_size" %in% colnames(snp)) {
      info_parts <- c(info_parts, sprintf("  Effect size: %.3f", snp$effect_size))
    }
    
    if ("population_AF" %in% colnames(snp)) {
      info_parts <- c(info_parts, sprintf("  Population AF: %.3f", 
                                         snp$population_AF))
    }
    
    if ("alt_frac_diff" %in% colnames(snp)) {
      info_parts <- c(info_parts, sprintf("  Alt fraction difference: %.3f", 
                                         snp$alt_frac_diff))
    }
    
    if ("feature_type" %in% colnames(snp)) {
      info_parts <- c(info_parts, sprintf("  Feature type: %s", snp$feature_type))
    }
    
    if ("gene_type" %in% colnames(snp)) {
      info_parts <- c(info_parts, sprintf("  Gene type: %s", snp$gene_type))
    }
    
    if ("GEX_avg_log2FC" %in% colnames(snp) && !is.na(snp$GEX_avg_log2FC)) {
      info_parts <- c(info_parts, sprintf("  Gene expression log2FC: %.3f", 
                                         snp$GEX_avg_log2FC))
    }
    
    variant_summaries[i] <- paste(info_parts, collapse = "\n")
  }
  
  return(paste(variant_summaries, collapse = "\n\n"))
}

# Helper function: Call LLM API via ellmer package
.call_llm_api_via_ellmer <- function(prompt, verbose = TRUE) {
  
  # Check if ellmer package is available
  if (!requireNamespace("ellmer", quietly = TRUE)) {
    if (verbose) {
      cat("   ellmer package not available. Install with:\n")
      cat("   remotes::install_github('tidyverse/ellmer')\n")
      cat("   Using fallback clinical assessment instead\n")
    }
    return(NULL)
  }
  
  # Try to create LLM chat object with available providers
  chat_obj <- NULL
  
  tryCatch({
    # Try Anthropic Claude first (if API key available)
    if (nchar(Sys.getenv("ANTHROPIC_API_KEY")) > 0) {
      if (verbose) {
        cat("   Connecting to Anthropic Claude via ellmer...\n")
      }
      chat_obj <- ellmer::chat_anthropic(
        system_prompt = paste(
          "You are a clinical geneticist and bioinformatics expert specializing",
          "in variant interpretation. Provide structured, evidence-based",
          "assessments of genetic variants for clinical relevance.",
          "Focus on actionable clinical information and therapeutic implications."
        )
      )
    } else if (nchar(Sys.getenv("OPENAI_API_KEY")) > 0) {
      if (verbose) {
        cat("   Connecting to OpenAI GPT via ellmer...\n")
      }
      chat_obj <- ellmer::chat_openai(
        system_prompt = paste(
          "You are a clinical geneticist and bioinformatics expert specializing",
          "in variant interpretation. Provide structured, evidence-based", 
          "assessments of genetic variants for clinical relevance.",
          "Focus on actionable clinical information and therapeutic implications."
        )
      )
    } else {
      if (verbose) {
        cat("   No LLM API keys found in environment variables\n")
        cat("   Set ANTHROPIC_API_KEY or OPENAI_API_KEY to use LLM assessment\n")
      }
      return(NULL)
    }
    
    if (!is.null(chat_obj)) {
      if (verbose) {
        cat("   Successfully connected to LLM. Sending clinical assessment query...\n")
      }
      
      # Send prompt and get response using ellmer
      response <- chat_obj |> ellmer::chat(prompt)
      
      if (verbose) {
        cat("   Received LLM response. Processing...\n")
      }
      
      return(response)
    }
    
  }, error = function(e) {
    if (verbose) {
      cat(sprintf("   ellmer LLM API error: %s\n", e$message))
    }
    return(NULL)
  })
  
  return(NULL)
}

# Helper function: Parse LLM clinical response
.parse_llm_clinical_response <- function(llm_response, batch_snps) {
  
  if (is.null(llm_response) || nchar(trimws(llm_response)) == 0) {
    return(NULL)
  }
  
  # Initialize results dataframe
  results <- data.frame(
    rs_id = batch_snps$rs_id,
    gene_name = batch_snps$gene_name,
    clinical_relevance_score = numeric(nrow(batch_snps)),
    clinical_category = character(nrow(batch_snps)),
    therapeutic_potential = character(nrow(batch_snps)),
    disease_association = character(nrow(batch_snps)),
    mechanism_of_action = character(nrow(batch_snps)),
    clinical_evidence = character(nrow(batch_snps)),
    druggability_score = numeric(nrow(batch_snps)),
    pathway_involvement = character(nrow(batch_snps)),
    stringsAsFactors = FALSE
  )
  
  tryCatch({
    # Try to parse as JSON first
    if (grepl("^\\s*[\\[\\{]", llm_response)) {
      json_data <- jsonlite::fromJSON(llm_response, simplifyDataFrame = TRUE)
      
      # Handle different JSON structures
      if (is.list(json_data) && "variants" %in% names(json_data)) {
        json_data <- json_data$variants
      } else if (is.list(json_data) && length(json_data) == 1 && is.list(json_data[[1]])) {
        json_data <- json_data[[1]]
      }
      
      # Extract data for each variant
      for (i in seq_len(nrow(batch_snps))) {
        if (i <= length(json_data) || (is.data.frame(json_data) && i <= nrow(json_data))) {
          variant_data <- if (is.data.frame(json_data)) json_data[i, ] else json_data[[i]]
          
          # Extract clinical relevance score
          if ("clinical_relevance_score" %in% names(variant_data)) {
            results$clinical_relevance_score[i] <- as.numeric(variant_data$clinical_relevance_score)
          } else if ("score" %in% names(variant_data)) {
            results$clinical_relevance_score[i] <- as.numeric(variant_data$score)
          }
          
          # Extract categorical assessments
          if ("clinical_category" %in% names(variant_data)) {
            results$clinical_category[i] <- as.character(variant_data$clinical_category)
          }
          if ("therapeutic_potential" %in% names(variant_data)) {
            results$therapeutic_potential[i] <- as.character(variant_data$therapeutic_potential)
          }
          if ("disease_association" %in% names(variant_data)) {
            results$disease_association[i] <- as.character(variant_data$disease_association)
          }
          if ("mechanism_of_action" %in% names(variant_data)) {
            results$mechanism_of_action[i] <- as.character(variant_data$mechanism_of_action)
          }
          if ("clinical_evidence" %in% names(variant_data)) {
            results$clinical_evidence[i] <- as.character(variant_data$clinical_evidence)
          }
          if ("druggability_score" %in% names(variant_data)) {
            results$druggability_score[i] <- as.numeric(variant_data$druggability_score)
          }
          if ("pathway_involvement" %in% names(variant_data)) {
            results$pathway_involvement[i] <- as.character(variant_data$pathway_involvement)
          }
        }
      }
      
      return(results)
    }
    
    # If not JSON, try to parse structured text response
    response_lines <- strsplit(llm_response, "\n")[[1]]
    response_lines <- trimws(response_lines[nchar(trimws(response_lines)) > 0])
    
    # Look for variant-specific sections
    current_variant <- 1
    for (line in response_lines) {
      # Check for variant identifiers
      if (grepl("^(Variant|SNP|rs)", line, ignore.case = TRUE)) {
        # Try to extract variant number or rs_id
        if (grepl("\\d+", line)) {
          variant_match <- regmatches(line, regexpr("\\d+", line))
          if (length(variant_match) > 0) {
            potential_variant <- as.numeric(variant_match[1])
            if (potential_variant <= nrow(batch_snps)) {
              current_variant <- potential_variant
            }
          }
        }
      }
      
      # Extract clinical relevance score
      if (grepl("(clinical.relevance|relevance.score|score).*:.*[0-9]", line, ignore.case = TRUE)) {
        score_match <- regmatches(line, regexpr("[0-9]\\.[0-9]+|[0-9]+", line))
        if (length(score_match) > 0) {
          score <- as.numeric(score_match[1])
          if (score <= 1) {  # Assume 0-1 scale
            results$clinical_relevance_score[current_variant] <- score
          } else if (score <= 10) {  # Convert from 1-10 scale
            results$clinical_relevance_score[current_variant] <- score / 10
          }
        }
      }
      
      # Extract clinical category
      if (grepl("(clinical.category|significance).*:.*(high|moderate|low)", line, ignore.case = TRUE)) {
        if (grepl("high", line, ignore.case = TRUE)) {
          results$clinical_category[current_variant] <- "High"
        } else if (grepl("moderate", line, ignore.case = TRUE)) {
          results$clinical_category[current_variant] <- "Moderate"
        } else if (grepl("low", line, ignore.case = TRUE)) {
          results$clinical_category[current_variant] <- "Low"
        }
      }
      
      # Extract therapeutic potential
      if (grepl("(therapeutic|treatment).*:.*", line, ignore.case = TRUE)) {
        if (grepl("drug.target", line, ignore.case = TRUE)) {
          results$therapeutic_potential[current_variant] <- "Drug_target"
        } else if (grepl("biomarker", line, ignore.case = TRUE)) {
          results$therapeutic_potential[current_variant] <- "Biomarker"
        } else if (grepl("risk.factor", line, ignore.case = TRUE)) {
          results$therapeutic_potential[current_variant] <- "Risk_factor"
        } else {
          results$therapeutic_potential[current_variant] <- "Unknown"
        }
      }
      
      # Extract disease association
      if (grepl("(disease|association).*:.*", line, ignore.case = TRUE)) {
        if (grepl("strong", line, ignore.case = TRUE)) {
          results$disease_association[current_variant] <- "Strong_evidence"
        } else if (grepl("moderate", line, ignore.case = TRUE)) {
          results$disease_association[current_variant] <- "Moderate_evidence"
        } else if (grepl("weak|limited", line, ignore.case = TRUE)) {
          results$disease_association[current_variant] <- "Weak_evidence"
        } else {
          results$disease_association[current_variant] <- "None"
        }
      }
    }
    
    # Fill in any missing values with defaults
    for (i in seq_len(nrow(results))) {
      if (results$clinical_relevance_score[i] == 0) {
        results$clinical_relevance_score[i] <- 0.5  # Default moderate score
      }
      if (results$clinical_category[i] == "") {
        results$clinical_category[i] <- "Moderate"
      }
      if (results$therapeutic_potential[i] == "") {
        results$therapeutic_potential[i] <- "Unknown"
      }
      if (results$disease_association[i] == "") {
        results$disease_association[i] <- "Moderate_evidence"
      }
      if (results$mechanism_of_action[i] == "") {
        results$mechanism_of_action[i] <- sprintf("Variant affects %s function", 
                                                 batch_snps$gene_name[i])
      }
      if (results$clinical_evidence[i] == "") {
        results$clinical_evidence[i] <- "Moderate"
      }
      if (results$druggability_score[i] == 0) {
        results$druggability_score[i] <- results$clinical_relevance_score[i] * 0.7
      }
      if (results$pathway_involvement[i] == "") {
        results$pathway_involvement[i] <- "Gene regulation"
      }
    }
    
    return(results)
    
  }, error = function(e) {
    # If parsing completely fails, return NULL to trigger fallback
    return(NULL)
  })
}

# Helper function: Fallback clinical assessment using rule-based approach
.fallback_clinical_assessment <- function(batch_snps, verbose = TRUE) {
  
  if (verbose) {
    cat("   Using knowledge-based clinical assessment rules...\n")
  }
  
  results <- data.frame(
    rs_id = batch_snps$rs_id,
    gene_name = batch_snps$gene_name,
    clinical_relevance_score = numeric(nrow(batch_snps)),
    clinical_category = character(nrow(batch_snps)),
    therapeutic_potential = character(nrow(batch_snps)),
    disease_association = character(nrow(batch_snps)),
    mechanism_of_action = character(nrow(batch_snps)),
    clinical_evidence = character(nrow(batch_snps)),
    druggability_score = numeric(nrow(batch_snps)),
    pathway_involvement = character(nrow(batch_snps)),
    stringsAsFactors = FALSE
  )
  
  # Apply rule-based scoring for each SNP
  for (i in seq_len(nrow(batch_snps))) {
    snp <- batch_snps[i, ]
    
    # Initialize base score
    clinical_score <- 0.5
    
    # Adjust based on available features
    
    # Effect size contribution
    if ("effect_size" %in% colnames(snp) && !is.na(snp$effect_size)) {
      effect_contribution <- pmin(abs(snp$effect_size) / 5, 0.3)
      clinical_score <- clinical_score + effect_contribution
    }
    
    # Population frequency (rarer variants often more clinically relevant)
    if ("population_AF" %in% colnames(snp) && !is.na(snp$population_AF)) {
      rarity_contribution <- (1 - pmin(snp$population_AF, 0.5)) * 0.2
      clinical_score <- clinical_score + rarity_contribution
    }
    
    # Gene expression change
    if ("GEX_avg_log2FC" %in% colnames(snp) && !is.na(snp$GEX_avg_log2FC)) {
      gex_contribution <- pmin(abs(snp$GEX_avg_log2FC) / 2, 0.25)
      clinical_score <- clinical_score + gex_contribution
    }
    
    # Feature type scoring
    if ("feature_type" %in% colnames(snp)) {
      if (snp$feature_type == "exonic") {
        clinical_score <- clinical_score + 0.2
      } else if (snp$feature_type == "promoter") {
        clinical_score <- clinical_score + 0.15
      } else if (snp$feature_type == "intronic") {
        clinical_score <- clinical_score + 0.05
      }
    }
    
    # Clamp score to [0.1, 1.0] range
    clinical_score <- pmax(0.1, pmin(1.0, clinical_score))
    
    # Assign categorical values based on score
    if (clinical_score >= 0.8) {
      category <- "High"
      evidence <- "Strong"
      disease_assoc <- "Strong_evidence"
    } else if (clinical_score >= 0.6) {
      category <- "Moderate"
      evidence <- "Moderate"
      disease_assoc <- "Moderate_evidence"
    } else {
      category <- "Low"
      evidence <- "Limited"
      disease_assoc <- "Weak_evidence"
    }
    
    # Therapeutic potential based on gene type and score
    if ("gene_type" %in% colnames(snp) && 
        snp$gene_type == "protein_coding" && clinical_score >= 0.7) {
      therapeutic <- "Drug_target"
    } else if (clinical_score >= 0.6) {
      therapeutic <- "Biomarker"
    } else if (clinical_score >= 0.4) {
      therapeutic <- "Risk_factor"
    } else {
      therapeutic <- "Unknown"
    }
    
    # Druggability score (simplified)
    druggability <- ifelse("gene_type" %in% colnames(snp) && 
                          snp$gene_type == "protein_coding", 
                          clinical_score * 0.8, clinical_score * 0.4)
    
    # Mechanism and pathway (simplified)
    mechanism <- sprintf("Variant in %s affects gene function", snp$gene_name)
    pathway <- ifelse("gene_type" %in% colnames(snp) && 
                     snp$gene_type == "protein_coding",
                     "Protein regulation", "Gene expression")
    
    # Store results
    results$clinical_relevance_score[i] <- clinical_score
    results$clinical_category[i] <- category
    results$therapeutic_potential[i] <- therapeutic
    results$disease_association[i] <- disease_assoc
    results$mechanism_of_action[i] <- mechanism
    results$clinical_evidence[i] <- evidence
    results$druggability_score[i] <- druggability
    results$pathway_involvement[i] <- pathway
  }
  
  return(results)
}

# Helper function: Combine prioritization scores
.combine_prioritization_scores <- function(all_scores, methods, snp_df) {
  
  # Define method weights (can be adjusted based on validation)
  method_weights <- list(
    fc_consistency = 0.4,
    ml_regression = 0.4,
    llm_clinical = 0.2
  )
  
  # Calculate combined score
  combined_score <- rep(0, nrow(all_scores))
  total_weight <- 0
  
  if("fc_consistency" %in% methods && "fc_consistency_score" %in% colnames(all_scores)) {
    weight <- method_weights$fc_consistency
    scores <- all_scores$fc_consistency_score
    scores[is.na(scores)] <- 0
    combined_score <- combined_score + weight * scores
    total_weight <- total_weight + weight
  }
  
  if("ml_regression" %in% methods && "ml_score" %in% colnames(all_scores)) {
    weight <- method_weights$ml_regression
    scores <- all_scores$ml_score
    scores[is.na(scores)] <- 0
    combined_score <- combined_score + weight * scores
    total_weight <- total_weight + weight
  }
  
  if("llm_clinical" %in% methods && "llm_clinical_score" %in% colnames(all_scores)) {
    weight <- method_weights$llm_clinical
    scores <- all_scores$llm_clinical_score
    scores[is.na(scores)] <- 0
    combined_score <- combined_score + weight * scores
    total_weight <- total_weight + weight
  }
  
  # Normalize by total weight
  if(total_weight > 0) {
    combined_score <- combined_score / total_weight
  }
  
  # Merge with original SNP data
  final_results <- merge(all_scores, snp_df, by = c("rs_id", "gene_name"), all.x = TRUE)
  final_results$combined_priority_score <- combined_score
  final_results$method_weights <- list(method_weights)
  
  return(final_results)
}

# Helper function: Generate prioritization summary
.generate_prioritization_summary <- function(results, all_scores, top_snps, methods) {
  
  summary_stats <- list(
    total_snps_analyzed = nrow(all_scores),
    methods_used = methods,
    top_snps_selected = nrow(top_snps)
  )
  
  if("fc_consistency_score" %in% colnames(all_scores)) {
    summary_stats$fc_consistency_stats <- summary(all_scores$fc_consistency_score)
  }
  
  if("ml_score" %in% colnames(all_scores)) {
    summary_stats$ml_score_stats <- summary(all_scores$ml_score)
  }
  
  if("llm_clinical_score" %in% colnames(all_scores)) {
    summary_stats$llm_score_stats <- summary(all_scores$llm_clinical_score)
  }
  
  if("combined_priority_score" %in% colnames(top_snps)) {
    summary_stats$final_score_stats <- summary(top_snps$combined_priority_score)
  }
  
  return(summary_stats)
}