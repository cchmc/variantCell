#' @title aggregateByGroup: Aggregate SNP data by cell groups
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
#' \dontrun{
#' # Basic usage - aggregate by cell type
#' collapsed <- project$aggregateByGroup(
#'   group_by = "cell_type",
#'   use_normalized = TRUE
#' )
#'
#' # Analyze only donor cells with stricter filtering
#' donor_agg <- project$aggregateByGroup(
#'   group_by = "cell_type",
#'   donor_type = "Donor",
#'   min_cells_per_group = 5
#' )
#'
#' # Aggregate by disease status
#' disease_agg <- project$aggregateByGroup(
#'   group_by = "disease_status",
#'   use_normalized = TRUE
#' )
#' }
variantCell$set("public",  "aggregateByGroup", function(group_by,
                                                        donor_type = NULL,
                                                        min_cells_per_group = 3,
                                                        use_normalized = TRUE) {

  # Input validation
  if(!group_by %in% colnames(self$snp_database$cell_metadata)) {
    stop(sprintf("Column '%s' not found in metadata", group_by))
  }

  if(use_normalized && !("dp_matrix_normalized" %in% names(self$snp_database))) {
    warning("Normalized counts requested but not available. Using raw counts.")
    use_normalized <- FALSE
  }

  # Get data
  meta <- self$snp_database$cell_metadata
  ad_matrix <- self$snp_database$ad_matrix
  dp_matrix <- self$snp_database$dp_matrix
  dp_matrix_norm <- if(use_normalized) self$snp_database$dp_matrix_normalized else NULL

  # Check for non-transplant mode
  non_transplant_mode <- FALSE
  if(length(unique(meta$donor_type)) == 1 && unique(meta$donor_type)[1] == "donor0") {
    non_transplant_mode <- TRUE
    cat("\nNon-transplant mode detected (single donor type 'donor0')")
  }

  # Apply donor type filter if specified and not in non-transplant mode
  if(!is.null(donor_type) && !non_transplant_mode) {
    # Add explicit NA handling
    cells_to_use <- !is.na(meta$donor_type) & meta$donor_type == donor_type

    if(sum(cells_to_use) == 0) {
      stop(sprintf("No cells found for donor_type: %s", donor_type))
    }

    # Print some diagnostic information
    cat(sprintf("\nFiltering for donor type: %s", donor_type))
    cat(sprintf("\nTotal cells before filter: %d", nrow(meta)))
    cat(sprintf("\nCells passing filter: %d", sum(cells_to_use)))
    cat(sprintf("\nNA values in donor_type: %d", sum(is.na(meta$donor_type))))

    meta <- meta[cells_to_use, ]
    ad_matrix <- ad_matrix[, cells_to_use]
    dp_matrix <- dp_matrix[, cells_to_use]
    if(use_normalized) dp_matrix_norm <- dp_matrix_norm[, cells_to_use]
  } else if(!is.null(donor_type) && non_transplant_mode) {
    # In non-transplant mode but donor_type was specified
    cat(sprintf("\nIgnoring donor_type parameter in non-transplant mode (all cells are 'donor0')"))
  }

  # Get unique groups and sort them
  groups <- sort(unique(meta[[group_by]]))
  n_groups <- length(groups)

  # Initialize matrices
  collapsed_ad <- Matrix::Matrix(0, nrow=nrow(ad_matrix), ncol=n_groups, sparse=TRUE)
  collapsed_dp <- Matrix::Matrix(0, nrow=nrow(dp_matrix), ncol=n_groups, sparse=TRUE)
  if(use_normalized) {
    collapsed_dp_norm <- Matrix::Matrix(0, nrow=nrow(dp_matrix), ncol=n_groups, sparse=TRUE)
  }
  colnames(collapsed_ad) <- groups
  colnames(collapsed_dp) <- groups
  if(use_normalized) colnames(collapsed_dp_norm) <- groups

  # Create metadata
  collapsed_meta <- data.frame(
    group = groups,
    n_cells = sapply(groups, function(g) sum(meta[[group_by]] == g)),
    mean_depth = sapply(groups, function(g) {
      cells <- meta[[group_by]] == g
      if(sum(cells) > 0) {
        mean(Matrix::colSums(dp_matrix[, cells, drop=FALSE]))
      } else {
        0
      }
    }),
    stringsAsFactors = FALSE
  )

  # Process each group
  for(i in seq_along(groups)) {
    group <- groups[i]
    group_cells <- meta[[group_by]] == group

    if(collapsed_meta$n_cells[i] >= min_cells_per_group) {
      # Sum counts for all cells in group
      collapsed_ad[, i] <- Matrix::rowSums(ad_matrix[, group_cells, drop=FALSE])
      collapsed_dp[, i] <- Matrix::rowSums(dp_matrix[, group_cells, drop=FALSE])
      if(use_normalized) {
        collapsed_dp_norm[, i] <- Matrix::rowSums(dp_matrix_norm[, group_cells, drop=FALSE])
      }
    }
  }

  # Add filter status
  collapsed_meta$filter_status <- ifelse(
    collapsed_meta$n_cells >= min_cells_per_group,
    "included",
    "filtered_low_cells"
  )

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

  return(list(
    ad_matrix = collapsed_ad,
    dp_matrix = collapsed_dp,
    dp_matrix_normalized = if(use_normalized) collapsed_dp_norm else NULL,
    metadata = collapsed_meta,
    group_by = group_by,
    parameters = list(
      min_cells_per_group = min_cells_per_group,
      donor_type = donor_type,
      normalized = use_normalized,
      non_transplant_mode = non_transplant_mode
    ),
    snp_info = self$snp_database$snp_info,
    snp_annotations = self$snp_database$snp_annotations
  ))
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
                                                 include_rs_ids = TRUE) {  
  
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
#' @seealso
#' \code{\link{aggregateByGroup}} for preparing input data
#' \code{\link{findDESNPs}} for cell-level differential analysis
#' \code{\link{plotSNPs}} for visualizing the identified SNPs
#'
variantCell$set("public", "findSNPsByGroup", function(ident.1,
                                                      ident.2 = NULL,
                                                      aggregated_data,
                                                      min_depth = 10,          # Minimum depth in group with SNP
                                                      min_alt_frac = 0.2,      # Minimum alt allele fraction
                                                      max_alt_frac_other = 0.1, # Maximum alt fraction in other group
                                                      return_all = TRUE,
                                                      include_rs_ids = TRUE) {  # NEW PARAMETER
  
  # Input validation
  if(!all(c("ad_matrix", "dp_matrix", "metadata") %in% names(aggregated_data))) {
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
  
  # Get data
  ad_matrix <- aggregated_data$ad_matrix
  dp_matrix <- aggregated_data$dp_matrix
  meta <- aggregated_data$metadata
  
  # Check for non-transplant mode in the aggregated data
  non_transplant_mode <- FALSE
  if("parameters" %in% names(aggregated_data) &&
     "non_transplant_mode" %in% names(aggregated_data$parameters)) {
    non_transplant_mode <- aggregated_data$parameters$non_transplant_mode
  }
  
  # Create group masks and get cells
  group1_mask <- meta$group == ident.1 & meta$filter_status == "included"
  if(is.null(ident.2)) {
    # For "rest", combine all other groups
    group2_mask <- meta$filter_status == "included" & meta$group != ident.1
    group2_name <- "rest"
    cells1 <- sum(meta$n_cells[group1_mask])
    cells2 <- sum(meta$n_cells[group2_mask])
  } else {
    group2_mask <- meta$group == ident.2 & meta$filter_status == "included"
    group2_name <- ident.2
    cells1 <- meta$n_cells[group1_mask]
    cells2 <- meta$n_cells[group2_mask]
  }
  
  # Print initial summary
  cat("\n=== Starting SNP Analysis ===")
  cat(sprintf("\nComparing SNP presence between groups:"))
  cat(sprintf("\n%s: %d cells", ident.1, cells1))
  cat(sprintf("\n%s: %d cells", group2_name, cells2))
  if(non_transplant_mode) {
    cat("\nRunning in non-transplant mode (single donor)")
  } else if(!is.null(aggregated_data$parameters$donor_type)) {
    cat(sprintf("\nFiltered to donor type: %s", aggregated_data$parameters$donor_type))
  } else {
    cat("\nUsing all donor types")
  }
  
  if(has_rs_ids && include_rs_ids) {
    cat("\nrs# identifiers will be included in results")
  } else if(include_rs_ids) {
    cat("\nrs# identifiers requested but not available")
  }
  
  # Get aggregate counts for groups
  ad1 <- rowSums(ad_matrix[, group1_mask, drop=FALSE])
  dp1 <- rowSums(dp_matrix[, group1_mask, drop=FALSE])
  ad2 <- rowSums(ad_matrix[, group2_mask, drop=FALSE])
  dp2 <- rowSums(dp_matrix[, group2_mask, drop=FALSE])
  
  # Calculate allele fractions
  alt_frac1 <- ad1/dp1
  alt_frac2 <- ad2/dp2
  
  # Initialize results storage
  results <- vector("list", nrow(ad_matrix))
  results_count <- 0
  
  # Find SNPs meeting criteria
  for(i in seq_len(nrow(ad_matrix))) {
    # Look for SNPs present in group1 but not group2
    snp_in_group1 <- dp1[i] >= min_depth && alt_frac1[i] >= min_alt_frac && !is.na(alt_frac1[i])
    snp_not_in_group2 <- dp2[i] >= min_depth && alt_frac2[i] <= max_alt_frac_other && !is.na(alt_frac2[i])
    
    # Also look for SNPs present in group2 but not group1
    snp_in_group2 <- dp2[i] >= min_depth && alt_frac2[i] >= min_alt_frac && !is.na(alt_frac2[i])
    snp_not_in_group1 <- dp1[i] >= min_depth && alt_frac1[i] <= max_alt_frac_other && !is.na(alt_frac1[i])
    
    if((snp_in_group1 && snp_not_in_group2) || (snp_in_group2 && snp_not_in_group1)) {
      # Calculate presence score (0-1)
      # Higher score means stronger evidence for group-specific presence
      if(snp_in_group1 && snp_not_in_group2) {
        presence_score <- (alt_frac1[i] - alt_frac2[i]) *
          (dp1[i]/min_depth) *
          (1 - alt_frac2[i]/min_alt_frac)
        presence_score <- min(1, presence_score)
      } else {
        presence_score <- (alt_frac2[i] - alt_frac1[i]) *
          (dp2[i]/min_depth) *
          (1 - alt_frac1[i]/min_alt_frac)
        presence_score <- min(1, presence_score)
      }
      
      results_count <- results_count + 1
      
      # Create base results data frame
      result_row <- data.frame(
        snp_idx = i,
        chromosome = aggregated_data$snp_info$CHROM[i],
        position = aggregated_data$snp_info$POS[i],
        ref = aggregated_data$snp_info$REF[i],
        alt = aggregated_data$snp_info$ALT[i],
        feature_type = aggregated_data$snp_annotations$feature_type[i],
        gene_name = aggregated_data$snp_annotations$gene_name[i],
        gene_type = aggregated_data$snp_annotations$gene_type[i],
        
        # Coverage metrics
        depth_1 = dp1[i],
        depth_2 = dp2[i],
        alt_count_1 = ad1[i],
        alt_count_2 = ad2[i],
        alt_frac_1 = alt_frac1[i],
        alt_frac_2 = alt_frac2[i],
        n_cells_1 = cells1,
        n_cells_2 = cells2,
        
        # Presence score and classification
        presence_score = presence_score,
        presence = if(snp_in_group1 && snp_not_in_group2) sprintf("Present in %s", ident.1)
        else sprintf("Present in %s", group2_name),
        
        stringsAsFactors = FALSE
      )
      
      # Add rs# identifier if available and requested
      if(has_rs_ids && include_rs_ids) {
        result_row$rs_id <- self$snp_database$snp_info$rs_id[i]
      }
      
      # Store results
      results[[results_count]] <- result_row
    }
  }
  
  # Process final results
  if(results_count > 0) {
    all_results <- do.call(rbind, results[1:results_count])
    
    # Sort results by presence score
    all_results <- all_results[order(-all_results$presence_score), ]
    
    # Create summary
    summary <- list(
      total_tested = nrow(ad_matrix),
      passed_filters = results_count,
      present_in_group1 = sum(all_results$presence == sprintf("Present in %s", ident.1)),
      present_in_group2 = sum(all_results$presence == sprintf("Present in %s", group2_name)),
      parameters = list(
        min_depth = min_depth,
        min_alt_frac = min_alt_frac,
        max_alt_frac_other = max_alt_frac_other,
        non_transplant_mode = non_transplant_mode,
        include_rs_ids = include_rs_ids,
        rs_ids_available = has_rs_ids
      ),
      patterns = table(all_results$presence)
    )
    
    # Print summary
    cat("\n\n=== Analysis Summary ===")
    cat(sprintf("\nTotal SNPs tested: %d", summary$total_tested))
    cat(sprintf("\nSNPs passing filters: %d", summary$passed_filters))
    cat(sprintf("\n - Present in %s: %d", ident.1, summary$present_in_group1))
    cat(sprintf("\n - Present in %s: %d", group2_name, summary$present_in_group2))
    
    if(has_rs_ids && include_rs_ids) {
      rs_count <- sum(!is.na(all_results$rs_id))
      cat(sprintf("\n - SNPs with rs# identifiers: %d (%.1f%%)", 
                  rs_count, (rs_count/nrow(all_results))*100))
    }
    
    cat("\n\nPresence distribution:")
    print(summary$patterns)
    
    return(list(
      results = all_results,
      summary = summary
    ))
  }
  
  cat("\nNo SNPs meeting criteria found")
  return(NULL)
})