#' @title Add a sample with SNP data to the variantCell project
#' @name addSampleData
#'
#' @description
#' Processes and adds a single-cell sample to the variantCell project by integrating metadata,
#' donor assignments from Vireo, and SNP data from CellSNP. This is the main function
#' for adding samples to a variantCell project and supports different input data types and
#' both transplant and non-transplant experimental designs.
#'
#' @param sample_id Character. Unique identifier for the sample.
#' @param vireo_path Character or NULL. Path to the Vireo donor_ids.tsv file. Required for
#'   transplant mode, can be NULL for non-transplant mode.
#' @param cellsnp_path Character. Path to the directory containing CellSNP output files.
#' @param cell_data Object. Cell data in one of three forms: Seurat object, SingleCellExperiment
#'   object, or a data frame with cell metadata (with cell identifiers as row names).
#' @param data_type Character. Type of cell_data object: "seurat", "sce", or "dataframe".
#' @param prefix_text Character. Text to prepend to cell identifiers in the Vireo data to match
#'   the cell barcodes in the input data.
#' @param donor_type Named character vector or NULL. Mapping between Vireo donor_id values and
#'   biological roles, e.g. c(donor0 = "Donor", donor1 = "Recipient"). Required for transplant mode.
#' @param non_transplant_mode Logical. Whether this sample is from a non-transplant experiment
#'   (TRUE) or a transplant experiment (FALSE).
#' @param min_cells Integer. Minimum number of cells with alternative allele required for a SNP
#'   to be included in the filtered dataset.
#' @param min_alt_frac Numeric. Minimum alternative allele fraction required when counting cells
#'   for the min_cells filter.
#' @param normalize Logical. Whether to calculate normalized SNP counts (TRUE) or not (FALSE).
#' @param scale.factor Numeric. Scaling factor for normalization. Only used if normalize=TRUE.
#' @param sample_metadata Data frame or NULL. Additional sample-level metadata.
#'
#' @return Invisibly returns self (the variantCell object) with the sample added to the samples list.
#'
#' @details
#' This function performs several steps:
#' 1. Processes the cell data and donor assignments based on data_type and mode
#' 2. Reads and processes the CellSNP data (SNP info and count matrices)
#' 3. Matches cell barcodes between the data sources
#' 4. Optionally normalizes the count data
#' 5. Filters SNPs based on minimum criteria
#' 6. Integrates metadata and stores the sample in the project
#'
#' The function supports two modes:
#' - Transplant mode: Uses Vireo to assign cells to donors in transplantation scenarios
#' - Non-transplant mode: Treats all cells as coming from a single donor
#'
#' @note
#' - For transplant mode, vireo_path and donor_type parameters are required
#' - For non-transplant mode, vireo_path can be NULL and donor_type defaults to c(donor0 = "donor0")
#' - The `normalize` parameter controls whether normalized expression values are calculated,
#'   which is useful for expression-based analyses
#' - The function expects specific file structure for CellSNP output: a base VCF file,
#'   and AD/DP count matrices in Matrix Market format
#'
#' @examples
#' \dontrun{
#' # Initialize a variantCell project
#' project <- variantCell$new()
#'
#' # Example 1: Add a sample in transplant mode using a Seurat object
#' project$addSampleData(
#'   sample_id = "transplant_sample1",
#'   vireo_path = "path/to/vireo/donor_ids.tsv",
#'   cellsnp_path = "path/to/cellsnp/output/",
#'   cell_data = seurat_object,
#'   data_type = "seurat",
#'   prefix_text = "Patient1_Sample1_",
#'   donor_type = c(donor0 = "Donor", donor1 = "Recipient"),
#'   normalize = TRUE
#' )
#'
#' # Example 2: Add a sample in non-transplant mode using a metadata data frame
#' project$addSampleData(
#'   sample_id = "non_transplant_sample1",
#'   vireo_path = NULL,
#'   cellsnp_path = "path/to/cellsnp/output/",
#'   cell_data = metadata_df,
#'   data_type = "dataframe",
#'   prefix_text = "Normal_Sample1_",
#'   non_transplant_mode = TRUE,
#'   normalize = TRUE
#' )
#' }
variantCell$set("public", "addSampleData", function(sample_id,
                                                    vireo_path = NULL,
                                                    cellsnp_path,
                                                    cell_data,
                                                    data_type = "seurat",   # Seurat, dataframe or sce
                                                    prefix_text,
                                                    donor_type = NULL,
                                                    non_transplant_mode = FALSE,
                                                    min_cells = 0,
                                                    min_alt_frac = 0,
                                                    normalize = TRUE,
                                                    scale.factor = 10000,
                                                    sample_metadata = NULL) {

  # Input validation
  if(non_transplant_mode) {
    # For non-transplant mode, vireo_path can be NULL
    if(is.null(donor_type)) {
      donor_type <- c(donor0 = "donor0")
    }
  } else {
    # For transplant mode, vireo_path is required
    if(is.null(vireo_path)) {
      stop("vireo_path must be provided when non_transplant_mode is FALSE")
    }

    # Validate donor_type
    if(is.null(donor_type) || !all(names(donor_type) %in% c("donor0", "donor1"))) {
      stop("donor_type must be named vector with 'donor0' and 'donor1' when non_transplant_mode is FALSE")
    }
  }

  # Helper function for matrix conversion
  convertToSparse <- function(mat, type = "CsparseMatrix") {
    if(!inherits(mat, "CsparseMatrix")) {
      mat <- as(mat, "CsparseMatrix")
    }
    return(mat)
  }

  if(!is.logical(normalize)) {
    stop("normalize must be TRUE or FALSE")
  }
  if(is.null(sample_id) || !is.character(sample_id)) {
    stop("sample_id must be a non-NULL character string")
  }
  if(sample_id %in% names(self$samples)) {
    stop(sprintf("Sample '%s' already exists in database", sample_id))
  }

  if(!is.null(vireo_path) && !file.exists(vireo_path)) {
    stop(sprintf("Vireo file not found: %s", vireo_path))
  }
  if(!file.exists(cellsnp_path)) {
    stop(sprintf("CellSNP directory not found: %s", cellsnp_path))
  }

  # Process cells and metadata differently based on data_type

  cat(sprintf("\nProcessing sample: %s\n", sample_id))

  if(non_transplant_mode) {
    # For non-transplant mode, skip Vireo processing
    # Just use the cell identifiers directly from the data
    if(data_type == "seurat") {
      # For Seurat objects
      matching_cells <- colnames(cell_data)
      metadata <- cell_data@meta.data
    } else if(data_type == "dataframe") {
      # For dataframes
      matching_cells <- rownames(cell_data)
      metadata <- cell_data
    } else if(data_type == "sce") {
      # For SingleCellExperiment
      matching_cells <- colnames(cell_data)
      metadata <- as.data.frame(colData(cell_data))
    } else {
      stop(sprintf("Unsupported data_type: %s", data_type))
    }

    # Create a simple donor data assignment (all cells to a single donor)
    donor_data <- data.frame(
      donor_id = rep("donor0", length(matching_cells)),
      row.names = matching_cells,
      stringsAsFactors = FALSE
    )
  } else {
    # Normal Vireo processing for transplant mode
    if(data_type == "seurat") {
      processed_data <- self$process_vireo_seurat(cell_data, vireo_path, prefix_text)
      matching_cells <- processed_data$matching_cells
      donor_data <- processed_data$donor_data
      seurat_obj <- processed_data$seurat_object
    } else if(data_type == "dataframe") {
      processed_data <- self$process_vireo_dataframe(cell_data, vireo_path, prefix_text)
      matching_cells <- processed_data$matching_cells
      donor_data <- processed_data$donor_data
      metadata <- processed_data$metadata
    } else if(data_type == "sce") {
      processed_data <- self$process_vireo_sce(cell_data, vireo_path, prefix_text)
      matching_cells <- processed_data$matching_cells
      donor_data <- processed_data$donor_data
      sce_obj <- processed_data$sce_object
    } else {
      stop(sprintf("Unsupported data_type: %s", data_type))
    }
  }

  if(length(matching_cells) == 0) {
    stop("No cells found for processing")
  }



  # Read CellSNP data
  cat("\nReading CellSNP data...")
  vcf_file <- file.path(cellsnp_path, "cellSNP.base.vcf")
  if(!file.exists(vcf_file)) {
    stop(sprintf("VCF file not found: %s", vcf_file))
  }

  vcf_lines <- readLines(vcf_file)
  header_lines <- grep("^#", vcf_lines)
  snp_data <- fread(text = vcf_lines[!startsWith(vcf_lines, "#")])
  colnames(snp_data)[1:8] <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")

  # Read matrices
  cat("\nReading allele count matrices...")
  ad_mtx <- convertToSparse(Matrix::readMM(file.path(cellsnp_path, "cellSNP.tag.AD.mtx")))
  dp_mtx <- convertToSparse(Matrix::readMM(file.path(cellsnp_path, "cellSNP.tag.DP.mtx")))

  # Process barcodes and find matching indices
  barcodes_file <- file.path(cellsnp_path, "cellSNP.samples.tsv")
  barcodes <- read.delim(barcodes_file, header=FALSE, stringsAsFactors=FALSE)[[1]]
  processed_barcodes <- sub(prefix_text, "", matching_cells)
  matching_indices <- match(processed_barcodes, barcodes)
  valid_matches <- !is.na(matching_indices)

  if(sum(valid_matches) == 0) {
    stop("No valid cell barcode matches found")
  }

  # Update matching cells and indices
  matching_cells <- matching_cells[valid_matches]
  matching_indices <- matching_indices[valid_matches]

  # Get matched matrices
  ad_matched <- convertToSparse(ad_mtx[, matching_indices])
  dp_matched <- convertToSparse(dp_mtx[, matching_indices])

  # Debug print
  cat("\nMatched Matrix Dimensions:")
  cat(sprintf("\nAD Matrix: %d rows x %d columns", nrow(ad_matched), ncol(ad_matched)))
  cat(sprintf("\nDP Matrix: %d rows x %d columns", nrow(dp_matched), ncol(dp_matched)))

  # Perform normalization before filtering if requested
  # Perform normalization if requested
  norm_result <- NULL
  if(normalize) {
    norm_result <- self$normalizeSnpCounts(convertToSparse(ad_matched),
                                           convertToSparse(dp_matched),
                                           scale.factor = scale.factor)
  }
  # Calculate alt fractions and apply filtering
  alt_fractions <- Matrix::sparseMatrix(i = integer(0), j = integer(0),
                                        dims = dim(dp_matched))
  valid_mask <- dp_matched > 0
  alt_fractions[valid_mask] <- ad_matched[valid_mask] / dp_matched[valid_mask]

  # Count cells meeting alt fraction threshold per SNP
  cells_passing <- Matrix::rowSums(alt_fractions >= min_alt_frac)

  # Calculate SNP metrics before filtering
  snp_metrics <- data.frame(
    snp_idx = 1:nrow(snp_data),
    total_depth = Matrix::rowSums(dp_matched),
    cells_passing = cells_passing,
    mean_dp_per_cell = Matrix::rowMeans(dp_matched),
    max_alt_frac = apply(alt_fractions, 1, max)
  )

  # Apply filter based on number of cells meeting alt fraction threshold
  keep_snps <- which(cells_passing >= min_cells)

  # Filter matrices and data
  ad_filtered <- convertToSparse(ad_matched[keep_snps, , drop=FALSE])
  dp_filtered <- convertToSparse(dp_matched[keep_snps, , drop=FALSE])
  snp_data_filtered <- snp_data[keep_snps, ]
  snp_metrics_filtered <- snp_metrics[keep_snps, ]

  # Filter normalized counts if they exist
  if(!is.null(norm_result)) {
    norm_result$norm_counts <- convertToSparse(norm_result$norm_counts[keep_snps, , drop=FALSE])
  }


  if(non_transplant_mode) {
    unified_metadata <- data.frame(
      cell_id = matching_cells,
      sample_id = sample_id,
      donor_id = donor_data[matching_cells, "donor_id"],
      donor_type = "donor0",  # Use a fixed value
      stringsAsFactors = FALSE
    )
  } else {
    unified_metadata <- data.frame(
      cell_id = matching_cells,
      sample_id = sample_id,
      donor_id = donor_data[matching_cells, "donor_id"],
      donor_type = donor_type[donor_data[matching_cells, "donor_id"]],
      stringsAsFactors = FALSE
    )
  }



  # Add source-specific metadata
  if(data_type == "seurat") {
    seurat_meta_cols <- setdiff(colnames(seurat_obj@meta.data), colnames(unified_metadata))
    unified_metadata[seurat_meta_cols] <- seurat_obj@meta.data[matching_cells, seurat_meta_cols]
  } else if(data_type == "dataframe") {
    df_meta_cols <- setdiff(colnames(metadata), colnames(unified_metadata))
    unified_metadata[df_meta_cols] <- metadata[matching_cells, df_meta_cols]
  } else if(data_type == "sce") {
    sce_meta_cols <- setdiff(colnames(colData(sce_obj)), colnames(unified_metadata))
    unified_metadata[sce_meta_cols] <- colData(sce_obj)[matching_cells, sce_meta_cols]
  }

  # Store sample data
  self$samples[[sample_id]] <- list(
    snp_data = snp_data_filtered,
    raw_metrics = list(
      ad_matrix = ad_filtered,  # raw AD for allele fractions
      dp_matrix = dp_filtered   # raw DP for allele fractions
    ),
    normalized_counts = if(normalize) norm_result$norm_counts else NULL,  # normalized DP for expression
    metadata = unified_metadata,
    cell_indices = matching_indices,
    snp_metrics = snp_metrics_filtered,
    filtering_info = list(
      min_cells = min_cells,
      min_alt_frac = min_alt_frac,
      total_snps = nrow(snp_data),
      passing_snps = length(keep_snps)
    ),
    normalization_info = if(normalize) list(
      scaling.factor = scale.factor,
      size_factors = norm_result$size_factors,
      metrics = norm_result$metrics
    ) else NULL
  )

  # Print summary
  cat("\n\nSample Addition Summary:")
  cat(sprintf("\nSample: %s", sample_id))
  cat(sprintf("\n - Total cells: %d", nrow(unified_metadata)))
  cat(sprintf("\n - Total SNPs: %d", nrow(snp_data)))
  cat(sprintf("\n - SNPs passing filters: %d", length(keep_snps)))
  cat(sprintf("\n - Filtering criteria:"))
  cat(sprintf("\n   * Minimum cells with alt allele: %d", min_cells))
  cat(sprintf("\n   * Minimum alt allele fraction: %.2f", min_alt_frac))
  cat(sprintf("\n - Donor distribution:"))
  print(table(unified_metadata$donor_type))


  invisible(self)
})
variantCell$set("public",  "normalizeSnpCounts", function(ad_matrix, dp_matrix, scale.factor = 10000) {

  # Helper function for matrix conversion - explicitly convert to dgCMatrix
  convertToSparse <- function(mat) {
    if(!inherits(mat, "dgCMatrix")) {
      # First convert to regular matrix then to dgCMatrix
      mat <- as(as.matrix(mat), "dgCMatrix")
    }
    return(mat)
  }

  # Convert input matrices to proper format
  dp_matrix <- convertToSparse(dp_matrix)

  # Calculate total depth per cell
  total_dp <- Matrix::colSums(dp_matrix)

  # Calculate scaling factor for each cell
  size_factors <- scale.factor / total_dp

  # Create sparse matrix for multiplication
  size_factor_mat <- Matrix::Diagonal(n = length(size_factors), x = size_factors)

  # Apply normalization and log transform, ensuring dgCMatrix output
  norm_dp <- convertToSparse(dp_matrix %*% size_factor_mat)
  norm_dp <- convertToSparse(log1p(norm_dp))

  # Calculate quality metrics
  metrics <- list(
    size_factors = summary(size_factors),
    pre_norm_depth = summary(Matrix::colSums(dp_matrix)),
    post_norm_depth = summary(Matrix::colSums(norm_dp)),
    scaling_factor = scale.factor
  )

  return(list(
    norm_counts = norm_dp,
    size_factors = size_factors,
    metrics = metrics
  ))
})
#' @title Build a unified SNP database from all samples
#' @name buildSNPDatabase
#'
#' @description
#' Integrates SNP data, annotations, and cell metadata from all samples in the variantCell project
#' into a unified database. This function creates merged sparse matrices for alternative allele (AD)
#' and depth (DP) counts across all cells, combines cell metadata, annotates SNPs with genomic
#' features, and calculates database-wide metrics. Optionally adds rsID numbers from VCF file reference.
#'
#'
#' @param add_rs_ids Logical. Whether to add rs# identifiers from a reference VCF file. Default: FALSE.
#' @param VCF_file_path Character. Path to reference VCF file (e.g., 1000 Genomes) for rs# annotation.
#'   Required if add_rs_ids = TRUE.
#'
#' @return Invisibly returns self (the variantCell object) with the unified SNP database constructed
#'   and stored in the snp_database field.
#'
#' @details
#' This function performs several key steps:
#' 1. Collects SNP information across all samples and identifies unique SNPs
#' 2. Retrieves genomic annotations for all SNPs (exonic, intronic, promoter, etc.)
#' 3. Combines metadata from all samples, handling missing columns appropriately
#' 4. Creates unified sparse matrices for AD and DP counts across all cells
#' 5. If available, also creates a matrix of normalized counts
#' 6. Calculates database-wide metrics for each SNP
#' 7. Generates a QC report with summary statistics
#' 8. Optionally, adds rsIDs
#'
#' The function handles the complexities of integrating data from multiple samples with
#' potentially different sets of SNPs and metadata columns. It manages matrix indexing,
#' column alignment, and other technical aspects needed to build a cohesive database.
#'
#' After running this function, all subsequent analyses (differential expression, plotting,
#' etc.) will use the unified database rather than individual sample data.
#'
#' @note
#' - This function must be called after adding all desired samples with `addSampleData()`
#' - The function requires at least one sample to be added
#' - SNP annotation may take significant time for large datasets
#' - The resulting database can use substantial memory for projects with many cells and SNPs
#'
#' @examples
#' \dontrun{
#' # Initialize a variantCell project
#' project <- variantCell$new()
#'
#' # Add samples
#' project$addSampleData(...)
#' project$addSampleData(...)
#'
#' # Build the unified SNP database
#' project$buildSNPDatabase()
#'
#' # Now the project is ready for analysis
#' project$setProjectIdentity("cell_type")
#' results <- project$findDESNPs(...)
#' }

variantCell$set("public", "buildSNPDatabase", function(add_rs_ids = FALSE, VCF_file_path = NULL) {
  cat("Building unified SNP database...\n")
  
  # Input validation
  if(length(self$samples) == 0) {
    stop("No samples added to database")
  }
  
  if(add_rs_ids) {
    if(is.null(VCF_file_path)) {
      stop("VCF_file_path must be provided when add_rs_ids = TRUE")
    }
    if(!file.exists(VCF_file_path)) {
      stop(sprintf("VCF file not found: %s", VCF_file_path))
    }
    cat(sprintf("rs# annotation will be added from: %s\n", VCF_file_path))
  }
  
  # First pass: collect SNP info and cell counts
  cat("\nCollecting SNP information across samples...")
  all_snps <- list()
  total_cells <- 0
  cell_ids <- character(0)
  has_normalized <- FALSE  # Track if any samples have normalized data
  
  # First loop: get SNP info and collect cell IDs
  for(sample_id in names(self$samples)) {
    sample <- self$samples[[sample_id]]
    snp_data <- sample$snp_data
    
    all_snps[[sample_id]] <- data.frame(
      CHROM = snp_data$CHROM,
      POS = snp_data$POS,
      REF = snp_data$REF,
      ALT = snp_data$ALT,
      sample = sample_id,
      snp_id = paste(snp_data$CHROM, snp_data$POS, snp_data$REF, snp_data$ALT, sep="_"),
      stringsAsFactors = FALSE
    )
    
    cell_ids <- c(cell_ids, sample$metadata$cell_id)
    total_cells <- total_cells + length(sample$metadata$cell_id)
    
    if(!is.null(sample$normalized_counts)) {
      has_normalized <- TRUE
    }
  }
  
  all_snps_df <- do.call(rbind, all_snps)
  unique_snps <- unique(all_snps_df[, c("CHROM", "POS", "REF", "ALT", "snp_id")])
  snp_sample_counts <- table(all_snps_df$snp_id)
  
  cat(sprintf("\nFound %d unique SNPs", nrow(unique_snps)))
  
  # Get annotations
  cat("\nAnnotating SNPs...")
  snp_annotations <- self$annotate_snps(unique_snps)
  
  # Initialize matrices
  combined_ad <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                                      dims = c(nrow(unique_snps), total_cells),
                                      dimnames = list(unique_snps$snp_id, cell_ids),
                                      giveCsparse = TRUE)
  combined_dp <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                                      dims = c(nrow(unique_snps), total_cells),
                                      dimnames = list(unique_snps$snp_id, cell_ids),
                                      giveCsparse = TRUE)
  
  combined_dp_norm <- if(has_normalized) {
    Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                         dims = c(nrow(unique_snps), total_cells),
                         dimnames = list(unique_snps$snp_id, cell_ids),
                         giveCsparse = TRUE)
  } else NULL
  
  # Create a list to store all metadata first
  metadata_list <- list()
  all_columns <- unique(unlist(lapply(self$samples, function(x) colnames(x$metadata))))
  
  cat("\nMetadata columns found across samples:")
  print(all_columns)
  
  # Process metadata for each sample
  for(sample_id in names(self$samples)) {
    sample <- self$samples[[sample_id]]
    current_meta <- sample$metadata
    
    # Add missing columns with NA
    missing_cols <- setdiff(all_columns, colnames(current_meta))
    if(length(missing_cols) > 0) {
      cat(sprintf("\nAdding missing columns for sample %s: %s",
                  sample_id, paste(missing_cols, collapse=", ")))
      for(col in missing_cols) {
        current_meta[[col]] <- NA
      }
    }
    
    # Ensure column order matches
    current_meta <- current_meta[, all_columns]
    metadata_list[[sample_id]] <- current_meta
  }
  
  # Combine metadata with error handling
  tryCatch({
    unified_metadata <- do.call(rbind, metadata_list)
    cat("\nSuccessfully combined metadata with dimensions:", dim(unified_metadata))
  }, error = function(e) {
    cat("\nError combining metadata:")
    cat("\nMetadata dimensions by sample:")
    for(sample_id in names(metadata_list)) {
      cat(sprintf("\n%s: %d x %d", sample_id,
                  nrow(metadata_list[[sample_id]]),
                  ncol(metadata_list[[sample_id]])))
    }
    stop("Failed to combine metadata: ", e$message)
  })
  
  # Process each sample's SNP data
  for(sample_id in names(self$samples)) {
    sample <- self$samples[[sample_id]]
    cat(sprintf("\nProcessing sample %s (%d cells)...",
                sample_id, ncol(sample$raw_metrics$ad_matrix)))
    
    # Create SNP mapping
    snp_data <- sample$snp_data
    snp_data$snp_id <- paste(snp_data$CHROM, snp_data$POS, snp_data$REF, snp_data$ALT, sep="_")
    snp_map <- match(unique_snps$snp_id, snp_data$snp_id)
    valid_snps <- !is.na(snp_map)
    
    if(sum(valid_snps) > 0) {
      # Get column indices
      col_idx <- match(sample$metadata$cell_id, colnames(combined_ad))
      
      if(any(is.na(col_idx))) {
        stop(sprintf("Some cell IDs from sample %s not found in combined matrix", sample_id))
      }
      
      # Update matrices
      tryCatch({
        combined_ad[which(valid_snps), col_idx] <- sample$raw_metrics$ad_matrix[snp_map[valid_snps], ]
        combined_dp[which(valid_snps), col_idx] <- sample$raw_metrics$dp_matrix[snp_map[valid_snps], ]
        
        if(!is.null(sample$normalized_counts)) {
          combined_dp_norm[which(valid_snps), col_idx] <- sample$normalized_counts[snp_map[valid_snps], ]
        }
      }, error = function(e) {
        stop(sprintf("Error updating matrices for sample %s: %s", sample_id, e$message))
      })
    }
  }
  
  # Calculate database-wide metrics
  cat("\nCalculating database-wide metrics...")
  db_metrics <- data.frame(
    snp_id = unique_snps$snp_id,
    chromosome = unique_snps$CHROM,
    position = unique_snps$POS,
    ref = unique_snps$REF,
    alt = unique_snps$ALT,
    feature_type = snp_annotations$feature_type,
    gene_name = snp_annotations$gene_name,
    mean_depth = Matrix::rowMeans(combined_dp),
    total_cells = Matrix::rowSums(combined_dp > 0),
    mean_alt_fraction = Matrix::rowSums(combined_ad) / Matrix::rowSums(combined_dp),
    n_samples = snp_sample_counts[unique_snps$snp_id],
    stringsAsFactors = FALSE
  )
  
  if(add_rs_ids) {
    cat("\n==================================================")
    cat("\nAdding rs# identifiers from reference VCF...")
    cat("\n==================================================")                       
    tryCatch({
      # Call the rs# matching function
      rs_ids <- private$match_rs_ids_internal(unique_snps, VCF_file_path)
      
      # Add rs# to snp_info, annotations, and metrics for consistency
      unique_snps$rs_id <- rs_ids          # ADD THIS LINE
      snp_annotations$rs_id <- rs_ids
      db_metrics$rs_id <- rs_ids
      
      # Report rs# matching results
      matched_rs <- sum(!is.na(rs_ids))
      total_snps <- length(rs_ids)
      cat(sprintf("\nrs# annotation results:"))
      cat(sprintf("\n - Total SNPs: %d", total_snps))
      cat(sprintf("\n - SNPs with rs# IDs: %d (%.1f%%)", matched_rs, (matched_rs/total_snps)*100))
      cat(sprintf("\n - SNPs without rs# IDs: %d (%.1f%%)", total_snps - matched_rs, ((total_snps - matched_rs)/total_snps)*100))
      
    }, error = function(e) {
      warning(sprintf("Failed to add rs# identifiers: %s", e$message))
      cat("\nContinuing without rs# annotation...")
      # Add empty rs_id columns
      unique_snps$rs_id <- NA_character_    # ADD THIS LINE
      snp_annotations$rs_id <- NA_character_
      db_metrics$rs_id <- NA_character_
    })
  } else {
    # Add empty rs_id columns for consistency
    unique_snps$rs_id <- NA_character_      # ADD THIS LINE
    snp_annotations$rs_id <- NA_character_
    db_metrics$rs_id <- NA_character_
  }
  
  # Store results
  self$snp_database <- list(
    ad_matrix = combined_ad,
    dp_matrix = combined_dp,
    dp_matrix_normalized = if(has_normalized) combined_dp_norm else NULL,
    cell_metadata = unified_metadata,
    snp_info = unique_snps,
    snp_annotations = snp_annotations,
    snp_metrics = db_metrics,
    qc_report = list(
      total_samples = length(self$samples),
      total_cells = total_cells,
      total_snps = nrow(unique_snps),
      has_normalized_data = has_normalized,
      has_rs_ids = add_rs_ids,
      rs_annotation_success = if(add_rs_ids) sum(!is.na(db_metrics$rs_id)) else 0,
      cells_per_sample = sapply(self$samples, function(x) ncol(x$raw_metrics$ad_matrix)),
      feature_distribution = table(snp_annotations$feature_type),
      biotype_distribution = table(snp_annotations$gene_type)
    )
  )
  
  # Print summary
  cat("\n\nSNP database summary:")
  cat(sprintf("\nTotal samples: %d", length(self$samples)))
  cat(sprintf("\nTotal cells: %d", total_cells))
  cat(sprintf("\nTotal unique SNPs: %d", nrow(unique_snps)))
  cat(sprintf("\nNormalized counts available: %s", if(has_normalized) "Yes" else "No"))
  if(add_rs_ids) {
    cat(sprintf("\nrs# identifiers: %d (%.1f%% coverage)", 
                self$snp_database$qc_report$rs_annotation_success,
                (self$snp_database$qc_report$rs_annotation_success/nrow(unique_snps))*100))
  }
  cat("\nFeature distribution:")
  print(table(snp_annotations$feature_type))
  
  invisible(self)
})

# Add the internal rs# matching helper function
variantCell$set("private", "match_rs_ids_internal", function(unique_snps, VCF_file_path) {
  
  require(data.table)
  
  # Create lookup using data.table for better performance
  cat("Creating position lookup...\n")
  project_positions <- unique(data.table(
    CHROM = as.character(unique_snps$CHROM),
    POS = as.integer(unique_snps$POS)
  ))
  
  cat(sprintf("Looking for %d unique positions\n", nrow(project_positions)))
  
  # Read VCF file
  cat("Reading VCF file...\n")
  
  if(grepl("\\.gz$", VCF_file_path)) {
    con <- gzfile(VCF_file_path, "r")
  } else {
    con <- file(VCF_file_path, "r")
  }
  
  # Read all lines
  all_lines <- readLines(con)
  close(con)
  
  # Remove header lines
  data_lines <- all_lines[!startsWith(all_lines, "#")]
  cat(sprintf("Found %d data lines in VCF\n", length(data_lines)))
  
  # Parse all VCF data at once using fread
  cat("Parsing VCF data...\n")
  vcf_data <- fread(text = paste(data_lines, collapse = "\n"))
  colnames(vcf_data) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
  
  # Ensure consistent data types
  vcf_data$CHROM <- as.character(gsub("^chr", "", vcf_data$CHROM))
  vcf_data$POS <- as.integer(vcf_data$POS)
  
  # Filter VCF to only positions of interest
  cat("Filtering VCF to positions of interest...\n")
  setkey(vcf_data, CHROM, POS)
  setkey(project_positions, CHROM, POS)
  
  # Inner join to get only matching positions
  matching_vcf <- vcf_data[project_positions, nomatch = 0]
  
  cat(sprintf("Found %d VCF records at target positions\n", nrow(matching_vcf)))
  
  if(nrow(matching_vcf) == 0) {
    cat("No matching positions found!\n")
    return(rep(NA_character_, nrow(unique_snps)))
  }
  
  # Handle multi-allelic sites
  cat("Processing multi-allelic sites...\n")
  expanded_vcf <- matching_vcf[rep(1:nrow(matching_vcf), sapply(strsplit(matching_vcf$ALT, ","), length))]
  alt_split <- unlist(strsplit(matching_vcf$ALT, ","))
  expanded_vcf$ALT <- alt_split
  
  # Create match keys for final matching
  expanded_vcf$match_key <- paste(expanded_vcf$CHROM, expanded_vcf$POS, 
                                  expanded_vcf$REF, expanded_vcf$ALT, sep = "_")
  unique_snps$match_key <- paste(unique_snps$CHROM, unique_snps$POS, 
                                 unique_snps$REF, unique_snps$ALT, sep = "_")
  
  # Final matching
  cat("Performing final SNP matching...\n")
  vcf_dt <- data.table(expanded_vcf)
  project_dt <- data.table(unique_snps)
  
  setkey(vcf_dt, match_key)
  setkey(project_dt, match_key)
  
  matched_data <- vcf_dt[project_dt, on = "match_key"]
  
  # Clean rs# identifiers
  matched_data$rs_id <- ifelse(matched_data$ID == "." | is.na(matched_data$ID), NA, matched_data$ID)
  
  # Return rs_ids in the same order as unique_snps
  rs_ids <- matched_data$rs_id[match(unique_snps$match_key, matched_data$match_key)]
  
  return(rs_ids)
})
variantCell$set("public", "annotate_snps", function(snp_info,
                                                    chunk_size = 5000,
                                                    promoter_upstream = 2000,
                                                    promoter_downstream = 200,
                                                    cached_features = NULL) {
  require(GenomicRanges)
  require(IRanges)
  require(EnsDb.Hsapiens.v86)  # Direct package instead of AnnotationHub
  
  cat("\nInitializing SNP annotation...")
  
  # Get EnsDb directly instead of through AnnotationHub
  if(is.null(cached_features)) {
    cat("\nAccessing EnsDb.Hsapiens.v86...")
    edb <- EnsDb.Hsapiens.v86  # Direct access
    
    # Get genomic features - rest of the code stays the same
    cat("\nRetrieving genomic features...")
    genes <- genes(edb)
    exons <- exons(edb)
    transcripts <- transcripts(edb)
    promoters <- promoters(genes, upstream = promoter_upstream,
                           downstream = promoter_downstream)
    
    # Store features for reuse
    cached_features <- list(
      genes = genes,
      exons = exons,
      transcripts = transcripts,
      promoters = promoters
    )
  } else {
    cat("\nUsing cached genomic features...")
    genes <- cached_features$genes
    exons <- cached_features$exons
    transcripts <- cached_features$transcripts
    promoters <- cached_features$promoters
  }
  
  # Create enhanced annotation data frame
  cat("\nInitializing annotation data frame...")
  annotations <- data.frame(
    snp_idx = 1:nrow(snp_info),
    chromosome = snp_info$CHROM,
    position = snp_info$POS,
    in_exon = FALSE,
    in_promoter = FALSE,
    in_gene = FALSE,
    in_transcript = FALSE,
    gene_id = NA_character_,
    gene_name = NA_character_,
    gene_type = NA_character_,
    transcript_ids = NA_character_,
    exon_ids = NA_character_,
    feature_type = "intergenic",
    promoter_distance = NA_integer_,
    strand = NA_character_,
    stringsAsFactors = FALSE
  )
  
  # Process in chunks
  n_chunks <- ceiling(nrow(snp_info) / chunk_size)
  cat(sprintf("\nProcessing %d SNPs in %d chunks...", nrow(snp_info), n_chunks))
  
  for(i in 1:n_chunks) {
    start_idx <- ((i-1) * chunk_size) + 1
    end_idx <- min(i * chunk_size, nrow(snp_info))
    
    cat(sprintf("\nProcessing chunk %d/%d (SNPs %d-%d)...",
                i, n_chunks, start_idx, end_idx))
    
    # Create GRanges for this chunk
    chunk_ranges <- GRanges(
      seqnames = snp_info$CHROM[start_idx:end_idx],
      ranges = IRanges(
        start = snp_info$POS[start_idx:end_idx],
        end = snp_info$POS[start_idx:end_idx]
      )
    )
    
    # Find overlaps for this chunk
    exon_overlaps <- suppressWarnings(findOverlaps(chunk_ranges, exons))
    promoter_overlaps <- suppressWarnings(findOverlaps(chunk_ranges, promoters))
    gene_overlaps <- suppressWarnings(findOverlaps(chunk_ranges, genes))
    transcript_overlaps <- suppressWarnings(findOverlaps(chunk_ranges, transcripts))
    
    # Process overlaps
    if(length(exon_overlaps) > 0) {
      chunk_indices <- queryHits(exon_overlaps) + start_idx - 1
      matched_exons <- exons[subjectHits(exon_overlaps)]
      
      annotations$in_exon[chunk_indices] <- TRUE
      annotations$feature_type[chunk_indices] <- "exonic"
      annotations$exon_ids[chunk_indices] <- mcols(matched_exons)$exon_id
    }
    
    if(length(promoter_overlaps) > 0) {
      chunk_indices <- queryHits(promoter_overlaps) + start_idx - 1
      matched_promoters <- promoters[subjectHits(promoter_overlaps)]
      
      annotations$in_promoter[chunk_indices] <- TRUE
      promoter_mask <- annotations$feature_type[chunk_indices] == "intergenic"
      annotations$feature_type[chunk_indices[promoter_mask]] <- "promoter"
      
      # Calculate distance to TSS
      tss_pos <- ifelse(strand(matched_promoters) == "+",
                        start(matched_promoters) + promoter_upstream,
                        end(matched_promoters) - promoter_upstream)
      annotations$promoter_distance[chunk_indices] <-
        abs(snp_info$POS[chunk_indices] - tss_pos)
    }
    
    if(length(transcript_overlaps) > 0) {
      chunk_indices <- queryHits(transcript_overlaps) + start_idx - 1
      matched_transcripts <- transcripts[subjectHits(transcript_overlaps)]
      
      annotations$in_transcript[chunk_indices] <- TRUE
      annotations$transcript_ids[chunk_indices] <- vapply(
        split(mcols(matched_transcripts)$tx_id, chunk_indices),
        paste, character(1), collapse=";"
      )
    }
    if(length(gene_overlaps) > 0) {
      chunk_indices <- queryHits(gene_overlaps) + start_idx - 1
      matched_genes <- genes[subjectHits(gene_overlaps)]
      
      annotations$in_gene[chunk_indices] <- TRUE
      annotations$gene_id[chunk_indices] <- mcols(matched_genes)$gene_id
      annotations$gene_name[chunk_indices] <- mcols(matched_genes)$symbol
      annotations$gene_type[chunk_indices] <- mcols(matched_genes)$gene_biotype
      annotations$strand[chunk_indices] <- as.character(strand(matched_genes))
      
      # Mark intronic SNPs
      intronic_mask <- annotations$feature_type[chunk_indices] == "intergenic" &
        !annotations$in_exon[chunk_indices] &
        !annotations$in_promoter[chunk_indices]
      annotations$feature_type[chunk_indices[intronic_mask]] <- "intronic"
    }
    
    if(i %% 10 == 0 || i == n_chunks) {
      cat("\nCurrent feature distribution:")
      print(table(annotations$feature_type))
    }
  }
  
  # Create detailed summary
  summary_stats <- list(
    total_snps = nrow(annotations),
    feature_counts = table(annotations$feature_type),
    gene_counts = sum(!is.na(annotations$gene_name)),
    unique_genes = length(unique(annotations$gene_name[!is.na(annotations$gene_name)])),
    biotype_distribution = table(annotations$gene_type[!is.na(annotations$gene_type)]),
    strand_distribution = table(annotations$strand[!is.na(annotations$strand)]),
    promoter_stats = if(any(!is.na(annotations$promoter_distance)))
      summary(annotations$promoter_distance[!is.na(annotations$promoter_distance)])
    else NULL,
    feature_overlap = list(
      exon_and_promoter = sum(annotations$in_exon & annotations$in_promoter),
      total_genic = sum(annotations$in_gene),
      total_exonic = sum(annotations$in_exon),
      total_intronic = sum(annotations$feature_type == "intronic"),
      total_promoter = sum(annotations$in_promoter)
    )
  )
  
  # Print final summary
  cat("\n\nAnnotation complete.")
  cat("\nFinal feature distribution:")
  print(summary_stats$feature_counts)
  cat(sprintf("\nUnique genes with SNPs: %d", summary_stats$unique_genes))
  cat("\nGene biotype distribution:")
  print(summary_stats$biotype_distribution)
  
  # Add summary as attribute
  attr(annotations, "summary") <- summary_stats
  # Add cached features as attribute for reuse
  attr(annotations, "cached_features") <- cached_features
  
  return(annotations)
})
