# Improved plotSNPHeatmap function with efficiency optimizations

plotSNPHeatmap_improved <- function(
  genes = NULL,
  snp_indices = NULL,
  group.by,
  split.by = NULL,
  min_alt_frac = 0.2,
  scale_data = TRUE,
  max_scale = 2,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 8,
  exclude_empty = TRUE,
  normalize_by_cells = TRUE,
  data_out = FALSE,
  use_rs_ids = TRUE,           
  rs_id_format = "mixed"       
) {
  
  # IMPROVEMENT 1: Hybrid color generation - dittoSeq up to 40, then HSV extension
  generate_smart_colors <- function(n) {
    if(n <= 0) return(character(0))
    
    if(n <= 40) {
      # Use dittoSeq colors for optimal distinction (up to 40)
      if(requireNamespace("dittoSeq", quietly = TRUE)) {
        return(dittoSeq::dittoColors()[1:n])
      } else {
        # Fallback if dittoSeq not available - first 40 colors manually
        ditto_fallback <- c(
          "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999",
          "#B3770A", "#4A8BC2", "#007A5E", "#C4C142", "#005A8E", "#B04700", "#A86B85", "#737373",
          "#F2B733", "#7AC4ED", "#33B08A", "#F7ED75", "#3392C4", "#E8824D", "#D9A5C4", "#B3B3B3",
          "#F5C866", "#A1D4F0", "#66C1A7", "#F9F0A8", "#66A8D1", "#ED9966", "#E4BDD4", "#CCCCCC",
          "#9E6B00", "#3D6E99", "#005E4A", "#A3A335", "#004671", "#8F3700", "#8A5569", "#595959"
        )
        return(ditto_fallback[1:n])
      }
    } else {
      # For >40 colors: dittoSeq base + HSV extension
      if(requireNamespace("dittoSeq", quietly = TRUE)) {
        base_colors <- dittoSeq::dittoColors()[1:40]
      } else {
        base_colors <- c(
          "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999",
          "#B3770A", "#4A8BC2", "#007A5E", "#C4C142", "#005A8E", "#B04700", "#A86B85", "#737373",
          "#F2B733", "#7AC4ED", "#33B08A", "#F7ED75", "#3392C4", "#E8824D", "#D9A5C4", "#B3B3B3",
          "#F5C866", "#A1D4F0", "#66C1A7", "#F9F0A8", "#66A8D1", "#ED9966", "#E4BDD4", "#CCCCCC",
          "#9E6B00", "#3D6E99", "#005E4A", "#A3A335", "#004671", "#8F3700", "#8A5569", "#595959"
        )
      }
      
      additional_needed <- n - 40
      
      # Generate additional colors using HSV for maximum distinction
      hues <- seq(0.05, 0.95, length.out = additional_needed)
      saturations <- rep(c(0.8, 0.6, 0.9), length.out = additional_needed)
      values <- rep(c(0.7, 0.9, 0.6), length.out = additional_needed)
      
      additional_colors <- hsv(h = hues, s = saturations, v = values)
      
      return(c(base_colors, additional_colors))
    }
  }
  
  if(is.null(self$snp_database)) {
    stop("SNP database not found. Run buildSNPDatabase first.")
  }
  
  if(is.null(genes) && is.null(snp_indices)) {
    stop("Must provide either genes or snp_indices")
  }
  
  if(!group.by %in% colnames(self$snp_database$cell_metadata)) {
    stop(sprintf("group.by column '%s' not found in metadata", group.by))
  }
  
  if(!is.null(split.by) && !split.by %in% colnames(self$snp_database$cell_metadata)) {
    stop(sprintf("split.by column '%s' not found in metadata", split.by))
  }
  
  # Get SNP indices if genes provided
  if(!is.null(genes)) {
    snp_indices <- which(self$snp_database$snp_annotations$gene_name %in% genes)
    if(length(snp_indices) == 0) {
      stop("No SNPs found for provided genes")
    }
  }
  
  # IMPROVEMENT 2: Subset matrices ONCE at the beginning for efficiency
  cat(sprintf("Subsetting matrices for %d SNPs...\n", length(snp_indices)))
  
  # Subset all matrices upfront
  dp_subset <- self$snp_database$dp_matrix[snp_indices, , drop = FALSE]
  ad_subset <- self$snp_database$ad_matrix[snp_indices, , drop = FALSE]
  
  # Also subset normalized matrix if available
  dp_norm_subset <- NULL
  if(!is.null(self$snp_database$dp_matrix_normalized)) {
    dp_norm_subset <- self$snp_database$dp_matrix_normalized[snp_indices, , drop = FALSE]
  }
  
  # Get metadata
  meta <- self$snp_database$cell_metadata
  
  # IMPROVEMENT 3: Optimized group stats calculation using subsetted matrices
  calculate_group_stats_optimized <- function(snp_row_idx, cell_mask) {
    # Use pre-subsetted matrices (snp_row_idx refers to row in subset, not original matrix)
    dp_vals <- dp_subset[snp_row_idx, cell_mask]
    ad_vals <- ad_subset[snp_row_idx, cell_mask]
    
    # Calculate alt fractions where depth > 0
    valid_mask <- dp_vals > 0
    if(sum(valid_mask) == 0) {
      return(list(mean_expr = 0, n_cells = 0, n_expr_cells = 0))
    }
    
    alt_frac <- rep(0, length(dp_vals))
    alt_frac[valid_mask] <- ad_vals[valid_mask] / dp_vals[valid_mask]
    
    # Get cells meeting alt fraction threshold
    alt_mask <- alt_frac >= min_alt_frac
    n_expr_cells <- sum(alt_mask)
    if(n_expr_cells == 0) {
      return(list(mean_expr = 0, n_cells = sum(cell_mask), n_expr_cells = 0))
    }
    
    # Get expression values from appropriate matrix
    if(is.null(dp_norm_subset)) {
      expr_vals <- dp_vals[alt_mask]
    } else {
      expr_vals <- dp_norm_subset[snp_row_idx, cell_mask][alt_mask]
    }
    
    # Normalize by cell count if requested
    if(normalize_by_cells) {
      mean_expr <- mean(expr_vals) * n_expr_cells / sum(cell_mask)
    } else {
      mean_expr <- mean(expr_vals)
    }
    
    return(list(
      mean_expr = mean_expr,
      n_cells = sum(cell_mask),
      n_expr_cells = n_expr_cells
    ))
  }
  
  # Create grouping factors
  if(!is.null(split.by)) {
    group_factor <- paste(meta[[group.by]], meta[[split.by]], sep = "_")
  } else {
    group_factor <- meta[[group.by]]
  }
  
  unique_groups <- unique(group_factor)
  n_groups <- length(unique_groups)
  n_snps <- length(snp_indices)
  
  cat(sprintf("Processing %d SNPs across %d groups...\n", n_snps, n_groups))
  
  # Create expression matrix
  group_matrix <- matrix(0, nrow = n_snps, ncol = n_groups)
  colnames(group_matrix) <- unique_groups
  
  # IMPROVEMENT 4: Process with progress indication and optimized loop
  for(i in 1:n_snps) {
    if(i %% 10 == 0) cat(sprintf("  Processing SNP %d/%d\n", i, n_snps))
    
    for(j in 1:n_groups) {
      cell_mask <- group_factor == unique_groups[j]
      if(sum(cell_mask) > 0) {
        stats <- calculate_group_stats_optimized(i, cell_mask)  # i is row index in subset
        group_matrix[i, j] <- stats$mean_expr
      }
    }
  }
  
  # Exclude empty rows/columns if requested
  if(exclude_empty) {
    row_sums <- rowSums(group_matrix)
    col_sums <- colSums(group_matrix)
    
    keep_rows <- row_sums > 0
    keep_cols <- col_sums > 0
    
    if(sum(keep_rows) == 0 || sum(keep_cols) == 0) {
      stop("No data remaining after filtering empty rows/columns")
    }
    
    group_matrix <- group_matrix[keep_rows, keep_cols, drop = FALSE]
    snp_indices <- snp_indices[keep_rows]  # Update SNP indices to match
  }
  
  # Check if rs# IDs are available and create row labels
  has_rs_ids <- FALSE
  if(use_rs_ids) {
    if(!is.null(self$snp_database) && 
       "snp_info" %in% names(self$snp_database) &&
       "rs_id" %in% colnames(self$snp_database$snp_info)) {
      has_rs_ids <- !all(is.na(self$snp_database$snp_info$rs_id))
    }
    
    if(!has_rs_ids && use_rs_ids) {
      cat("Warning: rs# identifiers requested but not available. Using chromosome:position format.\n")
    }
  }
  
  # Create row labels
  if(has_rs_ids && use_rs_ids) {
    rs_ids <- self$snp_database$snp_info$rs_id[snp_indices]
    
    if(rs_id_format == "rs_only") {
      row_labels <- ifelse(!is.na(rs_ids), rs_ids, 
                          paste(self$snp_database$snp_info$CHROM[snp_indices],
                               self$snp_database$snp_info$POS[snp_indices], sep=":"))
    } else if(rs_id_format == "chr_pos") {
      row_labels <- paste(self$snp_database$snp_info$CHROM[snp_indices],
                         self$snp_database$snp_info$POS[snp_indices], sep=":")
    } else { # "mixed" format (default)
      row_labels <- ifelse(!is.na(rs_ids), rs_ids,
                          paste(self$snp_database$snp_info$CHROM[snp_indices],
                               self$snp_database$snp_info$POS[snp_indices],
                               self$snp_database$snp_info$REF[snp_indices],
                               self$snp_database$snp_info$ALT[snp_indices], sep="_"))
    }
  } else {
    row_labels <- paste(self$snp_database$snp_info$CHROM[snp_indices],
                       self$snp_database$snp_info$POS[snp_indices],
                       self$snp_database$snp_info$REF[snp_indices],
                       self$snp_database$snp_info$ALT[snp_indices], sep="_")
  }
  
  rownames(group_matrix) <- row_labels
  
  # Scale data if requested
  if(scale_data) {
    scaled_matrix <- t(scale(t(group_matrix)))
    scaled_matrix[is.na(scaled_matrix)] <- 0
    scaled_matrix <- pmin(pmax(scaled_matrix, -max_scale), max_scale)
  } else {
    scaled_matrix <- group_matrix
  }
  
  # Return data if requested
  if(data_out) {
    snp_info_out <- data.frame(
      snp_idx = snp_indices,
      snp_id = row_labels,
      stringsAsFactors = FALSE
    )
    
    return(list(
      expression_matrix = scaled_matrix,
      snp_info = snp_info_out,
      group_names = colnames(scaled_matrix)
    ))
  }
  
  # IMPROVEMENT 5: Generate colors dynamically based on actual data
  gene_names <- self$snp_database$snp_annotations$gene_name[snp_indices]
  gene_types <- self$snp_database$snp_annotations$feature_type[snp_indices]
  
  unique_genes <- unique(gene_names)
  unique_features <- unique(gene_types)
  
  # Generate exactly the number of colors needed
  gene_colors <- setNames(generate_smart_colors(length(unique_genes)), unique_genes)
  feature_colors <- setNames(generate_smart_colors(length(unique_features)), unique_features)
  
  cat(sprintf("Generated %d gene colors and %d feature colors\n", 
              length(gene_colors), length(feature_colors)))
  
  # Create row split factor for genes
  row_split <- factor(gene_names)
  
  # Build annotations
  left_annotation <- ComplexHeatmap::rowAnnotation(
    Gene = gene_names,
    Feature = gene_types,
    col = list(
      Gene = gene_colors,
      Feature = feature_colors
    ),
    annotation_width = unit(c(1, 1), "cm"),
    annotation_name_gp = gpar(fontsize = 8)
  )
  
  # Create heatmap
  ht <- ComplexHeatmap::Heatmap(
    scaled_matrix,
    name = if(scale_data) "Scaled\nExpression" else "Expression",
    
    # Row parameters
    cluster_rows = cluster_rows,
    show_row_names = show_rownames,
    row_names_gp = gpar(fontsize = fontsize_row),
    row_split = if(length(unique_genes) > 1) row_split else NULL,
    
    # Column parameters  
    cluster_columns = cluster_cols,
    show_column_names = show_colnames,
    column_names_gp = gpar(fontsize = fontsize_col),
    
    # Annotations
    left_annotation = left_annotation,
    
    # Color scheme
    col = circlize::colorRamp2(
      c(-max_scale, 0, max_scale),
      c("blue", "white", "red")
    )
  )
  
  return(ht)
}

# Usage examples:
# 
# # Basic usage - automatically generates right number of colors
# heatmap <- project$plotSNPHeatmap(genes = c("GENE1", "GENE2", "GENE3"), group.by = "cell_type")
# 
# # Large gene list - efficiently handles many SNPs
# many_genes <- c("BRCA1", "TP53", "EGFR", "KRAS", "PIK3CA") # ... up to 100+ genes
# heatmap <- project$plotSNPHeatmap(genes = many_genes, group.by = "condition")
# 
# # The function will:
# # 1. Generate exactly the right number of colors (no waste, no shortage)
# # 2. Subset matrices once at the beginning (much faster)
# # 3. Process efficiently with progress updates