#' @title plotSNPs "Visualize SNP Distribution for a Gene"
#' @name plotSNPs
#'
#' @description
#' Creates genomic visualizations of SNP distribution for a specified gene, showing alternative allele
#' frequencies across different cell populations. This function generates interactive plots
#' displaying SNPs along the gene's coordinates with various grouping options.
#'
#' @param gene Character. Name of the gene to visualize.
#' @param group.by Character, optional. Primary grouping variable from metadata.
#'               If NULL, all cells are treated as one group.
#' @param split.by Character, optional. Secondary grouping variable for within-group comparisons.
#'               Useful for donor/recipient or condition comparisons within cell types.
#' @param idents Vector, optional. Specific identity values to include. Must be values from group.by column.
#' @param min_depth Integer. Minimum read depth required to include a SNP in visualization.
#' @param min_cells Integer. Minimum cells per group required to include a group in visualization.
#' @param min_alt_frac Numeric between 0 and 1. Minimum alternative allele fraction required to include a SNP.
#'                    Set to 0 to show all SNPs regardless of alt fraction.
#' @param flank_size Integer. Size of flanking regions (in bp) to include around the gene.
#' @param plot_density Logical. Whether to include density distribution plots below the main visualization.
#' @param data_out Logical. If TRUE, returns a data frame with SNP data instead of plots.
#' @param use_normalized Logical. Whether to use normalized read depth values.
#' @param color_scheme Vector. Named vector with "low" and "high" colors for alt fraction gradient.
#' @param point_size_range Numeric vector of length 2. Range of point sizes (minimum and maximum) for depth representation.
#'
#' @return
#' If data_out = FALSE (default): Plot grid with main SNP visualization and optional density plots.
#' If data_out = TRUE: Data frame containing SNP information for the gene including positions,
#' alternative allele frequencies, depths, and other metrics for each group.
#'
#' @details
#' This function creates a comprehensive visualization of SNPs within a gene region, showing:
#'
#' 1. A main plot with:
#'    - Gene structure with exons displayed as black rectangles
#'    - SNP positions represented as points along genomic coordinates
#'    - Alternative allele frequencies encoded by color
#'    - Read depth encoded by point size
#'    - Groups and splits organized vertically
#'
#' 2. Optional density plots showing:
#'    - Alternative allele fraction distribution
#'    - Read depth distribution (log10 scale)
#'
#' The function applies filters to ensure only high-quality data is displayed. SNPs with low
#' depth or low alternative allele frequency can be excluded. Groups with insufficient cells
#' are also not displayed.
#'
#' @note
#' - Requires gene annotation data from AnnotationHub (automatically accessed)
#' - May take longer to generate for genes with many SNPs or when using many groups
#' - Setting min_alt_frac = 0 shows all SNPs but may make plots very busy
#' - The data_out parameter is useful for further custom analysis or visualization
#'
#' @examples
#' # Basic visualization of a gene
#'
#'
#' \dontrun{
#'
#' # Split by cell type and donor type
#' project$plotSNPs(
#'   gene = "IL7R",
#'   group.by = "cell_type",
#'   split.by = "donor_type",
#'   min_alt_frac = 0.1
#' )
#'
#' # Focus on specific cell types
#' project$plotSNPs(
#'   gene = "PDCD1",
#'   group.by = "cell_type",
#'   idents = c("CD4_T", "CD8_T", "Treg"),
#'   flank_size = 10000
#' )
#'
#' # Return data frame instead of plot for custom visualization
#' snp_data <- project$plotSNPs(
#'   gene = "HLA-DRB1",
#'   group.by = "condition",
#'   data_out = TRUE
#' )
#'}
#' @seealso
#' \code{\link{findDESNPs}} for identifying differentially expressed SNPs
#' \code{\link{plotSNPHeatmap}} for heatmap visualization of multiple genes
variantCell$set("public",  "plotSNPs", function(gene,
                                                group.by = NULL,      # Primary grouping variable
                                                split.by = NULL,      # Optional splitting variable
                                                idents = NULL,        # Specific identities to include
                                                min_depth = 10,       # Minimum read depth
                                                min_cells = 3,        # Minimum cells per group
                                                min_alt_frac = 0.2,   # Minimum alt allele fraction
                                                flank_size = 5000,    # Genomic region padding
                                                plot_density = TRUE,  # Whether to plot density plots
                                                data_out = FALSE,     # Whether to return data frame instead of plot
                                                use_normalized = FALSE, # Whether to use normalized depth data
                                                color_scheme = c("low" = "blue", "high" = "red"),
                                                point_size_range = c(2, 8)) {


  require(ggplot2)
  require(cowplot)
  require(AnnotationHub)

  # Input validation
  if(is.null(self$snp_database)) {
    stop("SNP database not found. Run buildSNPDatabase first.")
  }

  # Parameter validation
  params <- list(
    min_depth = min_depth,
    min_cells = min_cells,
    min_alt_frac = min_alt_frac,
    flank_size = flank_size
  )

  # Get metadata and clean
  meta <- self$snp_database$cell_metadata

  # Remove NAs from donor_type if it's the split variable
  if(!is.null(split.by) && split.by == "donor_type") {
    meta <- meta[!is.na(meta[[split.by]]),]
  }

  # Validate grouping parameters
  if(!is.null(group.by) && !group.by %in% colnames(meta)) {
    stop(sprintf("group.by column '%s' not found in metadata", group.by))
  }
  if(!is.null(split.by) && !split.by %in% colnames(meta)) {
    stop(sprintf("split.by column '%s' not found in metadata", split.by))
  }
  # Check for normalized data if requested
  if(use_normalized && is.null(self$snp_database$dp_matrix_normalized)) {
    warning("Normalized counts requested but not available. Using raw counts.")
    use_normalized <- FALSE
  }

  # Set up factors for proper ordering
  if(!is.null(group.by)) {
    meta[[group.by]] <- factor(meta[[group.by]])
  }
  if(!is.null(split.by)) {
    if(split.by == "donor_type") {
      meta[[split.by]] <- factor(meta[[split.by]], levels = c("Donor", "Recipient"))
    } else {
      meta[[split.by]] <- factor(meta[[split.by]])
    }
  }

  # Add idents filtering after initial factor setup
  if(!is.null(idents)) {
    if(is.null(group.by)) {
      stop("idents parameter requires group.by to be specified")
    }
    if(!all(idents %in% unique(meta[[group.by]]))) {
      missing_idents <- setdiff(idents, unique(meta[[group.by]]))
      stop(sprintf("Some specified identities not found in %s: %s",
                   group.by, paste(missing_idents, collapse=", ")))
    }
    meta <- meta[meta[[group.by]] %in% idents,]
    meta[[group.by]] <- factor(meta[[group.by]], levels = idents)
  }

  # Get gene coordinates and structure
  ah <- AnnotationHub()
  edb_query <- query(ah, c("EnsDb", "Homo sapiens", "104"))
  edb <- edb_query[[1]]
  gene_info <- genes(edb, filter = SymbolFilter(gene))
  exon_info <- exons(edb, filter = SymbolFilter(gene))

  # Validate gene info
  if(length(gene_info) == 0) {
    stop(sprintf("Gene '%s' not found in annotation database", gene))
  }

  # Create grouping structure
  if(is.null(group.by)) {
    meta$group <- "all"
    group.by <- "group"
  }


  # Define plotting region
  plot_start <- start(gene_info) - params$flank_size
  plot_end <- end(gene_info) + params$flank_size
  gene_chr <- as.character(seqnames(gene_info))

  # Get relevant SNPs for the gene
  gene_snps <- which(self$snp_database$snp_metrics$gene_name == gene)
  if(length(gene_snps) == 0) {
    stop(sprintf("No SNPs found for gene: %s", gene))
  }

  # Get unique groups
  unique_groups <- levels(factor(meta[[group.by]]))

  # Create data frame for plotting with progress reporting
  cat(sprintf("\nProcessing SNP data for gene %s...", gene))
  plot_data <- list()
  total_groups <- length(unique_groups)

  for(group_idx in seq_along(unique_groups)) {
    group <- unique_groups[group_idx]
    cat(sprintf("\nProcessing group %d/%d: %s", group_idx, total_groups, group))

    # Get base group cells
    base_group_cells <- which(!is.na(meta[[group.by]]) & meta[[group.by]] == group)

    if(!is.null(split.by)) {
      # Process each split within the group
      split_values <- if(split.by == "donor_type") {
        c("Donor", "Recipient")
      } else {
        levels(factor(meta[[split.by]]))
      }

      for(split_val in split_values) {
        # Get cells for this group-split combination
        split_cells <- base_group_cells[
          !is.na(meta[[split.by]][base_group_cells]) &
            meta[[split.by]][base_group_cells] == split_val
        ]

        if(length(split_cells) >= params$min_cells) {
          # Get SNP data for these cells
          ad_data <- self$snp_database$ad_matrix[gene_snps, split_cells, drop=FALSE]
          dp_data <- if(use_normalized) {
            self$snp_database$dp_matrix_normalized[gene_snps, split_cells, drop=FALSE]
          } else {
            self$snp_database$dp_matrix[gene_snps, split_cells, drop=FALSE]
          }


          # Calculate metrics
          depth_sums <- Matrix::rowSums(dp_data)
          ad_sums <- Matrix::rowSums(ad_data)
          depth_mask <- depth_sums >= params$min_depth

          if(any(depth_mask)) {
            # Apply depth filter first
            valid_snps <- gene_snps[depth_mask]
            depth_filtered <- depth_sums[depth_mask]
            alt_fractions <- ad_sums[depth_mask] / depth_filtered

            # Apply alt fraction filter if specified
            if(!is.null(params$min_alt_frac)) {
              alt_mask <- alt_fractions >= params$min_alt_frac
              if(any(alt_mask)) {
                valid_snps <- valid_snps[alt_mask]
                alt_fractions <- alt_fractions[alt_mask]
                depth_filtered <- depth_filtered[alt_mask]

                if(length(valid_snps) > 0) {
                  # Calculate cells with data for filtered SNPs
                  n_cells <- Matrix::rowSums(dp_data[depth_mask, , drop=FALSE] > 0)[alt_mask]

                  # Debug output
                  cat(sprintf("\nCreating data frame with: %d valid SNPs", length(valid_snps)))
                  cat(sprintf("\n - alt_fractions length: %d", length(alt_fractions)))
                  cat(sprintf("\n - depth_filtered length: %d", length(depth_filtered)))

                  plot_data[[paste(group, split_val, sep="_")]] <- data.frame(
                    group = rep(group, length(valid_snps)),
                    split = rep(split_val, length(valid_snps)),
                    snp_idx = valid_snps,
                    chromosome = self$snp_database$snp_info$CHROM[valid_snps],  # Add this line
                    position = self$snp_database$snp_info$POS[valid_snps],
                    ref = self$snp_database$snp_info$REF[valid_snps],
                    alt = self$snp_database$snp_info$ALT[valid_snps],
                    alt_fraction = alt_fractions,
                    depth = depth_filtered,
                    n_cells = n_cells,
                    feature_type = self$snp_database$snp_annotations$feature_type[valid_snps],
                    stringsAsFactors = FALSE
                  )
                }
              }
            } else {
              # If no alt fraction filter, use depth filtered data directly
              if(length(valid_snps) > 0) {
                n_cells <- Matrix::rowSums(dp_data[depth_mask, , drop=FALSE] > 0)

                plot_data[[paste(group, split_val, sep="_")]] <- data.frame(
                  group = rep(group, length(valid_snps)),
                  split = rep(split_val, length(valid_snps)),
                  snp_idx = valid_snps,
                  position = self$snp_database$snp_info$POS[valid_snps],
                  ref = self$snp_database$snp_info$REF[valid_snps],
                  alt = self$snp_database$snp_info$ALT[valid_snps],
                  alt_fraction = alt_fractions,
                  depth = depth_filtered,
                  n_cells = n_cells,
                  feature_type = self$snp_database$snp_annotations$feature_type[valid_snps],
                  stringsAsFactors = FALSE
                )
              }
            }
          }
        }
      }
    } else {
      # Process group without splitting
      if(length(base_group_cells) >= params$min_cells) {
        ad_data <- self$snp_database$ad_matrix[gene_snps, base_group_cells, drop=FALSE]
        dp_data <- if(use_normalized) {
          self$snp_database$dp_matrix_normalized[gene_snps, split_cells, drop=FALSE]
        } else {
          self$snp_database$dp_matrix[gene_snps, split_cells, drop=FALSE]
        }


        depth_sums <- Matrix::rowSums(dp_data)
        ad_sums <- Matrix::rowSums(ad_data)
        depth_mask <- depth_sums >= params$min_depth

        if(any(depth_mask)) {
          # Apply depth filter first
          valid_snps <- gene_snps[depth_mask]
          depth_filtered <- depth_sums[depth_mask]
          alt_fractions <- ad_sums[depth_mask] / depth_filtered

          # Apply alt fraction filter if specified
          if(!is.null(params$min_alt_frac)) {
            alt_mask <- alt_fractions >= params$min_alt_frac
            if(any(alt_mask)) {
              valid_snps <- valid_snps[alt_mask]
              alt_fractions <- alt_fractions[alt_mask]
              depth_filtered <- depth_filtered[alt_mask]

              if(length(valid_snps) > 0) {
                n_cells <- Matrix::rowSums(dp_data[depth_mask, , drop=FALSE] > 0)[alt_mask]

                plot_data[[as.character(group)]] <- data.frame(
                  group = group,
                  split = NA,
                  snp_idx = valid_snps,
                  position = self$snp_database$snp_info$POS[valid_snps],
                  ref = self$snp_database$snp_info$REF[valid_snps],
                  alt = self$snp_database$snp_info$ALT[valid_snps],
                  alt_fraction = alt_fractions,
                  depth = depth_filtered,
                  n_cells = n_cells,
                  feature_type = self$snp_database$snp_annotations$feature_type[valid_snps],
                  stringsAsFactors = FALSE
                )
              }
            }
          } else {
            # If no alt fraction filter, use depth filtered data directly
            if(length(valid_snps) > 0) {
              n_cells <- Matrix::rowSums(dp_data[depth_mask, , drop=FALSE] > 0)

              plot_data[[as.character(group)]] <- data.frame(
                group = group,
                split = NA,
                snp_idx = valid_snps,
                position = self$snp_database$snp_info$POS[valid_snps],
                ref = self$snp_database$snp_info$REF[valid_snps],
                alt = self$snp_database$snp_info$ALT[valid_snps],
                alt_fraction = alt_fractions,
                depth = depth_filtered,
                n_cells = n_cells,
                feature_type = self$snp_database$snp_annotations$feature_type[valid_snps],
                stringsAsFactors = FALSE
              )
            }
          }
        }
      }
    }
  }


  if(length(plot_data) == 0) {
    stop("No data meeting criteria found for plotting")
  }

  # Combine all data
  plot_df <- do.call(rbind, plot_data)

  # Print summary statistics
  cat("\n\nPlotting Summary:")
  cat(sprintf("\nTotal SNPs: %d", length(unique(plot_df$snp_idx))))
  cat(sprintf("\nFeature types:"))
  print(table(plot_df$feature_type))

  # Create plot grouping structure
  if(!is.null(split.by)) {
    unique_splits <- if(split.by == "donor_type") {
      c("Donor", "Recipient")
    } else {
      levels(factor(plot_df$split))
    }
    plot_df$plot_group <- paste(plot_df$group, plot_df$split, sep="_")
    group_levels <- unlist(lapply(unique_groups, function(g) {
      paste(g, unique_splits, sep="_")
    }))
    plot_df$plot_group <- factor(plot_df$plot_group, levels = group_levels)
  } else {
    plot_df$plot_group <- plot_df$group
  }

  # Create gene track data
  gene_track_data <- data.frame(
    start = start(gene_info),
    end = end(gene_info),
    type = "gene"
  )

  exon_track_data <- data.frame(
    start = start(exon_info),
    end = end(exon_info),
    type = "exon"
  )

  # Calculate optimal range
  range_info <- self$calculate_optimal_range(
    positions = plot_df$position,
    gene_start = gene_track_data$start,
    gene_end = gene_track_data$end,
    min_padding = params$flank_size
  )

  # Calculate group spacing
  y_positions <- self$calculate_y_positions(plot_df, split.by)

  if(data_out) {
    # Create output data frame with specified columns
    output_df <- data.frame(
      group = plot_df$group,
      split = plot_df$split,
      snp_idx = plot_df$snp_idx,
      chromosome = plot_df$chromosome,
      position = plot_df$position,
      ref = plot_df$ref,
      alt = plot_df$alt,
      feature_type = plot_df$feature_type,
      gene_name = rep(gene, nrow(plot_df)),
      gene_type = self$snp_database$snp_annotations$gene_type[plot_df$snp_idx],
      depth = plot_df$depth,
      alt_fraction = plot_df$alt_fraction,
      n_cells = plot_df$n_cells,
      stringsAsFactors = FALSE
    )
    return(output_df)
  }

  # Create main plot
  p_main <- self$create_main_plot(
    plot_df, gene_track_data, exon_track_data, range_info,
    y_positions, gene, gene_chr, color_scheme, point_size_range,
    flank_size = flank_size
  )
  # Create distribution plots

  if(plot_density) {
    p_dist <- self$create_distribution_plots(
      plot_df, split.by, group.by)

    return(plot_grid(
      p_main,
      plot_grid(p_dist$alt_frac, p_dist$depth, ncol = 2),
      ncol = 1,
      rel_heights = c(2, 1)
    ))
  }
  else {
    # Return main plot only
    return(p_main)
  }
})
variantCell$set("public",  "calculate_optimal_range", function(positions, gene_start, gene_end, min_padding = 1000) {
  # Get actual range of SNP positions
  snp_regions <- range(positions)

  # Calculate reasonable padding - use 10% of gene length or min_padding, whichever is larger
  gene_length <- gene_end - gene_start
  padding <- max(min_padding, gene_length * 0.1)

  # Set boundaries to include gene and nearby SNPs with minimal extra space
  plot_start <- min(gene_start, snp_regions[1]) - padding
  plot_end <- max(gene_end, snp_regions[2]) + padding

  return(list(start = plot_start, end = plot_end))
})


variantCell$set("public",  "calculate_y_positions", function(plot_df, split.by) {
  if(!is.null(split.by)) {
    # Create position mapping for groups and their splits
    y_positions <- list()
    current_pos <- 0
    group_spacing <- 0.5
    within_group_spacing <- 0.2

    for(group in unique(plot_df$group)) {
      splits <- unique(plot_df$split[plot_df$group == group])
      y_positions[[group]] <- list()
      for(split in splits) {
        y_positions[[group]][[split]] <- current_pos
        current_pos <- current_pos + within_group_spacing
      }
      current_pos <- current_pos + group_spacing
    }

    # Create label data frames
    group_labels <- data.frame(
      group = names(y_positions),
      y_position = sapply(y_positions, function(x) mean(unlist(x))),
      label_type = "group"
    )

    split_labels <- do.call(rbind, lapply(names(y_positions), function(group) {
      splits <- names(y_positions[[group]])
      data.frame(
        group = group,
        split = splits,
        y_position = unlist(y_positions[[group]]),
        label_type = "split"
      )
    }))

    all_labels <- rbind(
      transform(group_labels, split = NA),
      split_labels
    )

    return(list(
      positions = y_positions,
      labels = all_labels
    ))
  } else {
    positions <- as.numeric(factor(plot_df$group)) - 1
    labels <- data.frame(
      group = unique(plot_df$group),
      split = NA,
      y_position = sort(unique(positions)),
      label_type = "group"
    )
    return(list(
      positions = positions,
      labels = labels
    ))
  }
})

variantCell$set("public",  "create_main_plot", function(plot_df, gene_track_data, exon_track_data, range_info,
                                                        y_positions, gene, gene_chr, color_scheme, point_size_range, flank_size) {
  if(is.list(y_positions$positions)) {
    plot_df$y_position <- sapply(1:nrow(plot_df), function(i) {
      y_positions$positions[[plot_df$group[i]]][[plot_df$split[i]]]
    })
  } else {
    plot_df$y_position <- y_positions$positions
  }

  # Increase left padding for labels
  x_padding_left <- (range_info$end - range_info$start) * 0.15
  x_padding_right <- (range_info$end - range_info$start) * 0.02
  x_min <- range_info$start - x_padding_left
  x_max <- range_info$end + x_padding_right

  # Apply flank size to gene region
  plot_start <- range_info$start - flank_size
  plot_end <- range_info$end + flank_size

  gene_track_y <- -1
  exon_height <- 0.4
  gene_height <- 0.2
  label_offset <- 0.8

  p <- ggplot() +
    geom_segment(data = gene_track_data,
                 aes(x = start, xend = end,
                     y = gene_track_y, yend = gene_track_y),
                 linewidth = gene_height,
                 color = "grey50") +
    annotate("text",
             x = mean(c(gene_track_data$start, gene_track_data$end)),
             y = gene_track_y - label_offset,
             label = gene,
             fontface = "bold") +
    geom_rect(data = exon_track_data,
              aes(xmin = start, xmax = end,
                  ymin = gene_track_y - exon_height/2,
                  ymax = gene_track_y + exon_height/2),
              fill = "black") +
    geom_point(data = plot_df,
               aes(x = position,
                   y = y_position,
                   size = depth,
                   color = alt_fraction)) +
    # Adjust label positions
    geom_text(data = subset(y_positions$labels, label_type == "group"),
              aes(x = x_min + x_padding_left * 0.1,  # Adjusted position
                  y = y_position,
                  label = group),
              hjust = 1,
              fontface = "bold",
              size = 4) +
    geom_text(data = subset(y_positions$labels, label_type == "split" & !is.na(split)),
              aes(x = x_min + x_padding_left * 0.4,  # Adjusted position
                  y = y_position,
                  label = split),
              hjust = 1,
              size = 3.5) +
    scale_color_gradient(low = color_scheme["low"],
                         high = color_scheme["high"],
                         limits = c(0, 1)) +
    scale_size_continuous(range = point_size_range) +
    scale_x_continuous(expand = c(0.02, 0),  # Add small expansion
                       limits = c(x_min, x_max)) +
    coord_cartesian(ylim = c(gene_track_y - label_offset - 0.2,
                             max(plot_df$y_position) + 0.5)) +
    labs(title = sprintf("SNP Distribution in %s", gene),
         subtitle = sprintf("Chr%s:%.0f-%.0f", gene_chr, plot_start, plot_end),
         x = "Genomic Position",
         y = NULL,
         color = "Alt Allele\nFraction",
         size = "Read\nDepth") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.x = element_line(color = "grey90"),
          panel.grid.major.y = element_blank())

  return(p)
})

variantCell$set("public",  "create_distribution_plots", function(plot_df, split.by, group.by) {
  fill_var <- if(!is.null(split.by)) "split" else "group"
  facet_var <- if(!is.null(split.by)) "group" else NULL

  # Alt fraction distribution
  p_alt <- ggplot(plot_df, aes_string(x = "alt_fraction", fill = fill_var)) +
    geom_density(alpha = 0.5) +
    labs(x = "Alternative Allele Fraction",
         y = "Density",
         fill = if(!is.null(split.by)) split.by else group.by) +
    theme_minimal()

  # Depth distribution
  p_depth <- ggplot(plot_df, aes_string(x = "log10(depth)", fill = fill_var)) +
    geom_density(alpha = 0.5) +
    labs(x = "Log10(Read Depth)",
         y = "Density",
         fill = if(!is.null(split.by)) split.by else group.by) +
    theme_minimal()

  if(!is.null(facet_var)) {
    p_alt <- p_alt + facet_wrap(as.formula(paste("~", facet_var)))
    p_depth <- p_depth + facet_wrap(as.formula(paste("~", facet_var)))
  }

  return(list(
    alt_frac = p_alt,
    depth = p_depth
  ))
})
#' @title plotSNPHeatmap: Create a heatmap visualization of SNP expression across cell groups
#' @name plotSNPHeatmap
#'
#' @description
#' This function generates a heatmap visualization of SNP expression data across different cell
#' groups, with optional splitting by an additional variable. It aggregates SNP read counts within
#' each group and can filter based on alternative allele frequency thresholds. The resulting
#' heatmap shows expression patterns with gene-based row annotation and clustering options.
#'
#' @param genes Character vector of gene names to include in the heatmap. Can be NULL if
#'   snp_indices is provided instead.
#' @param snp_indices Integer vector of specific SNP indices to include in the heatmap. Can be NULL
#'   if genes is provided instead.
#' @param group.by Character. Column name in metadata to use for primary grouping of cells.
#' @param split.by Character, optional. Column name in metadata to use for secondary grouping/splitting.
#' @param min_alt_frac Numeric between 0 and 1. Minimum alternative allele fraction required for a SNP to be
#'   counted as expressed in a cell. Default is 0.2.
#' @param scale_data Logical. Whether to scale data by row for visualization. Default is TRUE.
#' @param max_scale Numeric. Maximum value for scaled data (values will be capped at ±max_scale).
#'   Default is 2.
#' @param cluster_rows Logical. Whether to cluster rows (SNPs) in the heatmap. Default is TRUE.
#' @param cluster_cols Logical. Whether to cluster columns (cell groups) in the heatmap. Default is TRUE.
#' @param show_rownames Logical. Whether to display row names (SNP identifiers) in the heatmap.
#'   Default is TRUE.
#' @param show_colnames Logical. Whether to display column names (group identifiers) in the heatmap.
#'   Default is TRUE.
#' @param fontsize_row Numeric. Font size for row names. Default is 8.
#' @param fontsize_col Numeric. Font size for column names. Default is 8.
#' @param exclude_empty Logical. Whether to exclude rows and columns with no expression data.
#'   Default is TRUE.
#' @param normalize_by_cells Logical. Whether to normalize expression values by total cell count in
#'   each group (TRUE) or only by expressing cells (FALSE). Default is TRUE.
#' @param data_out Logical. Whether to return the underlying data matrices instead of the heatmap
#'   object. Default is FALSE.
#'
#' @return If data_out is FALSE (default), returns a ComplexHeatmap object that can be directly
#'   plotted or further customized. If data_out is TRUE, returns a list containing:
#'   \itemize{
#'     \item raw_matrix: Matrix of raw expression values
#'     \item scaled_matrix: Matrix of scaled expression values
#'     \item cell_counts: Matrix of total cell counts per group
#'     \item expr_cell_counts: Matrix of expressing cell counts per group
#'     \item snp_info: Data frame with SNP identifiers, gene names, and feature types
#'   }
#'
#' @details
#' This function calculates the mean expression of SNPs across cell groups, with filtering based
#' on minimum alternative allele frequency. For each SNP in each group, it computes:
#' 1. The number of cells with the SNP
#' 2. The number of cells expressing the SNP above the alt_frac threshold
#' 3. The mean expression value (normalized by total cells or expressing cells)
#'
#' The resulting heatmap includes annotation for gene names and feature types, with rows grouped
#' by gene. The heatmap uses a blue-to-red color scale for scaled expression values.
#'
#' @examples
#'
#' \dontrun{
#'
#' project$plotSNPHeatmap(genes = "BRCA1", group.by = "cell_type")
#'
#' # Plot multiple genes with custom settings
#' project$plotSNPHeatmap(
#'   genes = c("TP53", "KRAS", "EGFR"),
#'   group.by = "cell_type",
#'   split.by = "patient",
#'   min_alt_frac = 0.1,
#'   cluster_rows = FALSE
#' )
#'
#' # Return the underlying data for custom processing
#' snp_data <- project$plotSNPHeatmap(
#'   genes = "APOE",
#'   group.by = "condition",
#'   data_out = TRUE
#' )
#' }
#'
#' @seealso
#' \code{\link{aggregateByGroup}} for aggregating SNP data by groups
#' \code{\link{findSNPsByGroup}} for identifying differential SNPs between groups
#' \code{\link{plotSNPs}} for visualizing SNP distribution along genomic regions
#'
#' @importFrom ComplexHeatmap Heatmap rowAnnotation
#' @importFrom circlize colorRamp2
#' @importFrom grid unit gpar
variantCell$set("public",  "plotSNPHeatmap", function(
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
    data_out = FALSE
) {
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

  # Get metadata
  meta <- self$snp_database$cell_metadata

  # Function to calculate group stats for a specific SNP and cell group
  calculate_group_stats <- function(snp_idx, cell_mask) {
    # Get SNP data for this group
    dp_vals <- self$snp_database$dp_matrix[snp_idx, cell_mask]
    ad_vals <- self$snp_database$ad_matrix[snp_idx, cell_mask]

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

    # Get expression values
    if(is.null(self$snp_database$dp_matrix_normalized)) {
      expr_vals <- dp_vals[alt_mask]
    } else {
      expr_vals <- self$snp_database$dp_matrix_normalized[snp_idx, cell_mask][alt_mask]
    }

    # Calculate mean expression
    total_expr <- sum(expr_vals)
    if(normalize_by_cells) {
      mean_expr <- total_expr / sum(cell_mask)  # Normalize by total cells in group
    } else {
      mean_expr <- total_expr / n_expr_cells  # Normalize only by expressing cells
    }

    return(list(
      mean_expr = mean_expr,
      n_cells = sum(cell_mask),
      n_expr_cells = n_expr_cells
    ))
  }

  # Create group combinations and calculate stats
  if(!is.null(split.by)) {
    group_combinations <- unique(paste(meta[[group.by]], meta[[split.by]], sep="_"))
  } else {
    group_combinations <- unique(meta[[group.by]])
  }

  # Initialize matrices for expression and cell counts
  group_matrix <- matrix(0, nrow=length(snp_indices), ncol=length(group_combinations))
  cell_counts <- matrix(0, nrow=length(snp_indices), ncol=length(group_combinations))
  expr_cell_counts <- matrix(0, nrow=length(snp_indices), ncol=length(group_combinations))
  colnames(group_matrix) <- group_combinations
  colnames(cell_counts) <- group_combinations
  colnames(expr_cell_counts) <- group_combinations

  # Fill matrices
  for(i in seq_along(group_combinations)) {
    if(!is.null(split.by)) {
      combo <- strsplit(group_combinations[i], "_")[[1]]
      mask <- meta[[group.by]] == combo[1] & meta[[split.by]] == combo[2]
    } else {
      mask <- meta[[group.by]] == group_combinations[i]
    }

    if(sum(mask) > 0) {
      for(j in seq_along(snp_indices)) {
        stats <- calculate_group_stats(snp_indices[j], mask)
        group_matrix[j,i] <- stats$mean_expr
        cell_counts[j,i] <- stats$n_cells
        expr_cell_counts[j,i] <- stats$n_expr_cells
      }
    }
  }

  # Filter out empty rows and columns if requested
  if(exclude_empty) {
    # Find non-empty rows and columns
    non_empty_rows <- rowSums(expr_cell_counts) > 0
    non_empty_cols <- colSums(expr_cell_counts) > 0

    if(sum(non_empty_rows) == 0 || sum(non_empty_cols) == 0) {
      stop("No data remaining after filtering empty rows/columns")
    }

    group_matrix <- group_matrix[non_empty_rows, non_empty_cols, drop=FALSE]
    cell_counts <- cell_counts[non_empty_rows, non_empty_cols, drop=FALSE]
    expr_cell_counts <- expr_cell_counts[non_empty_rows, non_empty_cols, drop=FALSE]
    snp_indices <- snp_indices[non_empty_rows]
  }

  # Create row labels
  row_labels <- paste(
    self$snp_database$snp_info$CHROM[snp_indices],
    self$snp_database$snp_info$POS[snp_indices],
    self$snp_database$snp_info$REF[snp_indices],
    self$snp_database$snp_info$ALT[snp_indices],
    sep="_"
  )
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
    return(list(
      raw_matrix = group_matrix,
      scaled_matrix = scaled_matrix,
      cell_counts = cell_counts,
      expr_cell_counts = expr_cell_counts,
      snp_info = data.frame(
        snp_id = row_labels,
        gene = self$snp_database$snp_annotations$gene_name[snp_indices],
        feature_type = self$snp_database$snp_annotations$feature_type[snp_indices]
      )
    ))
  }

  # Create color mapping
  col_fun <- colorRamp2(
    seq(-max_scale, max_scale, length = 7),
    c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#F4A582", "#D6604D", "#B2182B")
  )

  # Create gene annotation
  gene_names <- self$snp_database$snp_annotations$gene_name[snp_indices]
  gene_types <- self$snp_database$snp_annotations$feature_type[snp_indices]

  # Create row split factor for genes
  row_split <- factor(gene_names)

  # Create the heatmap
  ht <- Heatmap(
    scaled_matrix,
    name = "SNP Expression",
    col = col_fun,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_cols,
    show_row_names = show_rownames,
    show_column_names = show_colnames,
    row_names_gp = gpar(fontsize = fontsize_row),
    column_names_gp = gpar(fontsize = fontsize_col),
    row_split = row_split,
    row_gap = unit(2, "mm"),
    left_annotation = rowAnnotation(
      Gene = gene_names,
      Feature = gene_types,
      "Cells" = anno_barplot(rowSums(expr_cell_counts)),
      show_annotation_name = TRUE,
      annotation_legend_param = list(
        Gene = list(nrow = length(unique(gene_names))),
        Feature = list(nrow = length(unique(gene_types)))
      )
    ),
    heatmap_legend_param = list(
      title = if(scale_data) "Z-score" else "Expression",
      title_position = "leftcenter-rot"
    )
  )

  return(ht)
})

