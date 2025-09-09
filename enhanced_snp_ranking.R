#' Enhanced SNP Ranking with Adaptive Scaling and Learning
#' 
#' Improved version with data-driven scaling and ensemble learning approaches
#'
#' @param snp_data Data frame with SNP analysis results
#' @param focus_group Character. Which group to prioritize. Default: "group1"  
#' @param top_n Integer. Number of top SNPs to return. Default: 50
#' @param method Character. Scoring method: "adaptive", "pca", "ensemble". Default: "ensemble"
#' @param verbose Logical. Print detailed progress. Default: TRUE
#' @return List with full_results, top_snps, and method_details

enhanced_snp_ranking <- function(snp_data, 
                                 focus_group = "group1",
                                 top_n = 50,
                                 method = "ensemble",
                                 verbose = TRUE) {
  
  if(verbose) cat("=== Enhanced SNP Ranking ===\n")
  
  # Validate inputs
  required_cols <- c("rs_id", "gene_name", "effect_size", "overall_quality", 
                    "alt_frac_diff", "population_AF", "presence")
  missing_cols <- setdiff(required_cols, colnames(snp_data))
  if(length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Prepare features with adaptive scaling
  features <- prepare_adaptive_features(snp_data, focus_group, verbose)
  
  # Apply selected scoring method
  if(method == "adaptive") {
    scores <- adaptive_weighted_scoring(features, verbose)
  } else if(method == "pca") {
    scores <- pca_scoring(features, verbose)  
  } else if(method == "ensemble") {
    scores <- ensemble_scoring(features, verbose)
  } else {
    stop("Unknown method. Choose: 'adaptive', 'pca', or 'ensemble'")
  }
  
  # Create results
  results <- create_enhanced_results(snp_data, features, scores, top_n, method, verbose)
  
  return(results)
}

# Adaptive feature preparation with data-driven scaling
prepare_adaptive_features <- function(snp_data, focus_group, verbose) {
  
  if(verbose) cat("Preparing features with adaptive scaling...\n")
  
  # Setup group-specific columns
  if(focus_group == "group1" || focus_group == "LAD") {
    enrichment_col <- "group1_fold_enrichment"
    target_patterns <- c("Present in LAD", "LAD_specific", "group1_specific")
  } else if(focus_group == "group2" || focus_group == "No_ACR") {
    enrichment_col <- "group2_fold_enrichment" 
    target_patterns <- c("Present in No_ACR", "No_ACR_specific", "group2_specific")
  } else {
    stop("focus_group must be 'group1', 'LAD', 'group2', or 'No_ACR'")
  }
  
  # Feature 1: Presence (categorical -> numeric)
  presence_score <- rep(0, nrow(snp_data))
  presence_score[snp_data$presence %in% target_patterns] <- 1.0
  presence_score[grepl("both|shared", snp_data$presence, ignore.case = TRUE)] <- 0.3
  
  # Feature 2: Enrichment (handle infinites, then rank-based scaling)
  enrichment_raw <- snp_data[[enrichment_col]]
  enrichment_raw[is.infinite(enrichment_raw) & enrichment_raw > 0] <- max(enrichment_raw[is.finite(enrichment_raw)]) * 2
  enrichment_raw[is.infinite(enrichment_raw) & enrichment_raw < 0] <- min(enrichment_raw[is.finite(enrichment_raw)]) * 2
  enrichment_raw[is.na(enrichment_raw)] <- 0
  # Rank-based scaling (robust to outliers)
  enrichment_score <- rank(enrichment_raw) / length(enrichment_raw)
  
  # Feature 3: Effect size (log transformation + quantile scaling)
  effect_size_raw <- abs(snp_data$effect_size)
  effect_size_raw[is.infinite(effect_size_raw)] <- max(effect_size_raw[is.finite(effect_size_raw)])
  effect_size_log <- log1p(effect_size_raw)  # log(1 + x) handles zeros
  effect_size_score <- pmin(1, (effect_size_log - quantile(effect_size_log, 0.05)) / 
                                (quantile(effect_size_log, 0.95) - quantile(effect_size_log, 0.05)))
  effect_size_score <- pmax(0, effect_size_score)
  
  # Feature 4: Gene expression (if available and not all missing)
  if(all(c("GEX_avg_log2FC", "GEX_p_val_adj") %in% colnames(snp_data))) {
    gex_magnitude <- abs(snp_data$GEX_avg_log2FC)
    gex_significance <- -log10(pmax(1e-15, snp_data$GEX_p_val_adj))
    
    # Check if we have valid (non-NA) gene expression data
    valid_gex_mag <- !is.na(gex_magnitude) & is.finite(gex_magnitude)
    valid_gex_sig <- !is.na(gex_significance) & is.finite(gex_significance)
    
    if(sum(valid_gex_mag) > 0 && sum(valid_gex_sig) > 0) {
      # Use valid data for scaling
      gex_mag_scaled <- rep(0.5, length(gex_magnitude))  # Default neutral
      gex_mag_scaled[valid_gex_mag] <- rank(gex_magnitude[valid_gex_mag]) / sum(valid_gex_mag)
      
      gex_sig_scaled <- rep(0.5, length(gex_significance))  # Default neutral
      if(sum(valid_gex_sig) > 1) {
        sig_90th <- quantile(gex_significance[valid_gex_sig], 0.9, na.rm=TRUE)
        gex_sig_scaled[valid_gex_sig] <- pmin(1, gex_significance[valid_gex_sig] / sig_90th)
      }
      
      gex_combined <- 0.7 * gex_mag_scaled + 0.3 * gex_sig_scaled
    } else {
      if(verbose) cat("  Warning: Gene expression data is all NA, using neutral scores\n")
      gex_combined <- rep(0.5, nrow(snp_data))  # Neutral if all missing
    }
  } else {
    if(verbose) cat("  Gene expression columns not found, using neutral scores\n")
    gex_combined <- rep(0.5, nrow(snp_data))  # Neutral if missing
  }
  
  # Feature 5: Quality (already 0-1, but handle missing)
  quality_score <- snp_data$overall_quality
  quality_score[is.na(quality_score)] <- median(quality_score, na.rm=TRUE)
  quality_score <- pmax(0, pmin(1, quality_score))
  
  # Feature 6: Alt fraction difference (sigmoid transformation for diminishing returns)
  alt_frac_raw <- abs(snp_data$alt_frac_diff)
  alt_frac_score <- 2 / (1 + exp(-4 * alt_frac_raw)) - 1  # Sigmoid: 0 to ~1
  
  # Feature 7: Population rarity (inverse frequency with plateau)
  pop_af <- snp_data$population_AF
  pop_af[is.na(pop_af)] <- 0.1  # Assume moderate frequency if missing
  rarity_score <- 1 / (1 + exp(10 * (pop_af - 0.1)))  # Sigmoid: rare variants favored
  
  # Combine into feature matrix
  feature_matrix <- data.frame(
    presence = presence_score,
    enrichment = enrichment_score,
    effect_size = effect_size_score,
    gex_impact = gex_combined,
    quality = quality_score,
    alt_frac = alt_frac_score,
    rarity = rarity_score
  )
  
  if(verbose) {
    cat(sprintf("  Features prepared for %d SNPs\n", nrow(feature_matrix)))
    cat("  Feature correlations:\n")
    print(round(cor(feature_matrix), 2))
  }
  
  return(feature_matrix)
}

# Method 1: Adaptive weighted scoring with correlation adjustment
adaptive_weighted_scoring <- function(features, verbose) {
  
  if(verbose) cat("Computing adaptive weighted scores...\n")
  
  # Calculate feature importance based on variance and orthogonality
  feature_vars <- apply(features, 2, var)
  feature_corr <- cor(features)
  
  # Reduce weights for highly correlated features
  base_weights <- feature_vars / sum(feature_vars)  # Variance-based weights
  
  # Correlation penalty: reduce weight if highly correlated with others
  for(i in 1:length(base_weights)) {
    max_corr <- max(abs(feature_corr[i, -i]))
    base_weights[i] <- base_weights[i] * (1 - 0.5 * max_corr)  # Penalty for correlation
  }
  
  # Normalize weights
  adaptive_weights <- base_weights / sum(base_weights)
  names(adaptive_weights) <- colnames(features)
  
  if(verbose) {
    cat("  Adaptive weights:\n")
    print(round(adaptive_weights, 3))
  }
  
  # Calculate weighted score
  scores <- as.numeric(as.matrix(features) %*% adaptive_weights)
  
  return(list(scores = scores, weights = adaptive_weights, method = "adaptive"))
}

# Method 2: PCA-based scoring
pca_scoring <- function(features, verbose) {
  
  if(verbose) cat("Computing PCA-based scores...\n")
  
  # Check for and handle missing/infinite values
  features_clean <- features
  for(i in 1:ncol(features_clean)) {
    col_data <- features_clean[, i]
    if(any(is.na(col_data) | is.infinite(col_data))) {
      # Replace NA/Inf with column median
      col_median <- median(col_data, na.rm = TRUE)
      if(is.na(col_median)) col_median <- 0.5  # If all NA, use neutral
      features_clean[is.na(col_data) | is.infinite(col_data), i] <- col_median
      if(verbose) cat(sprintf("  Fixed %d missing/infinite values in %s\n", 
                             sum(is.na(col_data) | is.infinite(col_data)), colnames(features)[i]))
    }
  }
  
  # Perform PCA on cleaned data
  pca_result <- prcomp(features_clean, center = TRUE, scale. = TRUE)
  
  # Use first PC as primary score, weight by explained variance
  pc1_scores <- pca_result$x[, 1]
  pc2_scores <- pca_result$x[, 2]
  
  # Combine first two PCs weighted by their explained variance
  var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  combined_scores <- pc1_scores * var_explained[1] + pc2_scores * var_explained[2]
  
  # Scale to 0-1 range
  final_scores <- (combined_scores - min(combined_scores)) / 
                  (max(combined_scores) - min(combined_scores))
  
  if(verbose) {
    cat(sprintf("  PC1 explains %.1f%%, PC2 explains %.1f%% of variance\n", 
                var_explained[1]*100, var_explained[2]*100))
    cat("  Feature loadings on PC1:\n")
    print(round(pca_result$rotation[, 1], 3))
  }
  
  return(list(scores = final_scores, 
              pca_result = pca_result, 
              var_explained = var_explained,
              method = "pca"))
}

# Method 3: Ensemble scoring (combines multiple approaches)
ensemble_scoring <- function(features, verbose) {
  
  if(verbose) cat("Computing ensemble scores...\n")
  
  # Get scores from different methods
  adaptive_result <- adaptive_weighted_scoring(features, verbose = FALSE)
  pca_result <- pca_scoring(features, verbose = FALSE)
  
  # Simple rank-based scoring for comparison
  rank_scores <- apply(features, 1, function(x) {
    ranks <- rank(x, na.last = "keep")
    mean(ranks / length(x), na.rm = TRUE)
  })
  # Handle case where all values are the same
  if(all(is.na(rank_scores)) || max(rank_scores, na.rm=TRUE) == min(rank_scores, na.rm=TRUE)) {
    rank_scores <- rep(0.5, length(rank_scores))
  } else {
    rank_scores[is.na(rank_scores)] <- 0.5
    rank_scores <- (rank_scores - min(rank_scores, na.rm=TRUE)) / (max(rank_scores, na.rm=TRUE) - min(rank_scores, na.rm=TRUE))
  }
  
  # Ensemble weights (can be tuned)
  ensemble_weights <- c(adaptive = 0.5, pca = 0.3, rank_based = 0.2)
  
  # Combine scores
  final_scores <- (
    ensemble_weights["adaptive"] * adaptive_result$scores +
    ensemble_weights["pca"] * pca_result$scores +  
    ensemble_weights["rank_based"] * rank_scores
  )
  
  # Normalize to 0-1
  final_scores <- (final_scores - min(final_scores)) / (max(final_scores) - min(final_scores))
  
  if(verbose) {
    cat("  Ensemble method weights:\n")
    print(ensemble_weights)
    cat("  Score correlation between methods:\n")
    score_matrix <- cbind(adaptive_result$scores, pca_result$scores, rank_scores)
    colnames(score_matrix) <- c("Adaptive", "PCA", "Rank")
    print(round(cor(score_matrix), 3))
  }
  
  return(list(scores = final_scores,
              component_scores = list(
                adaptive = adaptive_result$scores,
                pca = pca_result$scores, 
                rank_based = rank_scores
              ),
              adaptive_weights = adaptive_result$weights,
              pca_result = pca_result$pca_result,
              method = "ensemble"))
}

# Create enhanced results dataframe
create_enhanced_results <- function(snp_data, features, score_result, top_n, method, verbose) {
  
  # Create full results dataframe
  results_df <- data.frame(
    rank = NA,
    rs_id = snp_data$rs_id,
    gene_name = snp_data$gene_name,
    final_score = round(score_result$scores, 4),
    
    # Feature scores for transparency
    presence_score = round(features$presence, 3),
    enrichment_score = round(features$enrichment, 3), 
    effect_size_score = round(features$effect_size, 3),
    gex_impact_score = round(features$gex_impact, 3),
    quality_score = round(features$quality, 3),
    alt_frac_score = round(features$alt_frac, 3),
    rarity_score = round(features$rarity, 3),
    
    # Original values for reference
    presence_pattern = snp_data$presence,
    population_af = snp_data$population_AF,
    
    stringsAsFactors = FALSE
  )
  
  # Add component scores for ensemble method
  if(method == "ensemble") {
    results_df$adaptive_component = round(score_result$component_scores$adaptive, 3)
    results_df$pca_component = round(score_result$component_scores$pca, 3)  
    results_df$rank_component = round(score_result$component_scores$rank_based, 3)
  }
  
  # Sort and rank
  results_df <- results_df[order(results_df$final_score, decreasing = TRUE), ]
  results_df$rank <- 1:nrow(results_df)
  
  # Get top N
  top_results <- head(results_df, top_n)
  
  if(verbose) {
    cat(sprintf("Ranking complete! Top SNP: %s (%.3f)\n", 
                top_results$rs_id[1], top_results$final_score[1]))
  }
  
  # Return comprehensive results
  return(list(
    full_results = results_df,
    top_snps = top_results,
    method_details = score_result,
    feature_matrix = features,
    summary = list(
      method = method,
      total_snps = nrow(snp_data),
      top_n = top_n,
      score_range = range(score_result$scores)
    )
  ))
}

# Comparison function to evaluate different methods
compare_ranking_methods <- function(snp_data, focus_group = "group1", top_n = 20) {
  
  cat("=== Comparing Ranking Methods ===\n\n")
  
  # Run all three methods
  adaptive_results <- enhanced_snp_ranking(snp_data, focus_group, top_n, "adaptive", verbose = FALSE)
  pca_results <- enhanced_snp_ranking(snp_data, focus_group, top_n, "pca", verbose = FALSE)
  ensemble_results <- enhanced_snp_ranking(snp_data, focus_group, top_n, "ensemble", verbose = FALSE)
  
  # Compare top 10 from each method
  cat("Top 10 SNPs by method:\n")
  comparison <- data.frame(
    Rank = 1:10,
    Adaptive = adaptive_results$top_snps$rs_id[1:10],
    PCA = pca_results$top_snps$rs_id[1:10],
    Ensemble = ensemble_results$top_snps$rs_id[1:10]
  )
  print(comparison)
  
  # Calculate overlap
  adaptive_top10 <- adaptive_results$top_snps$rs_id[1:10]
  pca_top10 <- pca_results$top_snps$rs_id[1:10]
  ensemble_top10 <- ensemble_results$top_snps$rs_id[1:10]
  
  cat("\nMethod overlap (top 10):\n")
  cat(sprintf("Adaptive vs PCA: %d/10 shared\n", length(intersect(adaptive_top10, pca_top10))))
  cat(sprintf("Adaptive vs Ensemble: %d/10 shared\n", length(intersect(adaptive_top10, ensemble_top10))))
  cat(sprintf("PCA vs Ensemble: %d/10 shared\n", length(intersect(pca_top10, ensemble_top10))))
  
  return(list(adaptive = adaptive_results, pca = pca_results, ensemble = ensemble_results))
}

snp_data <- read.csv()

View(snp_data)

snp_results_compare <- compare_ranking_methods(snp_data=snp_data, focus_group = 'LAD')

snp_results <- enhanced_snp_ranking(snp_data = snp_data, focus_group = 'LAD')

View(snp_results$top_snps)

write.csv(snp_results$full_results, "")
