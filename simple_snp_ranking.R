#' Simple SNP Prioritization Function
#' 
#' Ranks SNPs based on key biological and statistical features
#' Focus on LAD-enriched variants with strong gene expression effects
#'
#' @param snp_data Data frame with SNP analysis results
#' @param focus_group Character. Which group to prioritize ("group1" or "group2"). Default: "group1" 
#' @param top_n Integer. Number of top SNPs to return. Default: 50
#' @param weights Named vector. Feature weights for scoring. See details.
#' @return List containing:
#'   \item{full_results}{Data frame with ALL SNPs ranked (includes rank column)}
#'   \item{top_snps}{Data frame with top N SNPs only}
#'   \item{summary}{List with analysis summary and parameters used}
#' 
#' @details
#' Optimized weights based on SNP differential analysis:
#' - presence: 0.35 (presence pattern in target group - most important)
#' - enrichment: 0.25 (population AF fold enrichment - key for group specificity)  
#' - gex_impact: 0.20 (gene expression magnitude & significance)
#' - effect_size: 0.15 (SNP statistical effect)
#' - quality: 0.05 (overall data quality)
#' - alt_frac: 0.00 (removed - redundant with effect_size)
#' - rarity: 0.00 (removed - group-relevant SNPs tend to be common variants)

simple_snp_ranking <- function(snp_data, 
                              focus_group = "group1",
                              top_n = 50,
                              verbose = FALSE, 
                              weights = c(presence = 0.35,      # Increased - most important
                                        enrichment = 0.25,     # High - population comparison key
                                        gex_impact = 0.20,     # Gene expression effects
                                        effect_size = 0.15,    # Statistical significance
                                        quality = 0.05,        # Data quality check
                                        alt_frac = 0.00,       # Remove - redundant with effect_size
                                        rarity = 0.00)) {      # Remove - data shows common SNPs matter
  
  # Validate inputs
  required_cols <- c("rs_id", "gene_name", "effect_size", "overall_quality", 
                    "alt_frac_diff", "population_AF", "GEX_avg_log2FC", "GEX_p_val_adj", "presence")
  
  missing_cols <- setdiff(required_cols, colnames(snp_data))
  if(length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Ensure weights sum to 1
  weights <- weights / sum(weights)
  
  # Detect available enrichment columns
  available_enrichment_cols <- grep("_fold_enrichment$", colnames(snp_data), value = TRUE)
  
  if(length(available_enrichment_cols) == 0) {
    stop("No fold enrichment columns found. Expected columns like 'group1_fold_enrichment'")
  }
  
  # Extract available group numbers/names from enrichment columns
  available_groups <- gsub("_fold_enrichment$", "", available_enrichment_cols)
  
  if(verbose) {
    cat("Available enrichment columns:", paste(available_enrichment_cols, collapse = ", "), "\n")
  }
  
  # Map focus_group to enrichment columns
  # For named groups like 'LAD', 'No_ACR', map to group1/group2 based on presence patterns
  if(focus_group %in% available_groups) {
    # Direct match (e.g., focus_group = "group1")
    enrichment_col <- paste0(focus_group, "_fold_enrichment")
    enrichment_level_col <- paste0(focus_group, "_enrichment_level")
    presence_search_term <- focus_group
  } else {
    # Named group (e.g., focus_group = "LAD") - need to map to group1/group2
    # Check presence column to see which group this represents
    if("presence" %in% colnames(snp_data)) {
      unique_patterns <- unique(snp_data$presence[!is.na(snp_data$presence)])
      focus_in_patterns <- unique_patterns[grepl(focus_group, unique_patterns, ignore.case = TRUE)]
      
      if(length(focus_in_patterns) > 0) {
        # Focus group found in presence patterns, now map to group1/group2
        if("group1_fold_enrichment" %in% available_enrichment_cols) {
          enrichment_col <- "group1_fold_enrichment"
          enrichment_level_col <- "group1_enrichment_level"
          presence_search_term <- focus_group
          if(verbose) cat(sprintf("Mapped focus group '%s' to group1 enrichment columns\n", focus_group))
        } else {
          stop("Found focus group in presence patterns but no group1_fold_enrichment column available")
        }
      } else {
        stop(sprintf("Focus group '%s' not found in presence patterns. Available patterns: %s", 
                    focus_group, paste(unique_patterns, collapse = ", ")))
      }
    } else {
      stop("No 'presence' column found to map named focus group")
    }
  }
  
  if(!enrichment_col %in% colnames(snp_data)) {
    stop("Missing enrichment column: ", enrichment_col)
  }
  
  cat("Calculating component scores for", nrow(snp_data), "SNPs...\n")
  
  # Component 1: Presence Score (0-1, higher = present in target group)
  # Automatically detect presence patterns for the focus group
  presence_score <- rep(0, nrow(snp_data))
  
  # Get unique presence patterns to understand the data
  unique_patterns <- unique(snp_data$presence[!is.na(snp_data$presence)])
  
  # Look for presence patterns that match the focus group
  # Try multiple pattern matching strategies using the presence search term
  focus_patterns <- unique_patterns[grepl(presence_search_term, unique_patterns, ignore.case = TRUE)]
  
  # Also look for patterns that contain the focus group name in "Present in X" format
  present_in_patterns <- unique_patterns[grepl(paste0("Present in.*", presence_search_term), unique_patterns, ignore.case = TRUE)]
  
  # Look for specific patterns
  specific_patterns <- unique_patterns[grepl(paste0(presence_search_term, ".*specific"), unique_patterns, ignore.case = TRUE)]
  
  # Combine all target patterns
  target_patterns <- unique(c(focus_patterns, present_in_patterns, specific_patterns))
  
  # Remove shared/both patterns from target patterns
  target_patterns <- target_patterns[!grepl("both|shared", target_patterns, ignore.case = TRUE)]
  
  # Score presence patterns
  presence_score[snp_data$presence %in% target_patterns] <- 1.0
  presence_score[grepl("both|shared", snp_data$presence, ignore.case = TRUE)] <- 0.3  # Lower score for shared
  
  if(verbose) {
    cat(sprintf("  Detected %d target presence patterns for %s:\n", length(target_patterns), focus_group))
    if(length(target_patterns) > 0) cat("   ", paste(target_patterns, collapse = ", "), "\n")
  }
  
  # Component 2: Enrichment Score (0-1, higher = more enriched in focus group)
  enrichment_raw <- snp_data[[enrichment_col]]
  # Handle infinite values (complete depletion in other group)
  enrichment_raw[is.infinite(enrichment_raw) & enrichment_raw > 0] <- 20  # Cap positive infinity
  enrichment_raw[is.infinite(enrichment_raw) & enrichment_raw < 0] <- -20 # Cap negative infinity
  enrichment_raw[is.na(enrichment_raw)] <- 0
  
  # Normalize enrichment: positive enrichment gets higher scores
  enrichment_score <- pmax(0, pmin(1, (enrichment_raw + 5) / 10))  # Scale -5 to +5 range to 0-1
  
  # Component 3: Gene Expression Impact (0-1, combines magnitude and significance)
  gex_magnitude <- abs(snp_data$GEX_avg_log2FC)
  gex_significance <- -log10(pmax(1e-10, snp_data$GEX_p_val_adj))  # Convert p-adj to -log10 scale
  
  # Normalize components
  gex_mag_norm <- pmin(1, gex_magnitude / 2)  # Cap at 2-fold change
  gex_sig_norm <- pmin(1, gex_significance / 10)  # Cap at 1e-10 p-value
  
  # Combine magnitude and significance (70% magnitude, 30% significance)
  gex_impact_score <- 0.7 * gex_mag_norm + 0.3 * gex_sig_norm
  
  # Component 4: SNP Effect Size (0-1)
  effect_size_raw <- abs(snp_data$effect_size)
  effect_size_raw[is.infinite(effect_size_raw)] <- 20  # Handle infinite effect sizes
  effect_size_score <- pmin(1, effect_size_raw / 20)  # Normalize to 0-1
  
  # Component 5: Quality Score (already 0-1 typically)
  quality_score <- pmin(1, pmax(0, snp_data$overall_quality))
  
  # Component 6: Alternative Allele Fraction Difference (0-1)
  alt_frac_score <- pmin(1, abs(snp_data$alt_frac_diff))
  
  # Component 7: Population Rarity Bonus (0-1, higher for rarer variants)
  # Common variants (AF > 0.3) get low scores, rare variants (AF < 0.05) get high scores
  rarity_score <- 1 - pmin(1, pmax(0, (snp_data$population_AF - 0.01) / 0.49))
  
  # Calculate final weighted score
  final_score <- (
    weights["presence"] * presence_score +
    weights["enrichment"] * enrichment_score +
    weights["gex_impact"] * gex_impact_score + 
    weights["effect_size"] * effect_size_score +
    weights["quality"] * quality_score +
    weights["alt_frac"] * alt_frac_score +
    weights["rarity"] * rarity_score
  )
  
  # Create results dataframe with component scores
  results <- data.frame(
    rank = NA,  # Will be filled after sorting
    rs_id = snp_data$rs_id,
    gene_name = snp_data$gene_name,
    final_score = round(final_score, 4),
    
    # Component scores (for transparency)
    presence_score = round(presence_score, 3),
    enrichment_score = round(enrichment_score, 3),
    gex_impact_score = round(gex_impact_score, 3), 
    effect_size_score = round(effect_size_score, 3),
    quality_score = round(quality_score, 3),
    alt_frac_score = round(alt_frac_score, 3),
    rarity_score = round(rarity_score, 3),
    
    # Original key values (for reference)
    presence_pattern = snp_data$presence,
    enrichment_raw = round(enrichment_raw, 2),
    gex_log2fc = round(snp_data$GEX_avg_log2FC, 3),
    gex_p_adj = snp_data$GEX_p_val_adj,
    population_af = snp_data$population_AF,
    alt_frac_diff = round(snp_data$alt_frac_diff, 3),
    
    stringsAsFactors = FALSE
  )
  
  # Sort by final score (descending) and add ranks
  results <- results[order(results$final_score, decreasing = TRUE), ]
  results$rank <- 1:nrow(results)
  
  # Return top N results
  top_results <- head(results, top_n)
  
  cat(sprintf("Top %d SNPs identified\n", nrow(top_results)))
  cat(sprintf("Score range: %.3f - %.3f\n", max(top_results$final_score), min(top_results$final_score)))
  cat(sprintf("Top SNP: %s (%s) - Score: %.3f\n", 
              top_results$rs_id[1], top_results$gene_name[1], top_results$final_score[1]))
  
  # Prepare return object with both full results and top N
  return_list <- list(
    # Full results dataframe with all SNPs ranked
    full_results = results,
    
    # Top N SNPs only
    top_snps = top_results,
    
    # Summary information
    summary = list(
      total_snps_analyzed = nrow(snp_data),
      top_n_selected = nrow(top_results),
      weights_used = weights,
      focus_group = focus_group,
      score_range = c(min = min(results$final_score), max = max(results$final_score))
    )
  )
  
  return(return_list)
}

# Example usage function
demonstrate_ranking <- function(snp_data) {
  cat("=== Simple SNP Ranking Demonstration ===\n\n")
  
  # Standard ranking focusing on group1 
  cat("1. Standard ranking (focus on group1 enrichment):\n")
  group1_results <- simple_snp_ranking(snp_data, focus_group = "group1", top_n = 20)
  
  # Show top 5
  cat("\nTop 5 group1-enriched candidates:\n")
  print(group1_results$top_snps[1:5, c("rank", "rs_id", "gene_name", "final_score", 
                                      "presence_score", "presence_pattern", "enrichment_score")])
  
  # Alternative weighting: prioritize gene expression impact
  cat("\n\n2. Gene expression focused ranking:\n") 
  gex_focused_weights <- c(presence = 0.20, enrichment = 0.15, gex_impact = 0.35, 
                          effect_size = 0.15, quality = 0.10, alt_frac = 0.05, rarity = 0.00)
  
  gex_results <- simple_snp_ranking(snp_data, focus_group = "group1", top_n = 10,
                                   weights = gex_focused_weights)
  
  cat("\nTop 5 gene expression impact candidates:\n")
  print(gex_results$top_snps[1:5, c("rank", "rs_id", "gene_name", "final_score",
                                   "presence_score", "gex_impact_score", "gex_log2fc")])
  
  return(list(standard = group1_results, gex_focused = gex_results))
}

snp_results <- simple_snp_ranking(snp_LAD, focus_group = 'LAD')

View(snp_results$top_snps)

setwd("~/Desktop/D-Drive/sc-analysis/integrated/11_7_23/SNP-paper")

write.csv(snp_results$full_results, 'simple_snp_results.csv')

read.csv()