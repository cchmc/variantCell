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
#' Default weights:
#' - presence: 0.25 (presence pattern in target group)
#' - enrichment: 0.20 (population AF fold enrichment)  
#' - gex_impact: 0.20 (gene expression magnitude & significance)
#' - effect_size: 0.15 (SNP statistical effect)
#' - quality: 0.10 (overall data quality)
#' - alt_frac: 0.05 (allele fraction difference)
#' - rarity: 0.05 (population rarity bonus)

simple_snp_ranking <- function(snp_data, 
                              focus_group = "group1",
                              top_n = 50,
                              weights = c(presence = 0.25,
                                        enrichment = 0.20,
                                        gex_impact = 0.20, 
                                        effect_size = 0.15,
                                        quality = 0.10,
                                        alt_frac = 0.05,
                                        rarity = 0.05)) {
  
  # Validate inputs
  required_cols <- c("rs_id", "gene_name", "effect_size", "overall_quality", 
                    "alt_frac_diff", "population_AF", "GEX_avg_log2FC", "GEX_p_val_adj", "presence")
  
  missing_cols <- setdiff(required_cols, colnames(snp_data))
  if(length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Ensure weights sum to 1
  weights <- weights / sum(weights)
  
  # Set up enrichment columns based on focus group
  if(focus_group == "group1") {
    enrichment_col <- "group1_fold_enrichment"
    enrichment_level_col <- "group1_enrichment_level"
  } else {
    enrichment_col <- "group2_fold_enrichment" 
    enrichment_level_col <- "group2_enrichment_level"
  }
  
  if(!enrichment_col %in% colnames(snp_data)) {
    stop("Missing enrichment column: ", enrichment_col)
  }
  
  cat("Calculating component scores for", nrow(snp_data), "SNPs...\n")
  
  # Component 1: Presence Score (0-1, higher = present in target group)
  # Map focus group to presence patterns
  if(focus_group == "group1") {
    target_patterns <- c("Present in LAD", "LAD_specific", "group1_specific") 
  } else {
    target_patterns <- c("Present in No_ACR", "No_ACR_specific", "group2_specific")
  }
  
  # Score presence patterns
  presence_score <- rep(0, nrow(snp_data))
  presence_score[snp_data$presence %in% target_patterns] <- 1.0
  presence_score[grepl("both|shared", snp_data$presence, ignore.case = TRUE)] <- 0.3  # Lower score for shared
  
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
  
  # Standard ranking focusing on group1 (LAD)
  cat("1. Standard ranking (focus on LAD enrichment):\n")
  lad_results <- simple_snp_ranking(snp_data, focus_group = "group1", top_n = 20)
  
  # Show top 5
  cat("\nTop 5 LAD-enriched candidates:\n")
  print(lad_results$top_snps[1:5, c("rank", "rs_id", "gene_name", "final_score", 
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
  
  return(list(standard = lad_results, gex_focused = gex_results))
}