# prioritizeSNPs: Comprehensive SNP Prioritization

Prioritizes SNPs using multiple methodologies including fold change
consistency with gene expression data, machine learning regression, and
LLM-based clinical relevance assessment. This function integrates
differential SNP analysis results with single-cell gene expression data
to identify the most biologically and clinically relevant variants.

## Usage

``` r
prioritizeSNPs(
  snp_df,
  gex_fc_df,
  method = c("fc_consistency", "ml_regression"),
  ml_features = c("alt_frac_diff", "effect_size", "presence_score", "overall_quality",
    "population_AF", "group1_fold_enrichment", "group2_fold_enrichment"),
  top_n_ml = 50,
  top_n_final = 20,
  llm_prompt_template = NULL,
  verbose = TRUE
)
```

## Arguments

- snp_df:

  Data frame. Output from findDESNPs or findSNPsByGroup with
  differential SNP analysis results. Must contain columns: rs_id,
  gene_name, GEX_avg_log2FC, GEX_avg_FC, GEX_p_val, GEX_p_val_adj

- gex_fc_df:

  Data frame. Gene expression fold change data with columns: gene_name
  and Average (average FC across clusters)

- method:

  Character vector. Prioritization methods to use: "fc_consistency" -
  Score based on GEX fold change consistency "ml_regression" - Machine
  learning regression prioritization "llm_clinical" - LLM-based clinical
  relevance scoring. Opt-in only: it issues billable calls to an
  external LLM API. Default: c("fc_consistency", "ml_regression")

- ml_features:

  Character vector. Features to use for ML regression. Default:
  c("alt_frac_diff", "effect_size", "presence_score", "overall_quality",
  "population_AF", "group1_fold_enrichment", "group2_fold_enrichment")

- top_n_ml:

  Integer. Number of top SNPs to select from ML prioritization for LLM
  assessment. Default: 50

- top_n_final:

  Integer. Final number of top SNPs to return. Default: 20

- llm_prompt_template:

  Character. Custom LLM prompt template for clinical assessment.
  Default: NULL (uses built-in template)

- verbose:

  Logical. Whether to print detailed progress information. Default: TRUE

## Value

A list containing:

- prioritized_snps:

  Data frame with top prioritized SNPs and all scoring methods

- fc_consistency_scores:

  Data frame with fold change consistency scores

- ml_scores:

  Data frame with ML regression scores

- llm_assessments:

  Data frame with LLM clinical relevance assessments

- method_weights:

  Weights used to combine different prioritization methods

- summary:

  Summary statistics of the prioritization process

## Details

This function implements a multi-stage SNP prioritization approach:

1.  **Fold Change Consistency**: Compares SNP-associated gene expression
    changes with overall condition-specific expression patterns to
    identify SNPs whose effects align with known disease biology.

2.  **Machine Learning Regression**: Uses multiple SNP characteristics
    (effect size, presence score, population frequency, etc.) to predict
    SNP importance using ensemble methods.

3.  **LLM Clinical Assessment**: Evaluates top-scoring SNPs for clinical
    relevance using structured queries to language models, considering
    gene function, disease association, and therapeutic implications.

The final prioritization combines scores from all methods using weighted
averaging, with higher weights given to methods showing better
concordance.

## Note

For LLM-based clinical assessment, this function requires:

- The ellmer package: install with
  remotes::install_github('tidyverse/ellmer')

- API keys set as environment variables: ANTHROPIC_API_KEY or
  OPENAI_API_KEY

- If LLM integration is unavailable, the function falls back to
  rule-based assessment

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage with all methods
prioritized <- project$prioritizeSNPs(
  snp_df = de_snp_results,
  gex_fc_df = gene_expression_changes
)

# Use only specific methods
prioritized <- project$prioritizeSNPs(
  snp_df = de_snp_results,
  gex_fc_df = gene_expression_changes,
  method = c("fc_consistency", "ml_regression"),
  top_n_final = 10
)

# Custom ML features
prioritized <- project$prioritizeSNPs(
  snp_df = de_snp_results,
  gex_fc_df = gene_expression_changes,
  ml_features = c("effect_size", "presence_score", "population_AF")
)
} # }
```
