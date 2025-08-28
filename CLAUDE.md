# variantCell Package - Claude Code Documentation

## Package Overview

The variantCell package is an R package for analyzing single-cell genetic variants, particularly focused on donor/recipient cell identification in transplant genomics. The package provides comprehensive tools for processing cellSNP data, integrating with vireo donor assignments, and performing differential SNP analysis.

## Core Architecture

### R6 Class Structure
The package is built around an R6 class (`variantCell`) that manages:
- SNP databases from multiple samples
- Cell metadata integration
- Normalization and analysis workflows
- Population allele frequency comparisons

### Key Components

#### 1. Data Processing (`R/02-process-vireo.R`, `R/03-snp-database.R`)
- **Vireo Integration**: Functions to process donor assignment data from vireo
  - `process_vireo_seurat()`: Integrate with Seurat objects
  - `process_vireo_sce()`: Integrate with SingleCellExperiment objects
  - `process_vireo_dataframe()`: Integrate with data frames
- **Sample Management**: `addSampleData()` function for processing individual samples
- **Database Building**: `buildSNPDatabase()` for unified SNP database creation

#### 2. SNP Analysis (`R/04-DE-SNP.R`)
- **Differential Analysis**: `findDESNPs()` for identifying differentially expressed SNPs
- **Group Analysis**: `findSNPsByGroup()` with GroupOnly and SampleStratified modes
- **Population AF Support**: Integration with population allele frequency data from VCF files
- **Fold Enrichment**: Log2 fold change calculations comparing observed vs population AF

#### 3. Visualization (`R/05-plots.R`)
- **Heatmaps**: `plotSNPHeatmap()` for visualizing SNP patterns
- **SNP Plots**: `plotSNPs()` with population AF comparison support

#### 4. Utility Functions (`R/06-utils.R`)
- **Identity Management**: `setProjectIdentity()` and `getCurrentIdentity()`
- **Cell Filtering**: `subsetVariantCell()` for cell subset operations
- **Downsampling**: `downsampleVariant()` for balancing cell numbers
- **Cell Extraction**: `getCellsForSNPs()` for identifying cells expressing specific SNPs

## Major Bug Fix - Memory Management for Large Samples (07/18/25)

### Problem Identified
Large samples (>10,000 cells) with high SNP density were causing memory failures during alt fraction calculation. The issue was particularly evident in:
- CLAD2: 324,330 SNPs × 13,124 cells → 23.3M non-zero entries
- CLAD4: 292,648 SNPs × 10,979 cells → 17.8M non-zero entries

### Root Cause
The sparse matrix operation `alt_fractions[valid_mask] <- ad_matched[valid_mask] / dp_matched[valid_mask]` was attempting to process millions of entries simultaneously, causing memory exhaustion and silent failures.

### Solution Implemented
Added intelligent downsampling based on matrix density rather than fixed cell thresholds:

1. **New Parameter**: `max_nonzero_entries = 12000000` in `addSampleData()`
2. **Dynamic Calculation**: Automatically calculates optimal cell count to stay under threshold
3. **Preserves SNP Coverage**: Preferentially downsamples cells over SNPs
4. **User Configurable**: Threshold can be adjusted per analysis needs

```r
# Usage examples:
project$addSampleData(...) # Uses default 12M threshold
project$addSampleData(..., max_nonzero_entries = 15000000) # More permissive
project$addSampleData(..., max_nonzero_entries = 8000000)  # More conservative
```

### Implementation Details
The fix is located in `R/03-snp-database.R` around lines 260-305:

```r
# Check matrix density before processing
valid_mask <- dp_matched > 0
n_nonzero <- sum(valid_mask)

if(n_nonzero > max_nonzero_entries) {
  # Calculate target downsampling
  target_cells <- min(n_cells, round(max_nonzero_entries * n_cells / n_nonzero))
  target_cells <- max(target_cells, 1000)  # Minimum safety threshold
  
  # Downsample matrices and update all related objects
  # ... (downsampling logic)
}
```

### Results
- **CLAD2**: 23.3M → ~6,100 cells → ~12M entries ✓
- **CLAD4**: 17.8M → ~7,400 cells → ~12M entries ✓
- **Processing Time**: Reduced from hours to minutes
- **Memory Usage**: Stays within manageable limits
- **Statistical Power**: Maintains thousands of cells for robust analysis

## New Feature - Cell Extraction by SNP Expression (07/21/25)

### Function Added: `getCellsForSNPs()`

A new utility function was added to extract cell IDs based on SNP expression criteria, enabling precise cell annotation in downstream analysis tools like Seurat.

#### Key Features:
1. **Flexible SNP Identification**: Accepts both chromosome:position format (e.g., "1:12345") and rs IDs
2. **Dual Threshold Control**: 
   - `min_alt_frac` and `max_alt_frac` for alternative allele fraction range
   - `min_dp` for minimum read depth requirement
3. **Quality Filtering**: Only includes cells meeting both expression and depth criteria
4. **Sample Restriction**: Optional filtering to specific samples
5. **Comprehensive Output**: Returns cell IDs plus summary statistics

#### Usage Examples:

```r
# Find cells expressing SNPs above 20% alt fraction
expressing_cells <- project$getCellsForSNPs(
  snp_ids = c("rs12345", "2:67890"),
  min_alt_frac = 0.2,
  min_dp = 5
)

# Find cells NOT expressing the alternative allele (reference only)
ref_cells <- project$getCellsForSNPs(
  snp_ids = "rs10314",
  min_alt_frac = 0,
  max_alt_frac = 0,
  min_dp = 5
)

# Use results to annotate Seurat object
seurat_obj$snp_status <- ifelse(
  colnames(seurat_obj) %in% expressing_cells$rs12345,
  "expressing", "not_expressing"
)
```

#### Output Structure:
- **Named list**: Each SNP returns a vector of qualifying cell IDs
- **Summary statistics**: Attached as attributes showing:
  - Total cells with data (any reads for the SNP)
  - Cells meeting criteria (passing both alt_frac and depth filters)
  - Mean alternative allele fraction and depth
  - Success rate percentage

#### Quality Control:
The function distinguishes between:
- **`total_cells_with_data`**: Cells with any reads for the SNP (DP > 0)
- **`cells_meeting_criteria`**: Cells passing both expression and minimum depth thresholds

This ensures reliable, high-quality cell annotations by filtering out low-coverage measurements that might yield unreliable alternative allele fractions.

#### Implementation Location:
The function is located in `R/06-utils.R` and integrates seamlessly with the existing SNP database structure, properly handling the `snp_info` metadata and `cell_metadata$cell_id` columns.

## Major New Feature - Comprehensive SNP Prioritization (07/28/25)

### Function Added: `prioritizeSNPs()`

A powerful new analysis function that prioritizes SNPs using multiple methodologies including fold change consistency with gene expression data, machine learning regression, and LLM-based clinical relevance assessment. This function integrates differential SNP analysis results with single-cell gene expression data to identify the most biologically and clinically relevant variants.

#### Multi-Modal Prioritization Approach

The function implements three complementary prioritization methods that can be used individually or in combination:

##### 1. **Fold Change Consistency Scoring**
- Compares SNP-associated gene expression changes with overall condition-specific expression patterns
- Identifies SNPs whose effects align with known disease biology
- Calculates directional matching and magnitude similarity
- Scores on 0-1 scale where higher values indicate better biological consistency

##### 2. **Advanced Machine Learning Regression**
- **Ensemble ML Architecture**: Combines four different modeling approaches:
  - **Linear Model** (35%): Weighted combination with biologically-informed feature weights
  - **Non-linear Model** (25%): Quadratic terms and synergistic interactions
  - **Principal Component Model** (20%): Variance-weighted composite scoring
  - **Rank-based Model** (20%): Robust rank-based prioritization

- **Advanced Feature Engineering**:
  - Interaction features (effect_size × presence_score)
  - Quality-weighted differences (alt_frac_diff × overall_quality)
  - Rare variant scoring based on population frequency
  - Minor allele frequency considerations

- **Robust Processing**:
  - Intelligent missing value imputation
  - Feature normalization and scaling
  - Controlled randomness for tie-breaking

##### 3. **LLM-Based Clinical Relevance Assessment**
- **Integration with ellmer Package**: Uses tidyverse's `ellmer` for structured LLM communication
- **Multi-Provider Support**: Automatically tries Anthropic Claude, then OpenAI GPT
- **Structured Assessment**: Evaluates clinical significance, therapeutic potential, disease associations
- **Batch Processing**: Efficient API usage with rate limiting and error handling
- **Graceful Fallback**: Rule-based clinical assessment when LLM APIs unavailable

#### Key Features

1. **Flexible Method Selection**: Choose any combination of the three approaches
2. **Intelligent API Integration**: 
   - Requires `ellmer` package: `remotes::install_github('tidyverse/ellmer')`
   - Uses environment variables: `ANTHROPIC_API_KEY` or `OPENAI_API_KEY`
   - Falls back to rule-based assessment if LLM unavailable

3. **Comprehensive Output**:
   - Top prioritized SNPs with combined scores
   - Detailed results from each prioritization method
   - Method weights and summary statistics
   - Clinical assessments with therapeutic implications

4. **Robust Error Handling**: Continues analysis even if individual methods fail

#### Usage Examples

```r
# Install required dependencies
remotes::install_github('tidyverse/ellmer')

# Set API key for LLM assessment (optional)
Sys.setenv(ANTHROPIC_API_KEY = "your_anthropic_api_key")
# OR
Sys.setenv(OPENAI_API_KEY = "your_openai_api_key")

# Basic usage with all methods
prioritized <- prioritizeSNPs(
  snp_df = de_snp_results,           # Output from findDESNPs/findSNPsByGroup
  gex_fc_df = gene_expression_fc     # Gene expression fold changes
)

# Use only specific methods
prioritized <- prioritizeSNPs(
  snp_df = de_results,
  gex_fc_df = gex_changes,
  method = c("fc_consistency", "ml_regression"),
  top_n_final = 10
)

# Custom ML feature selection
prioritized <- prioritizeSNPs(
  snp_df = de_results,
  gex_fc_df = gex_changes,
  ml_features = c("effect_size", "presence_score", "population_AF"),
  top_n_ml = 30,
  top_n_final = 15
)

# Custom LLM prompt for specific disease context
custom_prompt <- paste(
  "You are analyzing genetic variants in the context of chronic lung disease.",
  "Focus on variants with known associations to fibrosis, inflammation,",
  "and tissue remodeling pathways. {variant_info}",
  "Provide structured JSON output as specified."
)

prioritized <- prioritizeSNPs(
  snp_df = de_results,
  gex_fc_df = gex_changes,
  llm_prompt_template = custom_prompt
)
```

#### Output Structure

The function returns a comprehensive list containing:

```r
$prioritized_snps        # Top-ranked SNPs with combined priority scores
$fc_consistency_scores   # Detailed fold change consistency analysis
$ml_scores              # ML regression results with component breakdowns
$llm_assessments        # LLM clinical relevance evaluations
$method_weights         # Weighting scheme used for combination
$summary                # Comprehensive analysis summary statistics
```

#### Clinical Assessment Fields

When LLM assessment is used, each SNP receives detailed clinical annotations:

- **Clinical Relevance Score** (0-1): Overall clinical importance
- **Clinical Category**: High/Moderate/Low clinical significance
- **Therapeutic Potential**: Drug_target/Biomarker/Risk_factor/Unknown
- **Disease Association**: Strong_evidence/Moderate_evidence/Weak_evidence/None
- **Mechanism of Action**: Brief description of variant effects
- **Clinical Evidence**: Strong/Moderate/Limited/None
- **Druggability Score** (0-1): Therapeutic targetability
- **Pathway Involvement**: Key biological pathways affected

#### Integration with Analysis Workflow

The prioritization function integrates seamlessly with the existing variantCell workflow:

1. **Differential Analysis**: Run `findDESNPs()` or `findSNPsByGroup()`
2. **Expression Integration**: Use `analyze_snp_differential_expression()` helper
3. **Prioritization**: Apply `prioritizeSNPs()` for comprehensive ranking
4. **Clinical Interpretation**: Review LLM assessments for top variants
5. **Validation**: Use `getCellsForSNPs()` for cell-level validation

#### Performance Considerations

- **Memory Efficient**: Processes SNPs in batches to manage memory usage
- **API Optimization**: Batched LLM calls with rate limiting
- **Scalable**: Handles datasets with thousands of SNPs
- **Robust**: Continues analysis even with partial method failures

#### Implementation Location

The function and all helper functions are located in `R/04-DE-SNP.R` (lines 1977-3113), making it part of the core SNP analysis module. The implementation includes comprehensive error handling, progress reporting, and fallback mechanisms to ensure reliable operation across different computational environments.

## File Structure

```
variantCell/
├── R/
│   ├── 00-pkg-rqmts.r          # Package requirements
│   ├── 01-init.r               # R6 class initialization
│   ├── 02-process-vireo.R      # Vireo integration functions
│   ├── 03-snp-database.R       # Core database management
│   ├── 04-DE-SNP.R            # SNP analysis functions
│   ├── 05-plots.R             # Visualization functions
│   └── 06-utils.R             # Utility functions
├── vignettes/
│   ├── build_snp_database.Rmd # Database building tutorial
│   ├── donor_recipient.Rmd    # Donor/recipient identification
│   ├── cell_level_DE.Rmd      # Cell-level analysis
│   └── group_level_DE.Rmd     # Group-level analysis
└── CLAUDE.md                   # This documentation
```

## Usage Workflow

1. **Initialize Project**: `project <- variantCell$new()`
2. **Add Samples**: Use `addSampleData()` for each sample with cellSNP and vireo data
3. **Build Database**: `buildSNPDatabase()` to create unified SNP database
4. **Set Identity**: `setProjectIdentity()` to define cell groupings
5. **Analyze**: Use `findDESNPs()` or `findSNPsByGroup()` for SNP analysis
6. **Visualize**: Create plots with `plotSNPs()` or `plotSNPHeatmap()`

## Key Features

- **Multi-format Support**: Works with Seurat, SingleCellExperiment, and data frames
- **Population AF Integration**: Compare observed frequencies to population databases
- **Memory Efficient**: Intelligent downsampling for large datasets
- **Comprehensive Analysis**: From data import to statistical analysis and visualization
- **Transplant Focused**: Specialized tools for donor/recipient analysis

## Dependencies

The package integrates with the broader single-cell ecosystem including Seurat, Matrix (sparse matrices), data.table (fast I/O), and standard R statistical functions.