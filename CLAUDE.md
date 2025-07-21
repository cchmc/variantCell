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