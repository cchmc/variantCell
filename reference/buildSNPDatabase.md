# Build a unified SNP database from all samples

Integrates SNP data, annotations, and cell metadata from all samples in
the variantCell project into a unified database. This function creates
merged sparse matrices for alternative allele (AD) and depth (DP) counts
across all cells, combines cell metadata, annotates SNPs with genomic
features, and calculates database-wide metrics. Optionally adds rsID
numbers from VCF file reference.

## Arguments

- add_rs_ids:

  Logical. Whether to add rs# identifiers from a reference VCF file.
  Default: FALSE.

- VCF_file_path:

  Character. Path to reference VCF file (e.g., 1000 Genomes) for rs#
  annotation and/or population AF values. Required if add_rs_ids = TRUE
  or add_population_AF = TRUE.

- add_population_AF:

  Logical. Whether to add population allele frequencies from a reference
  VCF file. Default: FALSE.

## Value

Invisibly returns self (the variantCell object) with the unified SNP
database constructed and stored in the snp_database field.

## Details

This function performs several key steps:

1.  Collects SNP information across all samples and identifies unique
    SNPs

2.  Retrieves genomic annotations for all SNPs (exonic, intronic,
    promoter, etc.)

3.  Combines metadata from all samples, handling missing columns
    appropriately

4.  Creates unified sparse matrices for AD and DP counts across all
    cells

5.  If available, also creates a matrix of normalized counts

6.  Calculates database-wide metrics for each SNP

7.  Generates a QC report with summary statistics

8.  Optionally, adds rsIDs and/or population allele frequencies

The function handles the complexities of integrating data from multiple
samples with potentially different sets of SNPs and metadata columns. It
manages matrix indexing, column alignment, and other technical aspects
needed to build a cohesive database.

After running this function, all subsequent analyses (differential
expression, plotting, etc.) will use the unified database rather than
individual sample data.

## Note

- This function must be called after adding all desired samples with
  [`addSampleData()`](https://github.com/cchmc/variantCell/reference/addSampleData.md)

- The function requires at least one sample to be added

- SNP annotation may take significant time for large datasets

- The resulting database can use substantial memory for projects with
  many cells and SNPs

## Examples

``` r
if (FALSE) { # \dontrun{
# Initialize a variantCell project
project <- variantCell$new()

# Add samples
project$addSampleData(...)
project$addSampleData(...)

# Build the unified SNP database
project$buildSNPDatabase()

# Build with rs# IDs and population AF from 1000 Genomes VCF
project$buildSNPDatabase(add_rs_ids = TRUE, add_population_AF = TRUE,
                         VCF_file_path = "path/to/1000genomes.vcf")

# Now the project is ready for analysis
project$setProjectIdentity("cell_type")
results <- project$findDESNPs(...)
} # }
```
