# Annotate a SNP database against a reference VCF

Adds rs identifiers and population allele frequencies to an existing SNP
database by matching on chromosome, position, reference and alternative
allele. Runs standalone on a deserialised project, so a database built
before annotation was available can be annotated without rebuilding.

## Usage

``` r
annotate_snp_database(
  snp_database,
  VCF_file_path,
  add_rs_ids = TRUE,
  add_population_AF = TRUE,
  overwrite = FALSE
)
```

## Arguments

- snp_database:

  The `snp_database` element of a variantCell project.

- VCF_file_path:

  Path to a reference VCF, e.g. a 1000 Genomes panel. Multi-allelic
  records are expanded, and a leading `chr` on contig names is stripped
  so both naming conventions match.

- add_rs_ids:

  Logical. Add rs identifiers from the VCF ID column. Default TRUE.

- add_population_AF:

  Logical. Add population allele frequencies parsed from the INFO AF
  field. Default TRUE.

- overwrite:

  Logical. Replace existing non-missing annotations rather than filling
  only the gaps. Default FALSE.

## Value

The input `snp_database` with `rs_id` and/or `population_AF` columns
added to its `snp_info` element.

## Details

Matching is exact on all four of CHROM, POS, REF and ALT, so a site is
only annotated when the alleles agree as well as the position – a
position-only match would assign the wrong rs ID wherever two
alternative alleles exist at one locus.

This is the plain function rather than the `annotateSNPDatabase()` R6
method. Prefer it when working with a project restored from RDS: R6
methods are serialised with the object, so a deserialised project
carries whatever methods it was saved with and reinstalling the package
does not upgrade them.

On the 32-sample lung transplant cohort this assigns rs IDs to 99.7% of
738,596 sites and population AF to 100%, in about 11 seconds.

## See also

[`buildSNPDatabase`](https://github.com/cchmc/variantCell/reference/buildSNPDatabase.md),
which can annotate at build time.

## Examples

``` r
if (FALSE) { # \dontrun{
project$snp_database <- annotate_snp_database(
  project$snp_database,
  "genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf")

# Fill gaps only, leaving existing annotations untouched
project$snp_database <- annotate_snp_database(
  project$snp_database, vcf_path, overwrite = FALSE)
} # }
```
