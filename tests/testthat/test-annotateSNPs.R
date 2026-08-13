# Tests for annotate_snps(), added 2026-08-13.
#
# Found while making the vignettes executable: the transcript_ids assignment
# indexed the left-hand side with one element per *overlap* while the right-hand
# side had one element per distinct *SNP*, so R recycled the shorter vector and
# attached transcript IDs to the wrong SNPs. It surfaced only as a "number of
# items to replace is not a multiple of replacement length" warning on an
# otherwise passing run - the same way both cellsnp-import bugs surfaced.
#
# These run against EnsDb.Hsapiens.v86 at real coordinates, so they are skipped
# if the annotation package is not installed.

skip_if_no_ensdb <- function() {
  testthat::skip_if_not_installed("EnsDb.Hsapiens.v86")
  testthat::skip_if_not_installed("ensembldb")
}

# Positions taken from real exons of two genes on different chromosomes. Guessing
# evenly-spaced coordinates across a gene body lands most of them in introns,
# where there is no exon or transcript ID to misassign in the first place.
exonic_positions <- function(syms = c("JAK1", "CD74"), per_gene = 8) {
  edb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
  ex <- ensembldb::exons(edb, filter = ~ symbol %in% syms, columns = "symbol")
  ex <- ex[as.character(GenomicRanges::seqnames(ex)) %in% c(1:22, "X")]
  sym <- GenomicRanges::mcols(ex)$symbol
  do.call(rbind, lapply(syms, function(s) {
    e <- ex[sym == s]
    e <- e[order(-GenomicRanges::width(e))][seq_len(min(per_gene, length(e)))]
    data.frame(CHROM = as.character(GenomicRanges::seqnames(e)),
               POS = as.integer(GenomicRanges::start(e) + 5),
               gene = s, REF = "A", ALT = "G", stringsAsFactors = FALSE)
  }))
}


test_that("annotation emits no recycling warning", {
  skip_if_no_ensdb()
  snps <- exonic_positions()
  p <- variantCell$new()
  expect_no_warning(invisible(capture.output(p$annotate_snps(snps))))
})


test_that("transcript IDs land on the SNP they were derived from", {
  skip_if_no_ensdb()
  snps <- exonic_positions()
  p <- variantCell$new()
  invisible(capture.output(ann <- p$annotate_snps(snps)))

  expect_equal(nrow(ann), nrow(snps))
  expect_equal(ann$gene_name, snps$gene)

  # Every stored transcript must actually span the SNP it is attached to. Under
  # the recycling bug, IDs cycled onto SNPs on an entirely different chromosome.
  edb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
  tx <- ensembldb::transcripts(edb)
  has_tx <- which(!is.na(ann$transcript_ids))
  expect_gt(length(has_tx), 0)
  for (i in has_tx) {
    ids <- strsplit(ann$transcript_ids[i], ";")[[1]]
    hit <- tx[stats::na.omit(match(ids, GenomicRanges::mcols(tx)$tx_id))]
    expect_equal(length(hit), length(ids))
    expect_true(
      all(as.character(GenomicRanges::seqnames(hit)) == snps$CHROM[i] &
            GenomicRanges::start(hit) <= snps$POS[i] &
            GenomicRanges::end(hit) >= snps$POS[i]),
      info = sprintf("SNP %d (%s:%d, %s) got a non-overlapping transcript",
                     i, snps$CHROM[i], snps$POS[i], snps$gene[i]))
  }
})


test_that("exon IDs land on the SNP they were derived from", {
  skip_if_no_ensdb()
  snps <- exonic_positions()
  p <- variantCell$new()
  invisible(capture.output(ann <- p$annotate_snps(snps)))

  expect_true(all(ann$in_exon))
  expect_true(all(ann$feature_type == "exonic"))
  expect_true(all(!is.na(ann$exon_ids)))

  edb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
  ex <- ensembldb::exons(edb)
  for (i in seq_len(nrow(ann))) {
    ids <- strsplit(ann$exon_ids[i], ";")[[1]]
    hit <- ex[stats::na.omit(match(ids, GenomicRanges::mcols(ex)$exon_id))]
    expect_equal(length(hit), length(ids))
    expect_true(
      all(as.character(GenomicRanges::seqnames(hit)) == snps$CHROM[i] &
            GenomicRanges::start(hit) <= snps$POS[i] &
            GenomicRanges::end(hit) >= snps$POS[i]),
      info = sprintf("SNP %d (%s:%d, %s) got a non-overlapping exon",
                     i, snps$CHROM[i], snps$POS[i], snps$gene[i]))
  }
})


test_that("overlapping features are collapsed, not overwritten", {
  skip_if_no_ensdb()
  # A position inside a widely-shared exon belongs to many transcripts at once.
  # Keeping only the last one silently discarded most of the annotation.
  snps <- exonic_positions("CD74", per_gene = 4)
  p <- variantCell$new()
  invisible(capture.output(ann <- p$annotate_snps(snps)))

  n_tx <- vapply(strsplit(ann$transcript_ids, ";"), length, integer(1))
  expect_gt(max(n_tx), 1)
})
