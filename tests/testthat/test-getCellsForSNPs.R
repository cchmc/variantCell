# Tests for getCellsForSNPs(), added 2026-08-13.
#
# Two defects found by running the cell-level vignette for the first time:
#
#  1. The documentation has always advertised "chromosome:position format (e.g.
#     1:12345)" as an accepted identifier. It was never implemented - only the
#     full CHROM_POS_REF_ALT key and rs IDs matched. A caller following the docs
#     got an empty character vector and a warning.
#  2. An identifier that resolves to several rows - which any rs ID at a
#     multi-allelic position does - indexed the matrices without drop = FALSE.
#     The result was a matrix, and every downstream operation then ran over the
#     whole matrix instead of one site. Same family as the drop = FALSE defects
#     previously fixed in findDESNPs.
#
# Deterministic fixtures, no RNG.

# Three sites; positions 100 and 200 are the same position with two alternative
# alleles, which is the multi-row case.
make_snp_project <- function() {
  info <- data.frame(
    CHROM = c("1", "1", "2"),
    POS   = c(100L, 100L, 300L),
    REF   = c("A", "A", "C"),
    ALT   = c("G", "T", "T"),
    stringsAsFactors = FALSE
  )
  info$snp_id <- paste(info$CHROM, info$POS, info$REF, info$ALT, sep = "_")
  info$rs_id  <- c("rs1", "rs1", "rs2")   # one rs ID, two alleles

  ids <- paste0("cell", 1:4)
  #                    c1  c2  c3  c4
  ad <- rbind(c( 8,  0,  0,  0),   # 1:100 A>G
              c( 0,  9,  0,  0),   # 1:100 A>T
              c( 0,  0, 10,  0))   # 2:300
  dp <- rbind(c(10, 10, 10,  0),   # DP is a property of the position, so rows
              c(10, 10, 10,  0),   # 1 and 2 carry the same values
              c( 0,  0, 10, 10))
  dimnames(ad) <- dimnames(dp) <- list(info$snp_id, ids)

  p <- variantCell$new()
  p$snp_database <- list(
    ad_matrix = as(ad, "dgCMatrix"),
    dp_matrix = as(dp, "dgCMatrix"),
    snp_info  = info,
    cell_metadata = data.frame(cell_id = ids, sample_id = "S1",
                               stringsAsFactors = FALSE)
  )
  p$samples <- list(S1 = list())
  p
}

quiet <- function(expr) suppressWarnings(invisible(capture.output(res <- expr)))


test_that("chromosome:position is accepted, as documented", {
  p <- make_snp_project()
  quiet(res <- p$getCellsForSNPs("2:300", min_alt_frac = 0.5, min_dp = 5))

  expect_true(attr(res, "summary")$snp_found)
  expect_equal(res[["2:300"]], "cell3")
})


test_that("a chr prefix is tolerated on either side", {
  p <- make_snp_project()
  quiet(res <- p$getCellsForSNPs("chr2:300", min_alt_frac = 0.5, min_dp = 5))
  expect_equal(res[["chr2:300"]], "cell3")
})


test_that("the full key and rs IDs still work", {
  p <- make_snp_project()
  quiet(a <- p$getCellsForSNPs("2_300_C_T", min_alt_frac = 0.5, min_dp = 5))
  quiet(b <- p$getCellsForSNPs("rs2", min_alt_frac = 0.5, min_dp = 5))
  expect_equal(a[["2_300_C_T"]], "cell3")
  expect_equal(b[["rs2"]], "cell3")
})


test_that("an identifier matching several alleles collapses instead of returning a matrix", {
  p <- make_snp_project()
  # rs1 and 1:100 both cover two alternative alleles at one position. cell1
  # carries 8/10 of the G allele and cell2 carries 9/10 of the T allele, so both
  # qualify once alt reads are pooled. cell3 has depth but no alt reads.
  quiet(res <- p$getCellsForSNPs("rs1", min_alt_frac = 0.5, min_dp = 5))
  s <- attr(res, "summary")

  expect_true(s$snp_found)
  expect_setequal(res[["rs1"]], c("cell1", "cell2"))
  # Three cells have coverage at the position, not six - DP must be counted
  # once, not once per allele.
  expect_equal(s$total_cells_with_data, 3)
  expect_equal(s$cells_meeting_criteria, 2)
  expect_equal(s$mean_dp, 10)

  quiet(res2 <- p$getCellsForSNPs("1:100", min_alt_frac = 0.5, min_dp = 5))
  expect_setequal(res2[["1:100"]], c("cell1", "cell2"))
})


test_that("a reference-only query selects cells with coverage and no alt reads", {
  p <- make_snp_project()
  quiet(res <- p$getCellsForSNPs("1:100", min_alt_frac = 0, max_alt_frac = 0,
                                 min_dp = 5))
  expect_equal(res[["1:100"]], "cell3")
})


test_that("an unknown identifier warns and returns an empty vector", {
  p <- make_snp_project()
  expect_warning(
    invisible(capture.output(res <- p$getCellsForSNPs("9:999", min_dp = 1))),
    "not found"
  )
  expect_length(res[["9:999"]], 0)
  expect_false(attr(res, "summary")$snp_found)
})


test_that("sample restriction is honoured", {
  p <- make_snp_project()
  p$snp_database$cell_metadata$sample_id <- c("S1", "S1", "S2", "S2")
  p$samples <- list(S1 = list(), S2 = list())

  quiet(res <- p$getCellsForSNPs("1:100", min_alt_frac = 0.5, min_dp = 5,
                                 sample_ids = "S1"))
  expect_setequal(res[["1:100"]], c("cell1", "cell2"))

  quiet(res2 <- p$getCellsForSNPs("1:100", min_alt_frac = 0.5, min_dp = 5,
                                  sample_ids = "S2"))
  expect_length(res2[["1:100"]], 0)
})
