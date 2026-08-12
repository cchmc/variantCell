# Tests for computeAlleleFractionIndex() and findDEAlleleFraction(),
# added 2026-08-12.
#
# These test AD/DP -- the fraction of reads carrying the alternative base --
# which is the statistic for RNA-level variation (editing rate, ASE,
# heteroplasmy). findSNPsByGroup tests genotype and findDESNPs tests depth;
# neither tests this.
#
# The two properties that matter most and are easiest to get wrong:
#   1. the unit of observation is the SAMPLE, not the cell (pseudoreplication
#      would inflate significance by orders of magnitude), and
#   2. the permutation ceiling is reported honestly rather than buried.
#
# Fixtures are fully deterministic -- no RNG -- so the tests are stable and the
# package never disturbs the caller's random stream.

# Build a database in which each (sample, arm) pseudobulk has a known,
# exactly-representable alt fraction.
#   spec: data.frame(sample_id, group, n_cells, af)  -- af must be k/dp_per_cell
make_af_db <- function(spec, n_sites = 50, dp_per_cell = 10L) {
  cell_smp <- rep(spec$sample_id, spec$n_cells)
  cell_grp <- rep(spec$group,     spec$n_cells)
  cell_af  <- rep(spec$af,        spec$n_cells)
  n_cells  <- length(cell_smp)

  dp <- matrix(dp_per_cell, nrow = n_sites, ncol = n_cells)
  ad <- matrix(as.integer(round(cell_af * dp_per_cell)),
               nrow = n_sites, ncol = n_cells, byrow = TRUE)

  ids <- paste0("cell", seq_len(n_cells))
  dp <- as(dp, "dgCMatrix"); ad <- as(ad, "dgCMatrix")
  dimnames(dp) <- dimnames(ad) <- list(paste0("snp", seq_len(n_sites)), ids)

  list(ad_matrix = ad, dp_matrix = dp,
       snp_info = data.frame(CHROM = rep("1", n_sites), POS = seq_len(n_sites),
                             REF = "A", ALT = "G",
                             snp_id = paste0("snp", seq_len(n_sites)),
                             stringsAsFactors = FALSE),
       cell_metadata = data.frame(cell_id = ids, sample_id = cell_smp,
                                  donor_type = cell_grp,
                                  stringsAsFactors = FALSE))
}

new_project <- function(db, ident = "donor_type") {
  p <- variantCell$new()
  p$snp_database <- db
  p$current_project_ident <- ident
  p
}

# Six samples, each contributing a Donor and a Recipient pseudobulk. Recipient
# is 0.1 higher in every sample -- a perfectly consistent paired difference.
spec_paired <- data.frame(
  sample_id = rep(paste0("S", 1:6), each = 2),
  group     = rep(c("Donor", "Recipient"), times = 6),
  n_cells   = 40,
  af        = c(0.2, 0.3,  0.3, 0.4,  0.4, 0.5,
                0.2, 0.3,  0.5, 0.6,  0.3, 0.4),
  stringsAsFactors = FALSE
)


test_that("index is read-weighted and recovers the exact planted fraction", {
  p <- new_project(make_af_db(spec_paired))
  idx <- p$computeAlleleFractionIndex(verbose = FALSE)

  expect_equal(nrow(idx), 12)
  expect_setequal(unique(idx$group), c("Donor", "Recipient"))

  d <- idx[idx$group == "Donor" & idx$split == "S1", ]
  expect_equal(d$index, 0.2, tolerance = 1e-9)
  r <- idx[idx$group == "Recipient" & idx$split == "S1", ]
  expect_equal(r$index, 0.3, tolerance = 1e-9)

  # read-weighted index equals total_ad / total_dp by construction
  expect_equal(idx$index, idx$total_ad / idx$total_dp, tolerance = 1e-12)
})


test_that("the paired ceiling is reported and matches 2/2^n", {
  p <- new_project(make_af_db(spec_paired))
  idx <- p$computeAlleleFractionIndex(verbose = FALSE)
  expect_equal(attr(idx, "min_attainable_p_paired"), 2 / 2^6)

  # and it is a real limit: a perfectly consistent 6-pair difference cannot
  # produce a p below it
  d <- idx[idx$group == "Donor", ]
  r <- idx[idx$group == "Recipient", ]
  d <- d[order(d$split), ]; r <- r[order(r$split), ]
  pv <- suppressWarnings(wilcox.test(d$index, r$index, paired = TRUE)$p.value)
  expect_gte(pv, 2 / 2^6 - 1e-12)
})


test_that("common_sites forces both arms onto identical positions", {
  # Make one site uncovered in the Donor arm of S1 only.
  db <- make_af_db(spec_paired)
  donor_s1 <- which(db$cell_metadata$sample_id == "S1" &
                    db$cell_metadata$donor_type == "Donor")
  db$dp_matrix[1, donor_s1] <- 0
  db$ad_matrix[1, donor_s1] <- 0
  p <- new_project(db)

  with_common <- p$computeAlleleFractionIndex(common_sites = TRUE, verbose = FALSE)
  without     <- p$computeAlleleFractionIndex(common_sites = FALSE, verbose = FALSE)

  n_with <- with_common$n_sites[with_common$split == "S1"]
  expect_equal(length(unique(n_with)), 1L)   # both arms use the same count

  n_without <- without$n_sites[without$split == "S1"]
  expect_equal(length(unique(n_without)), 2L)  # arms differ when not forced
})


test_that("cells are not the unit of observation", {
  # Same fractions, but one spec has 10x the cells. If the function were
  # pseudoreplicating, the index would be unchanged but any downstream n would
  # scale with cells. Here n_cells scales and the index does not.
  small <- spec_paired
  big <- spec_paired; big$n_cells <- 400

  i_small <- new_project(make_af_db(small))$computeAlleleFractionIndex(verbose = FALSE)
  i_big   <- new_project(make_af_db(big))$computeAlleleFractionIndex(verbose = FALSE)

  expect_equal(i_small$index, i_big$index, tolerance = 1e-12)
  expect_equal(attr(i_small, "min_attainable_p_paired"),
               attr(i_big, "min_attainable_p_paired"))
  expect_true(all(i_big$n_cells == 10 * i_small$n_cells))
})


test_that("findDEAlleleFraction warns when the ceiling cannot clear Bonferroni", {
  # 3 pairs over 50 sites: floor is 2/2^3 = 0.25, Bonferroni is 0.001.
  spec3 <- spec_paired[spec_paired$sample_id %in% c("S1", "S2", "S3"), ]
  p <- new_project(make_af_db(spec3))

  expect_warning(
    res <- p$findDEAlleleFraction("Donor", "Recipient", min_units = 3, verbose = FALSE),
    "cannot clear Bonferroni"
  )
  expect_false(attr(res, "significance_reachable"))
  expect_equal(attr(res, "min_attainable_p"), 2 / 2^3)
})


test_that("findDEAlleleFraction recovers direction and magnitude", {
  p <- new_project(make_af_db(spec_paired))
  res <- suppressWarnings(
    p$findDEAlleleFraction("Donor", "Recipient", min_units = 3, verbose = FALSE))

  expect_equal(nrow(res), 50)
  expect_true(all(res$n_units == 6))
  # Recipient is uniformly 0.1 higher, so af_1 - af_2 is -0.1 everywhere
  expect_equal(unique(round(res$af_diff, 9)), -0.1)
  expect_true(all(res$p_adj <= 1 & res$p_adj >= 0))
})


test_that("a single testable site does not collapse to a vector", {
  db <- make_af_db(spec_paired, n_sites = 1)
  p <- new_project(db)
  res <- suppressWarnings(
    p$findDEAlleleFraction("Donor", "Recipient", min_units = 3, verbose = FALSE))
  expect_equal(nrow(res), 1)
  expect_false(is.na(res$af_diff[1]))
})


test_that("betabinom and wilcoxon agree on direction", {
  p <- new_project(make_af_db(spec_paired))
  w <- suppressWarnings(p$findDEAlleleFraction("Donor", "Recipient",
                                               test = "wilcoxon", verbose = FALSE))
  b <- suppressWarnings(p$findDEAlleleFraction("Donor", "Recipient",
                                               test = "betabinom", verbose = FALSE))
  key_w <- paste(w$CHROM, w$POS); key_b <- paste(b$CHROM, b$POS)
  expect_setequal(key_w, key_b)
  expect_equal(sign(w$af_diff[match(key_b, key_w)]), sign(b$af_diff))
})


test_that("site restriction accepts chr:pos and indices", {
  p <- new_project(make_af_db(spec_paired))
  by_idx <- p$computeAlleleFractionIndex(sites = 1:5, verbose = FALSE)
  by_key <- p$computeAlleleFractionIndex(sites = paste0("1:", 1:5), verbose = FALSE)
  expect_equal(by_idx$index, by_key$index, tolerance = 1e-12)
  expect_true(all(by_idx$n_sites == 5))
})


test_that("the caller's RNG is untouched", {
  p <- new_project(make_af_db(spec_paired))
  set.seed(1); before <- runif(3)
  set.seed(1)
  invisible(p$computeAlleleFractionIndex(verbose = FALSE))
  invisible(suppressWarnings(p$findDEAlleleFraction("Donor", "Recipient", verbose = FALSE)))
  after <- runif(3)
  expect_equal(before, after)
})


test_that("min_cells is applied before common_sites, not after", {
  # An arm that will be discarded for having too few cells must not constrain
  # the site intersection of the arm that survives. Otherwise the survivor is
  # left with a site set determined by a pseudobulk absent from the output.
  spec <- data.frame(
    sample_id = c("S1","S1"), group = c("Donor","Recipient"),
    n_cells = c(40, 2), af = c(0.2, 0.3), stringsAsFactors = FALSE)
  db <- make_af_db(spec, n_sites = 20)
  # the tiny Recipient arm covers only one site
  tiny <- which(db$cell_metadata$donor_type == "Recipient")
  db$dp_matrix[2:20, tiny] <- 0
  db$ad_matrix[2:20, tiny] <- 0
  p <- new_project(db)

  idx <- p$computeAlleleFractionIndex(min_dp_site = 10, min_cells = 20,
                                      common_sites = TRUE, verbose = FALSE)
  expect_equal(nrow(idx), 1)          # only Donor survives min_cells
  expect_equal(idx$group, "Donor")
  expect_equal(idx$n_sites, 20)       # NOT crippled to 1 by the dropped arm
  expect_equal(idx$index, 0.2, tolerance = 1e-9)
})
