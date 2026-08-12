# Regression tests for findDESNPs(), added 2026-08-11.
#
# Defect 1: with use_normalized = TRUE (the default), the per-SNP loop computed
#   alt_frac <- ad_raw / dp_NORMALIZED
# dp_matrix_normalized holds log1p(DP * 10000 / colSums(DP)) and there is no
# normalized AD matrix, so the result is not an allele fraction -- it reached
# 10.7 on the real 32-sample database. The min_alt_frac gate was therefore
# effectively inoperative and the reported mean_alt_frac columns were wrong.
#
# Defect 2: log2fc summed depth over gated cells only while the Wilcoxon test
# ran on all cells in the group, so effect size and p-value described different
# populations.
#
# Defect 3: a single qualifying SNP collapsed the matrix subsets to vectors and
# the function errored with "incorrect number of dimensions".

make_test_project <- function(n_snps = 200, n_a = 10, n_b = 10) {
  n_cells <- n_a + n_b

  # Background SNPs carry depth but no alt reads, so only the target SNP
  # survives pre-filtering. Their depth puts colSums in the same regime as the
  # real database (~1000), which is what makes dp_norm land around 4.
  dp <- matrix(5L, nrow = n_snps, ncol = n_cells)
  ad <- matrix(0L, nrow = n_snps, ncol = n_cells)

  # Target SNP (row 1):
  #   group A: 4 cells at alt_frac 0.5, 6 cells at alt_frac 0.1 -> 4 pass a 0.2 gate
  #   group B: 4 cells at alt_frac 0.5, 6 cells at alt_frac 0.0 -> 4 pass
  dp[1, seq_len(n_a)]              <- 10L
  ad[1, 1:4]                       <- 5L
  ad[1, 5:n_a]                     <- 1L
  dp[1, (n_a + 1):n_cells]         <- 4L
  ad[1, (n_a + 1):(n_a + 4)]       <- 2L

  dp <- as(dp, "dgCMatrix")
  ad <- as(ad, "dgCMatrix")
  dp_norm <- as(log1p(as.matrix(dp) %*% diag(10000 / Matrix::colSums(dp))),
                "dgCMatrix")

  cells <- paste0("cell", seq_len(n_cells))
  dimnames(dp) <- dimnames(ad) <- dimnames(dp_norm) <-
    list(paste0("snp", seq_len(n_snps)), cells)

  proj <- variantCell$new()
  proj$snp_database <- list(
    dp_matrix = dp, ad_matrix = ad, dp_matrix_normalized = dp_norm,
    snp_info = data.frame(CHROM = "1", POS = seq_len(n_snps),
                          REF = "A", ALT = "G", stringsAsFactors = FALSE),
    snp_annotations = data.frame(feature_type = "exon", gene_name = "GENE",
                                 gene_type = "protein_coding",
                                 stringsAsFactors = FALSE),
    cell_metadata = data.frame(cell_id = cells, donor_type = "Recipient",
                               grp = rep(c("A", "B"), c(n_a, n_b)),
                               stringsAsFactors = FALSE))
  proj$current_project_ident <- "grp"
  list(project = proj, dp_norm = dp_norm, n_a = n_a, n_b = n_b,
       n_cells = n_cells)
}

test_that("alt fractions are computed from raw depth, not normalized depth", {
  tp <- make_test_project()
  res <- suppressWarnings(capture.output(
    out <- tp$project$findDESNPs("A", "B", use_normalized = TRUE,
                                 min_alt_frac = 0.2, min_expr_cells = 3,
                                 use_parallel = FALSE)))
  row <- out$results[out$results$snp_idx == 1, ]

  af_a <- c(rep(0.5, 4), rep(0.1, 6))
  af_b <- c(rep(0.5, 4), rep(0.0, 6))

  # Only cells at or above min_alt_frac count as expressing. The old code
  # counted all 10 in group A because ad_raw/dp_norm exceeded 0.2 everywhere.
  expect_equal(row$expr_cells_1, sum(af_a >= 0.2))
  expect_equal(row$expr_cells_2, sum(af_b >= 0.2))

  # Reported mean alt fractions must be genuine fractions. The old code
  # returned 0.564 here, and values above 1 in general.
  expect_equal(row$mean_alt_frac_1, mean(af_a))
  expect_equal(row$mean_alt_frac_2, mean(af_b))
  expect_true(row$mean_alt_frac_1 <= 1 && row$mean_alt_frac_2 <= 1)
})

test_that("effect size and test statistic use the same cell population", {
  tp <- make_test_project()
  suppressWarnings(capture.output(
    out <- tp$project$findDESNPs("A", "B", use_normalized = TRUE,
                                 min_alt_frac = 0.2, min_expr_cells = 3,
                                 use_parallel = FALSE)))
  row <- out$results[out$results$snp_idx == 1, ]

  # Group means run over every cell in the group, matching wilcox.test(dp1, dp2).
  exp1 <- sum(tp$dp_norm[1, seq_len(tp$n_a)]) / tp$n_a
  exp2 <- sum(tp$dp_norm[1, (tp$n_a + 1):tp$n_cells]) / tp$n_b

  expect_equal(row$avg_expr_group1, exp1)
  expect_equal(row$avg_expr_group2, exp2)
  expect_equal(row$log2fc, log2((exp1 + 1) / (exp2 + 1)))
})

test_that("serial and parallel paths agree", {
  skip_if_not_installed("doParallel")
  skip_if_not_installed("foreach")
  tp <- make_test_project()
  suppressWarnings(capture.output({
    ser <- tp$project$findDESNPs("A", "B", use_parallel = FALSE)
    par <- tp$project$findDESNPs("A", "B", use_parallel = TRUE, n_cores = 2)
  }))
  cols <- c("log2fc", "avg_expr_group1", "avg_expr_group2", "expr_cells_1",
            "expr_cells_2", "mean_alt_frac_1", "mean_alt_frac_2", "pvalue")
  expect_equal(ser$results[, cols], par$results[, cols],
               ignore_attr = TRUE)
})

test_that("a single qualifying SNP does not collapse the matrix subsets", {
  # Regression for "incorrect number of dimensions": the fixture is built so
  # exactly one SNP survives pre-filtering.
  tp <- make_test_project()
  suppressWarnings(capture.output(
    out <- tp$project$findDESNPs("A", "B", use_parallel = FALSE)))
  expect_equal(nrow(out$results), 1L)
})

test_that("use_normalized = FALSE gives the same alt fractions", {
  tp <- make_test_project()
  suppressWarnings(capture.output(
    out <- tp$project$findDESNPs("A", "B", use_normalized = FALSE,
                                 use_parallel = FALSE)))
  row <- out$results[out$results$snp_idx == 1, ]
  expect_equal(row$mean_alt_frac_1, mean(c(rep(0.5, 4), rep(0.1, 6))))
  expect_equal(row$expr_cells_1, 4)
})
