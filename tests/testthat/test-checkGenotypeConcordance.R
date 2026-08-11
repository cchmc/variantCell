# Tests for checkGenotypeConcordance(), added 2026-08-11.
#
# Germline presence/absence is a genotype test: informative between different
# genomes, structurally empty within one genome, and confounded when a group
# pools several. The guard has to tell those three apart from the data alone.
#
# Fixtures are fully deterministic -- no RNG -- both so the tests are stable and
# because the package must never disturb the caller's random stream.

# Build a database from explicit per-cell genotypes.
#   genotypes: list of integer vectors (0/1/2), one per "individual"
#   assign:    data.frame(individual, sample_id, group, n_cells)
make_geno_db <- function(genotypes, assign, n_sites = 2000) {
  stopifnot(all(assign$individual %in% names(genotypes)))

  cells <- do.call(rbind, lapply(seq_len(nrow(assign)), function(k) {
    a <- assign[k, ]
    data.frame(individual = a$individual, sample_id = a$sample_id,
               group = a$group, n = a$n_cells, stringsAsFactors = FALSE)
  }))
  cell_ind <- rep(cells$individual, cells$n)
  cell_smp <- rep(cells$sample_id, cells$n)
  cell_grp <- rep(cells$group,      cells$n)
  n_cells  <- length(cell_ind)

  # Every cell carries depth 4 at every site; alt reads follow the genotype, so
  # pseudobulk alt fractions land exactly on 0, 0.5 and 1.
  dp <- matrix(4L, nrow = n_sites, ncol = n_cells)
  ad <- vapply(cell_ind, function(ind) as.integer(c(0L, 2L, 4L)[genotypes[[ind]] + 1L]),
               integer(n_sites))

  ids <- paste0("cell", seq_len(n_cells))
  dp <- as(dp, "dgCMatrix"); ad <- as(ad, "dgCMatrix")
  dimnames(dp) <- dimnames(ad) <- list(paste0("snp", seq_len(n_sites)), ids)

  list(ad_matrix = ad, dp_matrix = dp,
       snp_info = data.frame(CHROM = rep(as.character(1:22), length.out = n_sites),
                             POS = seq_len(n_sites), REF = "A", ALT = "G",
                             stringsAsFactors = FALSE),
       cell_metadata = data.frame(cell_id = ids, sample_id = cell_smp,
                                  group = cell_grp, stringsAsFactors = FALSE))
}

# Two individuals sharing 60% of their genotype calls -- below the 0.80
# different-genome threshold, and in the same range as the ~0.61 observed
# between real donor/recipient pairs.
n_sites <- 2000
A <- rep(c(0L, 1L, 2L), length.out = n_sites)
B <- A; B[1201:n_sites] <- (A[1201:n_sites] + 1L) %% 3L   # concordance 0.60
C <- A; C[801:n_sites]  <- (A[801:n_sites]  + 2L) %% 3L   # concordance 0.40

test_that("two distinct genomes are called different_genomes", {
  db <- make_geno_db(list(A = A, B = B),
                     data.frame(individual = c("A", "B"), sample_id = c("S1", "S1"),
                                group = c("g1", "g2"), n_cells = c(100, 100),
                                stringsAsFactors = FALSE))
  res <- checkGenotypeConcordance(db, which(db$cell_metadata$group == "g1"),
                                  which(db$cell_metadata$group == "g2"),
                                  verbose = FALSE)
  expect_equal(res$verdict, "different_genomes")
  expect_equal(res$concordance, 0.6, tolerance = 1e-6)
})

test_that("one genome on both sides is called same_genome", {
  db <- make_geno_db(list(A = A),
                     data.frame(individual = c("A", "A"), sample_id = c("S1", "S2"),
                                group = c("g1", "g2"), n_cells = c(100, 100),
                                stringsAsFactors = FALSE))
  res <- checkGenotypeConcordance(db, which(db$cell_metadata$group == "g1"),
                                  which(db$cell_metadata$group == "g2"),
                                  verbose = FALSE)
  expect_equal(res$verdict, "same_genome")
  expect_equal(res$concordance, 1, tolerance = 1e-6)
})

test_that("a group pooling several genomes is called heterogeneous_groups", {
  # group1 mixes individuals A and B; this is the failure mode of running
  # presence/absence on a pooled donor_type across patients.
  db <- make_geno_db(list(A = A, B = B, C = C),
                     data.frame(individual = c("A", "B", "C"),
                                sample_id  = c("S1", "S2", "S3"),
                                group      = c("g1", "g1", "g2"),
                                n_cells    = c(100, 100, 100),
                                stringsAsFactors = FALSE))
  res <- checkGenotypeConcordance(db, which(db$cell_metadata$group == "g1"),
                                  which(db$cell_metadata$group == "g2"),
                                  split_by = db$cell_metadata$sample_id,
                                  verbose = FALSE)
  expect_equal(res$verdict, "heterogeneous_groups")
  expect_lt(res$within1, 0.90)
  expect_true(is.na(res$within2))   # group2 holds a single sample
})

test_that("heterogeneity outranks a low between-group concordance", {
  # The between-group number alone would read as "different genomes"; the
  # within-group split is what reveals the contrast is not interpretable.
  db <- make_geno_db(list(A = A, B = B, C = C),
                     data.frame(individual = c("A", "B", "C"),
                                sample_id  = c("S1", "S2", "S3"),
                                group      = c("g1", "g1", "g2"),
                                n_cells    = c(100, 100, 100),
                                stringsAsFactors = FALSE))
  res <- checkGenotypeConcordance(db, which(db$cell_metadata$group == "g1"),
                                  which(db$cell_metadata$group == "g2"),
                                  verbose = FALSE)
  expect_lt(res$concordance, 0.80)
  expect_equal(res$verdict, "heterogeneous_groups")
})

test_that("sex chromosomes are excluded", {
  db <- make_geno_db(list(A = A, B = B),
                     data.frame(individual = c("A", "B"), sample_id = c("S1", "S1"),
                                group = c("g1", "g2"), n_cells = c(100, 100),
                                stringsAsFactors = FALSE))
  db$snp_info$CHROM[1:1000] <- "X"
  res <- checkGenotypeConcordance(db, which(db$cell_metadata$group == "g1"),
                                  which(db$cell_metadata$group == "g2"),
                                  verbose = FALSE)
  expect_lte(res$n_sites, 1000)
  expect_equal(res$verdict, "different_genomes")
})

test_that("input problems are rejected rather than silently mishandled", {
  db <- make_geno_db(list(A = A, B = B),
                     data.frame(individual = c("A", "B"), sample_id = c("S1", "S1"),
                                group = c("g1", "g2"), n_cells = c(100, 100),
                                stringsAsFactors = FALSE))
  g1 <- which(db$cell_metadata$group == "g1")
  g2 <- which(db$cell_metadata$group == "g2")

  expect_error(checkGenotypeConcordance(db, g1, c(g2, g1[1]), verbose = FALSE),
               "overlap")
  expect_error(checkGenotypeConcordance(db, g1, g2, min_cells = 1e6, verbose = FALSE),
               "min_cells")
  expect_error(checkGenotypeConcordance(db, db$cell_metadata$cell_id[g1],
                                        c("nosuchcell"), verbose = FALSE),
               "not found")
})

test_that("distinct genomes are counted from concordance, not from labels", {
  # Two individuals in group1 under three sample labels: A appears twice, so
  # the count must be 2, not 3. Label-counting would get this wrong.
  db <- make_geno_db(list(A = A, B = B, C = C),
                     data.frame(individual = c("A", "A", "B", "C"),
                                sample_id  = c("S1", "S2", "S3", "S4"),
                                group      = c("g1", "g1", "g1", "g2"),
                                n_cells    = c(100, 100, 100, 100),
                                stringsAsFactors = FALSE))
  res <- checkGenotypeConcordance(db, which(db$cell_metadata$group == "g1"),
                                  which(db$cell_metadata$group == "g2"),
                                  verbose = FALSE)
  expect_equal(res$n_genomes1, 2L)
  expect_equal(res$n_genomes2, 1L)
})

test_that("the attainable-evidence ceiling is reported for pooled contrasts", {
  db <- make_geno_db(list(A = A, B = B, C = C),
                     data.frame(individual = c("A", "B", "C"),
                                sample_id  = c("S1", "S2", "S3"),
                                group      = c("g1", "g1", "g2"),
                                n_cells    = c(100, 100, 100),
                                stringsAsFactors = FALSE))
  res <- checkGenotypeConcordance(db, which(db$cell_metadata$group == "g1"),
                                  which(db$cell_metadata$group == "g2"),
                                  verbose = FALSE)
  # 2 vs 1 individuals: best two-sided p for perfect separation is 2/choose(3,2)
  expect_equal(res$min_attainable_p, 2 / choose(3, 2))
  expect_false(res$significance_reachable)
  # The message must present this as an underpowered design, not an invalid one.
  expect_match(res$message, "association design")
  expect_match(res$message, "exploratory")
  expect_false(grepl("structurally empty", res$message))
})

test_that("a clean two-genome contrast reports no ceiling and one genome a side", {
  db <- make_geno_db(list(A = A, B = B),
                     data.frame(individual = c("A", "B"), sample_id = c("S1", "S1"),
                                group = c("g1", "g2"), n_cells = c(100, 100),
                                stringsAsFactors = FALSE))
  res <- checkGenotypeConcordance(db, which(db$cell_metadata$group == "g1"),
                                  which(db$cell_metadata$group == "g2"),
                                  verbose = FALSE)
  expect_equal(res$n_genomes1, 1L)
  expect_equal(res$n_genomes2, 1L)
  expect_equal(res$verdict, "different_genomes")
})

test_that("cell ids and integer indices give the same answer", {
  db <- make_geno_db(list(A = A, B = B),
                     data.frame(individual = c("A", "B"), sample_id = c("S1", "S1"),
                                group = c("g1", "g2"), n_cells = c(100, 100),
                                stringsAsFactors = FALSE))
  g1 <- which(db$cell_metadata$group == "g1")
  g2 <- which(db$cell_metadata$group == "g2")
  by_idx <- checkGenotypeConcordance(db, g1, g2, verbose = FALSE)
  by_id  <- checkGenotypeConcordance(db, db$cell_metadata$cell_id[g1],
                                     db$cell_metadata$cell_id[g2], verbose = FALSE)
  expect_equal(by_idx$concordance, by_id$concordance)
  expect_equal(by_idx$verdict, by_id$verdict)
})
