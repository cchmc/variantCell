# Tests for buildCellSNPDatabase(), added 2026-08-12.
#
# The import path for cellsnp-lite runs against a non-germline region list
# (editing sites, chrM, a candidate panel), which have no vireo output.
#
# The properties worth pinning: sites are unioned rather than intersected across
# samples, the union remap puts counts on the right rows, cell IDs are built so
# they can join external metadata, and OTH survives - it is the position-matched
# error floor that makes an editing signal interpretable.
#
# Fixtures are written to tempdir() and are fully deterministic - no RNG.

write_cellsnp <- function(dir, sites, barcodes, AD, DP, OTH = NULL, gzip_mtx = FALSE) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  vcf <- c("##fileformat=VCFv4.2",
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
           sprintf("%s\t%d\t.\t%s\t%s\t.\t.\t.", sites$CHROM, sites$POS,
                   sites$REF, sites$ALT))
  writeLines(vcf, file.path(dir, "cellSNP.base.vcf"))
  writeLines(barcodes, file.path(dir, "cellSNP.samples.tsv"))
  wm <- function(m, base) {
    p <- file.path(dir, base)
    Matrix::writeMM(methods::as(m, "CsparseMatrix"), p)
    if (gzip_mtx) {
      con <- gzfile(paste0(p, ".gz"), "w")
      writeLines(readLines(p, warn = FALSE), con); close(con); unlink(p)
    }
  }
  wm(DP, "cellSNP.tag.DP.mtx")
  wm(AD, "cellSNP.tag.AD.mtx")
  if (!is.null(OTH)) wm(OTH, "cellSNP.tag.OTH.mtx")
  dir
}

# Sample A covers sites 1-3; sample B covers sites 2-4. Union is 4 sites, and
# only site 2 and 3 are shared -- so an intersection bug shrinks to 2 and a
# naive rbind grows to 6.
sitesA <- data.frame(CHROM = c("1","1","1"), POS = c(100L,200L,300L),
                     REF = "A", ALT = "G", stringsAsFactors = FALSE)
sitesB <- data.frame(CHROM = c("1","1","1"), POS = c(200L,300L,400L),
                     REF = "A", ALT = "G", stringsAsFactors = FALSE)

# Distinct values per (site, sample) so a remap error is visible.
adA  <- Matrix::Matrix(c(1,0, 2,2, 3,3), nrow = 3, byrow = TRUE, sparse = TRUE)
dpA  <- Matrix::Matrix(c(10,10, 20,20, 30,30), nrow = 3, byrow = TRUE, sparse = TRUE)
adB  <- Matrix::Matrix(c(5,5, 6,6, 7,0), nrow = 3, byrow = TRUE, sparse = TRUE)
dpB  <- Matrix::Matrix(c(50,50, 60,60, 70,70), nrow = 3, byrow = TRUE, sparse = TRUE)
# Kept explicitly numeric: a 0/1 matrix through Matrix(sparse=TRUE) becomes a
# *pattern* matrix with no x slot, which is a different code path (covered
# separately below) rather than the one real cellsnp-lite output exercises.
othA <- methods::as(matrix(c(2,0, 0,1, 3,0), nrow = 3, byrow = TRUE), "dgCMatrix")
othB <- methods::as(matrix(c(0,1, 4,0, 0,2), nrow = 3, byrow = TRUE), "dgCMatrix")

make_dirs <- function(gzip_mtx = FALSE) {
  root <- file.path(tempdir(), paste0("cellsnp_fix_", as.integer(gzip_mtx)))
  unlink(root, recursive = TRUE)
  a <- write_cellsnp(file.path(root, "A"), sitesA, c("BC1-1","BC2-1"), adA, dpA, othA, gzip_mtx)
  b <- write_cellsnp(file.path(root, "B"), sitesB, c("BC1-1","BC3-1"), adB, dpB, othB, gzip_mtx)
  c(a, b)
}


test_that("sites are unioned, not intersected or concatenated", {
  d <- make_dirs()
  db <- buildCellSNPDatabase(d, c("A","B"), verbose = FALSE)

  expect_equal(nrow(db$snp_info), 4)
  expect_equal(db$snp_info$POS, c(100L, 200L, 300L, 400L))
  expect_equal(ncol(db$dp_matrix), 4)   # 2 cells per sample
})


test_that("the union remap puts counts on the correct rows", {
  d <- make_dirs()
  db <- buildCellSNPDatabase(d, c("A","B"), verbose = FALSE)
  DP <- as.matrix(db$dp_matrix)
  rownames(DP) <- db$snp_info$POS

  # Sample A cells are columns 1:2, sample B columns 3:4
  expect_equal(unname(DP["100", 1:2]), c(10, 10))
  expect_equal(unname(DP["100", 3:4]), c(0, 0))     # absent in B
  expect_equal(unname(DP["400", 1:2]), c(0, 0))     # absent in A
  expect_equal(unname(DP["400", 3:4]), c(70, 70))
  # shared site carries both samples' values
  expect_equal(unname(DP["200", ]), c(20, 20, 50, 50))
  expect_equal(unname(DP["300", ]), c(30, 30, 60, 60))
})


test_that("cell IDs are prefixed and unique across samples", {
  d <- make_dirs()
  db <- buildCellSNPDatabase(d, c("A","B"), verbose = FALSE)
  expect_equal(db$cell_metadata$cell_id,
               c("A_BC1-1","A_BC2-1","B_BC1-1","B_BC3-1"))
  expect_equal(db$cell_metadata$sample_id, c("A","A","B","B"))
  # BC1-1 occurs in both samples; without prefixing this would collide
  expect_false(anyDuplicated(db$cell_metadata$cell_id) > 0)
})


test_that("OTH survives and the error floor is computed", {
  d <- make_dirs()
  db <- buildCellSNPDatabase(d, c("A","B"), verbose = FALSE)
  expect_false(is.null(db$oth_matrix))
  expect_equal(dim(db$oth_matrix), dim(db$dp_matrix))
  expect_equal(sum(db$oth_matrix), sum(othA) + sum(othB))
  # per-base floor is OTH/DP halved (two non-signal bases)
  expect_equal(db$qc_report$per_base_error_floor,
               sum(db$oth_matrix) / sum(db$dp_matrix) / 2)
})


test_that("a pattern-type .mtx is read as counts of 1, not an error", {
  # writeMM emits a "pattern" header whenever every stored value is 1, and
  # readMM returns that as an n*Matrix with no x slot. An OTH matrix from a
  # shallow run really can be all ones, so this path is reachable in practice,
  # and reading @x blindly would error. Pattern matrices have bitten this
  # package before.
  root <- file.path(tempdir(), "cellsnp_pattern")
  unlink(root, recursive = TRUE)
  pat <- Matrix::Matrix(c(1,0, 0,1, 1,0), nrow = 3, byrow = TRUE, sparse = TRUE)
  probe <- tempfile(); Matrix::writeMM(pat, probe)
  expect_match(readLines(probe, n = 1), "pattern")
  expect_false(methods::.hasSlot(methods::as(Matrix::readMM(probe), "TsparseMatrix"), "x"))
  a <- write_cellsnp(file.path(root, "A"), sitesA, c("BC1-1","BC2-1"), adA, dpA, pat)
  db <- buildCellSNPDatabase(a, "A", verbose = FALSE)
  expect_equal(sum(db$oth_matrix), 3)   # three stored entries, each worth 1
})


test_that("gzipped and plain .mtx give identical results", {
  plain <- buildCellSNPDatabase(make_dirs(FALSE), c("A","B"), verbose = FALSE)
  gz    <- buildCellSNPDatabase(make_dirs(TRUE),  c("A","B"), verbose = FALSE)
  expect_equal(as.matrix(plain$dp_matrix), as.matrix(gz$dp_matrix))
  expect_equal(plain$snp_info, gz$snp_info)
})


test_that("metadata joins and a bad prefix is caught", {
  d <- make_dirs()
  meta <- data.frame(cell_id = c("A_BC1-1","A_BC2-1","B_BC1-1","B_BC3-1"),
                     donor_type = c("Donor","Recipient","Donor","Recipient"),
                     stringsAsFactors = FALSE)
  db <- buildCellSNPDatabase(d, c("A","B"), cell_metadata = meta, verbose = FALSE)
  expect_equal(db$cell_metadata$donor_type,
               c("Donor","Recipient","Donor","Recipient"))

  # a prefix that matches nothing must warn rather than silently produce NAs
  expect_warning(
    buildCellSNPDatabase(d, c("A","B"), cell_metadata = meta,
                         barcode_prefix = c("X_","Y_"), verbose = FALSE),
    "no cells matched"
  )
})


test_that("min_total_dp drops low-coverage sites", {
  d <- make_dirs()
  db <- buildCellSNPDatabase(d, c("A","B"), min_total_dp = 100, verbose = FALSE)
  # totals are 20, 140, 180, 140 -> site 100 drops
  expect_equal(db$snp_info$POS, c(200L, 300L, 400L))
  expect_equal(nrow(db$dp_matrix), 3)
})


test_that("require_metadata drops unmatched cells", {
  d <- make_dirs()
  meta <- data.frame(cell_id = c("A_BC1-1","B_BC1-1"),
                     cell_type = c("Mac","Mac"), stringsAsFactors = FALSE)
  db <- suppressWarnings(
    buildCellSNPDatabase(d, c("A","B"), cell_metadata = meta,
                         require_metadata = TRUE, verbose = FALSE))
  expect_equal(nrow(db$cell_metadata), 2)
  expect_equal(ncol(db$dp_matrix), 2)
  expect_true(all(!is.na(db$cell_metadata$cell_type)))
})


test_that("the imported database drives the allele-fraction functions", {
  d <- make_dirs()
  meta <- data.frame(cell_id = c("A_BC1-1","A_BC2-1","B_BC1-1","B_BC3-1"),
                     donor_type = c("Donor","Recipient","Donor","Recipient"),
                     stringsAsFactors = FALSE)
  p <- variantCellFromCellSNP(dirs = d, sample_ids = c("A","B"),
                              cell_metadata = meta, identity = "donor_type",
                              verbose = FALSE)
  idx <- p$computeAlleleFractionIndex(min_dp_site = 10, min_cells = 1, verbose = FALSE)
  expect_equal(nrow(idx), 4)                       # 2 samples x 2 arms
  expect_true(all(idx$index >= 0 & idx$index <= 1))
  expect_equal(idx$index, idx$total_ad / idx$total_dp, tolerance = 1e-12)
})


test_that("mismatched dirs and sample_ids error", {
  d <- make_dirs()
  expect_error(buildCellSNPDatabase(d, "A", verbose = FALSE), "same length")
  expect_error(buildCellSNPDatabase(d, c("A","A"), verbose = FALSE), "unique")
})


test_that("file connections are not leaked", {
  # readMM does not close a connection handed to it. R caps open connections at
  # ~128, so three leaked per sample would hard-fail a large import long before
  # the garbage collector reclaimed them.
  d <- make_dirs()
  before <- nrow(showConnections())
  for (i in 1:5) invisible(buildCellSNPDatabase(d, c("A","B"), verbose = FALSE))
  expect_equal(nrow(showConnections()), before)
})
