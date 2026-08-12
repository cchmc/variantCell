# ---------------------------------------------------------------------------
# Import cellsnp-lite output directly into a variantCell SNP database.
#
# addSampleData() builds a database from a cellSNP run paired with vireo donor
# assignments, which is the germline genotyping workflow. A run against a
# non-germline region list - RNA editing sites, chrM, a candidate panel - has no
# vireo output and no donor assignment, because there is no genotype to
# deconvolve. This module imports those runs.
#
# It also carries the OTH matrix through. At an A>G site the two remaining bases
# are pure sequencing error, so OTH/DP is an internal, position-matched error
# floor for the AD/DP signal. That comparison is what distinguishes a real
# editing signal from a base-calling artifact, and it is lost if OTH is dropped
# at import.
#
# Matrices are assembled once from accumulated triplets rather than by repeated
# cbind, which is quadratic in the number of samples.
# ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
#' Open a cellsnp-lite output file whether or not it was written with --gzip.
.cellsnp_open <- function(dir, base) {
  plain <- file.path(dir, base)
  gz    <- paste0(plain, ".gz")
  if (file.exists(plain)) return(list(path = plain, con = function() file(plain, "r")))
  if (file.exists(gz))    return(list(path = gz,    con = function() gzfile(gz, "r")))
  stop(sprintf("neither %s nor %s found", plain, gz), call. = FALSE)
}

#' @keywords internal
#' @noRd
#' Read one cellsnp-lite output directory. Returns triplets rather than a
#' matrix so the caller can remap row indices before a single assembly.
.read_cellsnp_sample <- function(dir, need_oth = TRUE) {
  if (!dir.exists(dir)) stop(sprintf("directory not found: %s", dir), call. = FALSE)

  vcf <- .cellsnp_open(dir, "cellSNP.base.vcf")
  con <- vcf$con(); on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
  lines <- lines[!startsWith(lines, "#")]
  if (length(lines) == 0) stop(sprintf("no variants in %s", vcf$path), call. = FALSE)
  f <- do.call(rbind, strsplit(lines, "\t", fixed = TRUE))
  sites <- data.frame(CHROM = f[, 1], POS = as.integer(f[, 2]),
                      REF = f[, 4], ALT = f[, 5], stringsAsFactors = FALSE)

  bc <- .cellsnp_open(dir, "cellSNP.samples.tsv")
  bcon <- bc$con(); barcodes <- readLines(bcon, warn = FALSE); close(bcon)
  barcodes <- barcodes[nzchar(barcodes)]

  # readMM does not close the connection it is handed. R caps open connections
  # at ~128, so leaking three per sample would hard-fail a large import.
  read_mtx <- function(base) {
    h <- .cellsnp_open(dir, base)
    mcon <- h$con()
    on.exit(close(mcon), add = TRUE)
    methods::as(Matrix::readMM(mcon), "TsparseMatrix")
  }
  AD <- read_mtx("cellSNP.tag.AD.mtx")
  DP <- read_mtx("cellSNP.tag.DP.mtx")
  OTH <- if (need_oth) tryCatch(read_mtx("cellSNP.tag.OTH.mtx"), error = function(e) NULL) else NULL

  if (nrow(DP) != nrow(sites)) {
    stop(sprintf("%s: matrix has %d rows but base VCF has %d variants",
                 dir, nrow(DP), nrow(sites)), call. = FALSE)
  }
  if (ncol(DP) != length(barcodes)) {
    stop(sprintf("%s: matrix has %d columns but %d barcodes",
                 dir, ncol(DP), length(barcodes)), call. = FALSE)
  }
  list(sites = sites, barcodes = barcodes, AD = AD, DP = DP, OTH = OTH)
}


#' @title buildCellSNPDatabase: Import cellsnp-lite Output as a SNP Database
#' @name buildCellSNPDatabase
#'
#' @description
#' Reads one or more cellsnp-lite output directories and merges them into the
#' \code{snp_database} structure used by variantCell, taking the union of sites
#' across samples.
#'
#' This is the import path for runs against a non-germline region list - RNA
#' editing sites, chrM, a candidate panel - which have no vireo output because
#' there is no genotype to deconvolve. For the germline genotyping workflow use
#' \code{addSampleData()} plus \code{buildSNPDatabase()}.
#'
#' @param dirs Character vector of cellsnp-lite output directories.
#' @param sample_ids Character vector of sample names, same length as `dirs`.
#' @param cell_metadata Optional data frame of per-cell metadata to join. Must
#'   contain a `cell_id` column (or have cell IDs as rownames) matching the
#'   constructed cell IDs. Columns such as `donor_type` and `cell_type` come
#'   from here; without them only `cell_id` and `sample_id` are available.
#' @param barcode_prefix Character vector, or NULL. Prefix prepended to each
#'   sample's barcodes to build cell IDs. Defaults to \code{paste0(sample_ids, "_")}.
#'   Must reproduce the cell naming used by whatever object supplies
#'   `cell_metadata`.
#' @param min_total_dp Numeric. Drop sites whose summed depth across all cells
#'   and samples falls below this. 0 keeps everything.
#' @param keep_oth Logical. Retain the OTH matrix, which gives a position-matched
#'   sequencing-error floor. Default TRUE.
#' @param require_metadata Logical. If TRUE, cells absent from `cell_metadata`
#'   are dropped. If FALSE (default) they are kept with NA metadata.
#' @param verbose Logical. Print progress.
#'
#' @return A list with `ad_matrix`, `dp_matrix`, `oth_matrix` (or NULL),
#'   `dp_matrix_normalized` (NULL - normalisation is a depth concept and is not
#'   applied here), `snp_info`, `cell_metadata` and `qc_report`.
#'
#' @examples
#' \dontrun{
#' db <- buildCellSNPDatabase(
#'   dirs       = file.path("editing_pilot/run", c("TBX3", "CLAD-2")),
#'   sample_ids = c("TBX3", "CLAD-2"),
#'   cell_metadata = seurat_meta)
#' project <- variantCell$new()
#' project$snp_database <- db
#' project$setProjectIdentity("donor_type")
#' }
#' @export
buildCellSNPDatabase <- function(dirs,
                                 sample_ids,
                                 cell_metadata = NULL,
                                 barcode_prefix = NULL,
                                 min_total_dp = 0,
                                 keep_oth = TRUE,
                                 require_metadata = FALSE,
                                 verbose = TRUE) {
  if (length(dirs) != length(sample_ids)) {
    stop("dirs and sample_ids must be the same length", call. = FALSE)
  }
  if (anyDuplicated(sample_ids)) stop("sample_ids must be unique", call. = FALSE)
  if (is.null(barcode_prefix)) barcode_prefix <- paste0(sample_ids, "_")
  if (length(barcode_prefix) == 1) barcode_prefix <- rep(barcode_prefix, length(dirs))

  # Pass 1: read each sample, collect the site universe.
  samples <- vector("list", length(dirs))
  keys <- character(0)
  for (k in seq_along(dirs)) {
    if (verbose) cat(sprintf("Reading %s ...\n", sample_ids[k]))
    s <- .read_cellsnp_sample(dirs[k], need_oth = keep_oth)
    s$key <- paste(s$sites$CHROM, s$sites$POS, s$sites$REF, s$sites$ALT, sep = ":")
    samples[[k]] <- s
    keys <- c(keys, s$key)
    if (verbose) cat(sprintf("  %d sites x %d cells\n", nrow(s$sites), length(s$barcodes)))
  }

  all_keys <- unique(keys)
  n_sites <- length(all_keys)
  if (verbose) cat(sprintf("Union: %d sites across %d samples\n", n_sites, length(dirs)))

  first <- match(all_keys, keys)
  site_tab <- do.call(rbind, lapply(samples, `[[`, "sites"))
  snp_info <- site_tab[first, , drop = FALSE]
  snp_info$snp_id <- all_keys
  rownames(snp_info) <- NULL

  # Pass 2: remap row indices and accumulate triplets. Assembling once avoids
  # the quadratic cost of cbind-ing a matrix per sample.
  cell_ids <- unlist(lapply(seq_along(samples), function(k)
    paste0(barcode_prefix[k], samples[[k]]$barcodes)), use.names = FALSE)
  cell_smp <- rep(sample_ids, vapply(samples, function(s) length(s$barcodes), integer(1)))
  if (anyDuplicated(cell_ids)) {
    stop("constructed cell IDs are not unique - check barcode_prefix", call. = FALSE)
  }

  col_off <- c(0L, cumsum(vapply(samples, function(s) length(s$barcodes), integer(1))))
  acc <- function(field) {
    i <- j <- integer(0); x <- numeric(0)
    for (k in seq_along(samples)) {
      m <- samples[[k]][[field]]
      if (is.null(m)) next
      rmap <- match(samples[[k]]$key, all_keys)
      i <- c(i, rmap[m@i + 1L])
      j <- c(j, m@j + 1L + col_off[k])
      # A MatrixMarket file declared "pattern" reads back as an n*Matrix with no
      # x slot; every stored entry then means 1. Reading @x blindly errors out.
      x <- c(x, if (methods::.hasSlot(m, "x")) m@x else rep(1, length(m@i)))
    }
    Matrix::sparseMatrix(i = i, j = j, x = x,
                         dims = c(n_sites, length(cell_ids)),
                         dimnames = list(all_keys, cell_ids))
  }

  if (verbose) cat("Assembling matrices...\n")
  DP <- acc("DP")
  AD <- acc("AD")
  OTH <- if (keep_oth && !is.null(samples[[1]]$OTH)) acc("OTH") else NULL

  if (min_total_dp > 0) {
    tot <- Matrix::rowSums(DP)
    keep <- tot >= min_total_dp
    if (verbose) cat(sprintf("Dropping %d sites below total DP %g; %d remain\n",
                             sum(!keep), min_total_dp, sum(keep)))
    DP <- DP[keep, , drop = FALSE]; AD <- AD[keep, , drop = FALSE]
    if (!is.null(OTH)) OTH <- OTH[keep, , drop = FALSE]
    snp_info <- snp_info[keep, , drop = FALSE]
    rownames(snp_info) <- NULL
  }

  meta <- data.frame(cell_id = cell_ids, sample_id = cell_smp, stringsAsFactors = FALSE)
  if (!is.null(cell_metadata)) {
    cm <- as.data.frame(cell_metadata, stringsAsFactors = FALSE)
    if (!"cell_id" %in% colnames(cm)) {
      if (is.null(rownames(cm))) stop("cell_metadata needs a cell_id column or rownames", call. = FALSE)
      cm$cell_id <- rownames(cm)
    }
    idx <- match(meta$cell_id, cm$cell_id)
    hit <- sum(!is.na(idx))
    if (verbose) cat(sprintf("Metadata joined for %d / %d cells (%.1f%%)\n",
                             hit, nrow(meta), 100 * hit / nrow(meta)))
    if (hit == 0) {
      warning("no cells matched cell_metadata - check barcode_prefix", call. = FALSE)
    } else if (hit < 0.5 * nrow(meta)) {
      warning(sprintf("only %.1f%% of cells matched cell_metadata - check barcode_prefix",
                      100 * hit / nrow(meta)), call. = FALSE)
    }
    extra <- setdiff(colnames(cm), c("cell_id", "sample_id"))
    for (col in extra) meta[[col]] <- cm[[col]][idx]
    if (require_metadata) {
      keep <- !is.na(idx)
      meta <- meta[keep, , drop = FALSE]
      DP <- DP[, keep, drop = FALSE]; AD <- AD[, keep, drop = FALSE]
      if (!is.null(OTH)) OTH <- OTH[, keep, drop = FALSE]
      rownames(meta) <- NULL
      if (verbose) cat(sprintf("Dropped %d cells absent from metadata\n", sum(!keep)))
    }
  }

  err <- if (!is.null(OTH)) sum(OTH) / max(1, sum(DP)) / 2 else NA_real_
  list(
    ad_matrix = AD,
    dp_matrix = DP,
    oth_matrix = OTH,
    dp_matrix_normalized = NULL,
    snp_info = snp_info,
    snp_annotations = NULL,
    snp_metrics = NULL,
    cell_metadata = meta,
    qc_report = list(
      source = "cellsnp-lite",
      total_samples = length(dirs),
      total_cells = nrow(meta),
      total_snps = nrow(snp_info),
      total_dp = sum(DP),
      total_ad = sum(AD),
      alt_fraction = sum(AD) / max(1, sum(DP)),
      per_base_error_floor = err,
      signal_to_error = if (is.na(err) || err == 0) NA_real_
                        else (sum(AD) / max(1, sum(DP))) / err,
      cells_per_sample = table(meta$sample_id)
    )
  )
}


#' @title variantCellFromCellSNP: Build a variantCell Object from cellsnp-lite Output
#' @name variantCellFromCellSNP
#'
#' @description
#' Convenience wrapper: calls \code{buildCellSNPDatabase()} and returns a
#' variantCell object with the database attached and an identity optionally set.
#'
#' @param ... Passed to \code{buildCellSNPDatabase()}.
#' @param identity Character or NULL. Metadata column to set as the project
#'   identity.
#'
#' @return A variantCell object.
#' @export
variantCellFromCellSNP <- function(..., identity = NULL) {
  db <- buildCellSNPDatabase(...)
  p <- variantCell$new()
  p$snp_database <- db
  if (!is.null(identity)) {
    if (!identity %in% colnames(db$cell_metadata)) {
      stop(sprintf("identity column '%s' not present in cell metadata", identity), call. = FALSE)
    }
    p$current_project_ident <- identity
  }
  p
}
