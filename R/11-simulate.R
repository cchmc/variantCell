# Synthetic cellSNP-lite / Vireo data generator.
#
# This exists so the vignettes, the examples and the plotting code can run end
# to end with no external data. It writes the same files a real cellsnp-lite
# plus vireo run produces - cellSNP.base.vcf, cellSNP.samples.tsv, the AD/DP/OTH
# MatrixMarket matrices and donor_ids.tsv - so a vignette built on it exercises
# the actual file-parsing import path rather than reaching into the object.
#
# Everything here is deterministic given `seed`, and the draw is made without
# disturbing the caller's random stream.


# Cell types and the compartment each belongs to. inferDonorType() keys on
# compartment, so this table also fixes what counts as "structural" downstream.
.sim_cell_types <- function() {
  data.frame(
    cell_type = c("Alveolar Epithelial Type 2", "Basal", "Ciliated",
                  "Capillary", "Alveolar Fibroblast",
                  "Macrophage", "Classical Monocyte",
                  "CD8+ T", "NK", "Memory B"),
    compartment = c("Epithelial", "Epithelial", "Epithelial",
                    "Endothelial", "Stromal",
                    "Myeloid", "Myeloid",
                    "Lymphoid", "Lymphoid", "Lymphoid"),
    structural = c(TRUE, TRUE, TRUE, TRUE, TRUE,
                   FALSE, FALSE, FALSE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
}

# Genes carrying the simulated sites. Every one returns exactly one record from
# EnsDb.Hsapiens.v86 on a primary autosome - genes with alt-contig or readthrough
# duplicates (B2M, HLA-A, PTPRC, EPCAM, ...) break plotSNPs(), which calls
# start()/end() on the result of a SymbolFilter lookup and would get a vector.
#
# `immune_log2` shifts read depth between immune and structural cells. That is
# the signal findDESNPs() is designed to find: it tests depth at a variant site,
# which is transcript abundance, not genotype.
.sim_genes <- function() {
  data.frame(
    gene = c("SFTPC", "SCGB1A1", "CLDN5", "DCN", "PECAM1", "VIM",
             "JAK1", "IFI6", "CD74", "FCN1", "VCAN", "MARCO"),
    base_depth = c(1.2, 1.0, 0.7, 0.8, 0.75, 1.1,
                   0.6, 0.65, 1.5, 0.9, 0.7, 0.6),
    immune_log2 = c(-2.5, -2.0, -1.8, -1.6, -1.4, 0.0,
                    0.2, 0.4, 2.2, 2.0, 1.6, 1.4),
    stringsAsFactors = FALSE
  )
}

# The default design. Two biopsies from P1 make the "same genome, different
# library" case available, which is the contrast that must come back empty.
.sim_design <- function() {
  data.frame(
    sample    = c("S1", "S2", "S3", "S4", "S5"),
    patient   = c("P1", "P1", "P2", "P3", "P4"),
    condition = c("Rejection", "Rejection", "Rejection", "Stable", "Stable"),
    stringsAsFactors = FALSE
  )
}


# Draw exon positions for each gene from EnsDb, so sites land where RNA-seq
# would actually cover them and annotate_snps() can assign a gene name.
.sim_site_positions <- function(genes_df, snps_per_gene, verbose) {
  if (!requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE) ||
      !requireNamespace("ensembldb", quietly = TRUE)) {
    stop("EnsDb.Hsapiens.v86 and ensembldb are required to simulate site positions.")
  }
  edb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86

  # One query for the whole panel. The formula interface avoids depending on
  # AnnotationFilter::SymbolFilter, which ensembldb does not re-export.
  syms <- genes_df$gene
  all_ex <- ensembldb::exons(edb, filter = ~ symbol %in% syms,
                             columns = c("symbol"))
  all_ex <- all_ex[as.character(GenomicRanges::seqnames(all_ex)) %in%
                     c(1:22, "X")]
  ex_sym <- GenomicRanges::mcols(all_ex)$symbol

  out <- vector("list", nrow(genes_df))
  for (i in seq_len(nrow(genes_df))) {
    sym <- genes_df$gene[i]
    ex <- all_ex[ex_sym == sym]
    if (length(ex) == 0) stop(sprintf("No exons found for gene '%s'.", sym))

    # Widest exons first: sites in a long 3' UTR exon are the ones a real 3'
    # library would actually see.
    ex <- ex[order(-GenomicRanges::width(ex))]
    pos <- integer(0)
    for (j in seq_along(ex)) {
      if (length(pos) >= snps_per_gene) break
      s <- GenomicRanges::start(ex[j]); e <- GenomicRanges::end(ex[j])
      take <- min(snps_per_gene - length(pos), max(1L, floor((e - s) / 10)))
      if (take > 0 && e > s) {
        pos <- c(pos, round(seq(s + 5, e - 5, length.out = take)))
      }
    }
    pos <- sort(unique(pos))[seq_len(min(snps_per_gene, length(unique(pos))))]

    out[[i]] <- data.frame(
      CHROM = as.character(GenomicRanges::seqnames(ex[1])),
      POS   = as.integer(pos),
      gene  = sym,
      base_depth  = genes_df$base_depth[i],
      immune_log2 = genes_df$immune_log2[i],
      stringsAsFactors = FALSE
    )
  }
  sites <- do.call(rbind, out)
  sites <- sites[order(as.integer(factor(sites$CHROM, levels = c(1:22, "X"))),
                       sites$POS), ]
  rownames(sites) <- NULL
  if (verbose) cat(sprintf("Simulated %d sites across %d genes\n",
                           nrow(sites), nrow(genes_df)))
  sites
}


# cellsnp-lite writes integer MatrixMarket. Writing it out by hand rather than
# via Matrix::writeMM is both more faithful and avoids the pattern-matrix trap:
# writeMM emits a "pattern" header whenever every stored value is 1, and readMM
# hands that back as an ngTMatrix with no x slot.
.sim_write_mtx <- function(file, i, j, x, nrow, ncol) {
  keep <- x > 0
  i <- i[keep]; j <- j[keep]; x <- x[keep]
  ord <- order(j, i)
  con <- file(file, "wt")
  on.exit(close(con), add = TRUE)
  writeLines(c("%%MatrixMarket matrix coordinate integer general",
               "%",
               sprintf("%d %d %d", nrow, ncol, length(i))), con)
  if (length(i) > 0) {
    writeLines(sprintf("%d %d %d", i[ord], j[ord], x[ord]), con)
  }
  invisible(file)
}


#' Simulate a cellSNP-lite plus Vireo dataset
#'
#' Writes a complete synthetic multi-genome experiment to disk in exactly the
#' layout `addSampleData()` and `inferDonorType()` expect, and returns the paths
#' and cell metadata needed to load it. Intended for vignettes, examples,
#' reproducible bug reports, and checking an installation end to end.
#'
#' Each simulated library contains two genetically distinct individuals. The
#' first supplies mostly structural cells and the second mostly immune cells,
#' reproducing the near-total lineage segregation seen in real transplant data -
#' which is what makes `inferDonorType()` work.
#'
#' Which individual Vireo calls `donor0` is randomised per sample, because that
#' is what Vireo actually does: the labels are arbitrary and are reassigned on
#' every run. The ground truth is returned in `$truth` so a vignette can show the
#' consequence rather than assert it.
#'
#' @param path Directory to write to. Default: a new session temp directory.
#' @param design Data frame with columns `sample`, `patient` and `condition`.
#'   Samples sharing a `patient` share both genomes. Default: five samples
#'   across four patients, two of them repeat biopsies of the same patient.
#' @param n_cells Integer. Cells per sample. Default 160.
#' @param snps_per_gene Integer. Sites per gene. Default 180.
#' @param genes Data frame of genes with `gene`, `base_depth` and `immune_log2`
#'   columns, or NULL for the built-in panel of twelve lung and immune genes.
#' @param seed Integer seed. The draw does not disturb the caller's RNG stream.
#' @param verbose Logical. Print progress. Default TRUE.
#'
#' @return A list with `root`, `cellsnp_paths`, `vireo_paths`, `vireo_dir`,
#'   `prefixes`, `metadata` (row names are prefixed cell IDs, ready for
#'   `addSampleData(data_type = "dataframe")`), `design`, `sites`, and `truth` -
#'   the genome behind each Vireo label, plus the per-genome genotype matrix.
#'
#' @examples
#' sim <- simulateVariantCellData(n_cells = 40, snps_per_gene = 10,
#'                                verbose = FALSE)
#' list.files(sim$cellsnp_paths[["S1"]])
#' head(sim$truth$labels)
#'
#' @export
simulateVariantCellData <- function(path = NULL,
                                    design = NULL,
                                    n_cells = 160,
                                    snps_per_gene = 180,
                                    genes = NULL,
                                    seed = 1,
                                    verbose = TRUE) {

  if (is.null(design)) design <- .sim_design()
  if (is.null(genes))  genes  <- .sim_genes()
  if (!all(c("sample", "patient", "condition") %in% colnames(design))) {
    stop("design must have columns: sample, patient, condition")
  }
  if (anyDuplicated(design$sample) > 0) stop("design$sample must be unique")
  if (is.null(path)) path <- file.path(tempdir(), "variantCell_sim")
  dir.create(path, recursive = TRUE, showWarnings = FALSE)

  # Reproducible, but the caller's stream is left exactly as it was found.
  if (exists(".Random.seed", envir = globalenv())) {
    old_seed <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = globalenv())), add = TRUE)
  }
  set.seed(seed)

  sites <- .sim_site_positions(genes, snps_per_gene, verbose)
  n_sites <- nrow(sites)
  bases <- c("A", "C", "G", "T")
  sites$REF <- sample(bases, n_sites, replace = TRUE)
  sites$ALT <- vapply(sites$REF, function(r) sample(setdiff(bases, r), 1), "")

  # One genotype per genome per site, drawn under Hardy-Weinberg from an allele
  # frequency above 0.05 - the same constraint the 1000G panel this package is
  # normally run against imposes, so every simulated variant is common.
  af <- stats::runif(n_sites, 0.05, 0.95)
  ctypes <- .sim_cell_types()
  patients <- unique(design$patient)
  genome_ids <- c(paste0(patients, "_graft"), paste0(patients, "_host"))
  genotypes <- vapply(genome_ids, function(g) {
    stats::rbinom(n_sites, 2, af)
  }, integer(n_sites))
  rownames(genotypes) <- paste(sites$CHROM, sites$POS, sites$REF, sites$ALT,
                               sep = "_")

  # Expected alt fraction. The 0.01/0.99 rather than 0/1 is sequencing error,
  # and it matters: without it a homozygous-reference site separates perfectly
  # and every downstream statistic looks better than it should.
  alt_p <- c(0.01, 0.5, 0.99)

  cellsnp_paths <- character(0)
  vireo_paths   <- character(0)
  prefixes      <- character(0)
  meta_list     <- list()
  label_rows    <- list()

  vireo_root <- file.path(path, "vireo")
  dir.create(vireo_root, showWarnings = FALSE)

  for (k in seq_len(nrow(design))) {
    sid <- design$sample[k]
    pid <- design$patient[k]
    prefix <- paste0(sid, "_")

    # Graft cells are mostly structural, host cells mostly immune.
    n_graft <- round(n_cells * 0.45)
    n_host  <- n_cells - n_graft
    pick <- function(n, structural) {
      pool <- ctypes$cell_type[ctypes$structural == structural]
      other <- ctypes$cell_type[ctypes$structural != structural]
      # ~8% of cells come from the opposite compartment, as in real data.
      n_off <- stats::rbinom(1, n, 0.08)
      c(sample(pool, n - n_off, replace = TRUE),
        sample(other, n_off, replace = TRUE))
    }
    cell_types <- c(pick(n_graft, TRUE), pick(n_host, FALSE))
    cell_genome <- c(rep(paste0(pid, "_graft"), n_graft),
                     rep(paste0(pid, "_host"), n_host))

    barcodes <- paste0(
      vapply(seq_len(n_cells), function(i)
        paste(sample(bases, 16, replace = TRUE), collapse = ""), ""),
      "-1")
    while (anyDuplicated(barcodes) > 0) {
      dup <- which(duplicated(barcodes))
      barcodes[dup] <- paste0(
        vapply(dup, function(i)
          paste(sample(bases, 16, replace = TRUE), collapse = ""), ""), "-1")
    }

    ord <- sample(n_cells)
    barcodes <- barcodes[ord]; cell_types <- cell_types[ord]
    cell_genome <- cell_genome[ord]
    compartment <- ctypes$compartment[match(cell_types, ctypes$cell_type)]
    is_immune <- !ctypes$structural[match(cell_types, ctypes$cell_type)]

    # ---- draw counts -------------------------------------------------------
    # lambda[s, c] = base_depth[s] * 2^(immune_log2[s] * is_immune[c]). The
    # shift is per site AND per cell, so it needs the outer product; recycling
    # two vectors of different lengths silently gives an n_sites x n_sites grid.
    lambda <- sites$base_depth *
      2^outer(sites$immune_log2, as.numeric(is_immune))
    dp <- matrix(stats::rpois(length(lambda), lambda),
                 nrow = n_sites, ncol = n_cells)
    nz <- which(dp > 0)
    ii <- ((nz - 1L) %% n_sites) + 1L
    jj <- ((nz - 1L) %/% n_sites) + 1L
    dpv <- dp[nz]
    gt <- genotypes[cbind(ii, match(cell_genome[jj], genome_ids))]
    adv <- stats::rbinom(length(dpv), dpv, alt_p[gt + 1L])
    # OTH is the two bases that are neither REF nor ALT: pure sequencing error,
    # and the position-matched floor the alt signal is judged against.
    othv <- stats::rbinom(length(dpv), dpv, 0.002)

    cdir <- file.path(path, paste0(sid, "_cellsnp"))
    dir.create(cdir, showWarnings = FALSE)
    .sim_write_mtx(file.path(cdir, "cellSNP.tag.DP.mtx"), ii, jj, dpv, n_sites, n_cells)
    .sim_write_mtx(file.path(cdir, "cellSNP.tag.AD.mtx"), ii, jj, adv, n_sites, n_cells)
    .sim_write_mtx(file.path(cdir, "cellSNP.tag.OTH.mtx"), ii, jj, othv, n_sites, n_cells)
    writeLines(barcodes, file.path(cdir, "cellSNP.samples.tsv"))

    tot_ad <- tapply(adv, ii, sum)[as.character(seq_len(n_sites))]
    tot_dp <- tapply(dpv, ii, sum)[as.character(seq_len(n_sites))]
    tot_ad[is.na(tot_ad)] <- 0; tot_dp[is.na(tot_dp)] <- 0
    writeLines(c(
      "##fileformat=VCFv4.2",
      "##source=simulateVariantCellData",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
      sprintf("%s\t%d\t.\t%s\t%s\t.\t.\tAD=%d;DP=%d",
              sites$CHROM, sites$POS, sites$REF, sites$ALT,
              as.integer(tot_ad), as.integer(tot_dp))),
      file.path(cdir, "cellSNP.base.vcf"))

    # ---- vireo output ------------------------------------------------------
    # The donor0/donor1 assignment is deliberately arbitrary per sample.
    flip <- stats::runif(1) < 0.5
    graft_label <- if (flip) "donor1" else "donor0"
    host_label  <- if (flip) "donor0" else "donor1"
    donor_id <- ifelse(grepl("_graft$", cell_genome), graft_label, host_label)
    # A small number of genotype doublets, as any real run has.
    n_dbl <- stats::rbinom(1, n_cells, 0.02)
    if (n_dbl > 0) donor_id[sample(n_cells, n_dbl)] <- "doublet"

    vdir <- file.path(vireo_root, paste0(sid, "_vireo"))
    dir.create(vdir, showWarnings = FALSE)
    vireo_df <- data.frame(
      cell = barcodes,
      donor_id = donor_id,
      prob_max = round(stats::runif(n_cells, 0.9, 1), 4),
      prob_doublet = round(stats::runif(n_cells, 0, 0.05), 4),
      n_vars = as.integer(Matrix::colSums(dp > 0)),
      stringsAsFactors = FALSE
    )
    utils::write.table(vireo_df, file.path(vdir, "donor_ids.tsv"),
                       sep = "\t", quote = FALSE, row.names = FALSE)

    meta_list[[sid]] <- data.frame(
      Barcode = barcodes,
      Biopsy = sid,
      patient = pid,
      condition = design$condition[k],
      cell_type = cell_types,
      compartment = compartment,
      row.names = paste0(prefix, barcodes),
      stringsAsFactors = FALSE
    )
    label_rows[[sid]] <- data.frame(
      sample = sid, patient = pid,
      donor0 = if (flip) "host" else "graft",
      donor1 = if (flip) "graft" else "host",
      stringsAsFactors = FALSE
    )

    cellsnp_paths[sid] <- cdir
    vireo_paths[sid]   <- file.path(vdir, "donor_ids.tsv")
    prefixes[sid]      <- prefix

    if (verbose) {
      cat(sprintf("  %s (%s, %s): %d cells, graft = %s\n",
                  sid, pid, design$condition[k], n_cells, graft_label))
    }
  }

  list(
    root = path,
    cellsnp_paths = cellsnp_paths,
    vireo_paths = vireo_paths,
    vireo_dir = vireo_root,
    prefixes = prefixes,
    metadata = do.call(rbind, unname(meta_list)),
    design = design,
    sites = sites,
    truth = list(labels = do.call(rbind, unname(label_rows)),
                 genotypes = genotypes,
                 allele_freq = af)
  )
}
