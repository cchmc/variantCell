#' @title inferDonorType: Infer Donor vs Recipient from Vireo Assignments
#' @name inferDonorType
#'
#' @description
#' Infers which Vireo genetic donor corresponds to the transplanted organ
#' ("Donor") and which corresponds to the transplant recipient ("Recipient"),
#' using the cell-type composition of each genotype cluster.
#'
#' In a transplant biopsy the graft contributes the structural tissue
#' (epithelium, endothelium, stroma), while the recipient contributes
#' infiltrating immune cells. The genotype cluster holding the larger share of
#' structural cells is therefore the donor.
#'
#' @details
#' **Why this function is necessary.** Vireo's `donor0` / `donor1` labels are
#' assigned arbitrarily and independently on every run. They carry no
#' biological meaning and are *not* stable across re-runs of the same sample:
#' re-processing a sample (for example after a CellRanger upgrade) will
#' frequently swap them. A donor/recipient mapping recorded against one Vireo
#' run must never be reused against another. Always re-infer.
#'
#' **Validation.** On a 22-sample lung transplant cohort with an independently
#' established mapping, this heuristic reproduced the known assignment for
#' 22/22 samples, with a median separation in structural fraction of >0.97.
#'
#' **Interpreting `margin`.** `margin` is the absolute difference in structural
#' fraction between the two genotype clusters. Well-separated samples score
#' near 1.0. Samples below `min_margin` are reported as `"ambiguous"` and are
#' assigned `NA` rather than guessed, because an inverted call silently swaps
#' donor and recipient for every downstream comparison. Always inspect the
#' returned `mapping` table before use.
#'
#' @param object A Seurat object, SingleCellExperiment, or a data frame of cell
#'   metadata.
#' @param vireo_dir Character. Directory holding the per-sample Vireo output
#'   folders.
#' @param sample_col Character. Metadata column identifying the sample. Must
#'   match the Vireo folder names (after `vireo_suffix`). Default: "Biopsy".
#' @param barcode_col Character. Metadata column holding the plain cell barcode
#'   as it appears in Vireo's `donor_ids.tsv` (e.g. "AAACCCAAGATGCGAC-1"). For a
#'   merged object this is usually *not* `colnames()`. Default: "Barcode".
#' @param structural_col Character. Metadata column used to identify structural
#'   cells. Default: "compartment".
#' @param structural_values Character vector. Values of `structural_col` that
#'   count as graft-derived structural tissue.
#'   Default: c("Epithelial", "Endothelial", "Stromal").
#' @param vireo_suffix Character. Appended to the sample name to locate its
#'   Vireo folder (e.g. "_vireo" gives "TBX1_vireo"). Default: "_vireo".
#' @param samples Character vector. Restrict to these samples. Default: NULL
#'   (all samples found in `sample_col`).
#' @param min_cells Integer. Minimum matched cells required to call a sample.
#'   Default: 50.
#' @param min_margin Numeric. Minimum separation in structural fraction needed
#'   to make a call; below this the sample is flagged "ambiguous" and left NA.
#'   Default: 0.1.
#' @param add_metadata Logical. Attach `vireo_donor` and `donor_type` columns to
#'   `object` and return it. Default: TRUE.
#' @param verbose Logical. Print progress and a summary table. Default: TRUE.
#'
#' @return A list with:
#'   \item{mapping}{Data frame, one row per sample: structural fractions, cell
#'     counts, `margin`, `call` ("ok"/"ambiguous"/"skipped"), and the inferred
#'     `donor0`/`donor1` roles.}
#'   \item{donor_types}{Named list, one entry per confidently called sample,
#'     each a named vector `c(donor0 = ..., donor1 = ...)` ready to pass
#'     straight to `addSampleData(donor_type = ...)`.}
#'   \item{assignments}{Data frame of per-cell assignments (sample, barcode,
#'     vireo donor, inferred donor_type).}
#'   \item{object}{The input object with metadata columns added, when
#'     `add_metadata = TRUE`.}
#'
#' @examples
#' \dontrun{
#' inf <- inferDonorType(merged_obj, vireo_dir = "Vireo/cellranger8")
#'
#' # ALWAYS review before use - especially any row flagged "ambiguous"
#' inf$mapping
#'
#' # feed a confidently-called sample straight into addSampleData()
#' project$addSampleData(
#'   sample_id   = "TBX1",
#'   vireo_path  = "Vireo/cellranger8/TBX1_vireo",
#'   cellsnp_path= "cellSNP/cellranger8/TBX1",
#'   cell_data   = inf$object,
#'   donor_type  = inf$donor_types[["TBX1"]],
#'   prefix_text = "TBX1_"
#' )
#' }
#'
#' @export
inferDonorType <- function(object,
                           vireo_dir,
                           sample_col = "Biopsy",
                           barcode_col = "Barcode",
                           structural_col = "compartment",
                           structural_values = c("Epithelial", "Endothelial", "Stromal"),
                           vireo_suffix = "_vireo",
                           samples = NULL,
                           min_cells = 50,
                           min_margin = 0.1,
                           add_metadata = TRUE,
                           verbose = TRUE) {

  # ---- extract metadata ----------------------------------------------------
  meta <- .extract_cell_metadata(object)

  for (col in c(sample_col, barcode_col, structural_col)) {
    if (!col %in% colnames(meta)) {
      stop(sprintf("Column '%s' not found in cell metadata. Available: %s",
                   col, paste(colnames(meta), collapse = ", ")))
    }
  }
  if (!dir.exists(vireo_dir)) {
    stop(sprintf("vireo_dir does not exist: %s", vireo_dir))
  }

  meta[[sample_col]]  <- as.character(meta[[sample_col]])
  meta[[barcode_col]] <- as.character(meta[[barcode_col]])

  all_samples <- sort(unique(meta[[sample_col]]))
  if (!is.null(samples)) {
    missing <- setdiff(samples, all_samples)
    if (length(missing) > 0) {
      warning("Samples not present in metadata: ", paste(missing, collapse = ", "))
    }
    all_samples <- intersect(samples, all_samples)
  }
  if (length(all_samples) == 0) stop("No samples to process.")

  if (verbose) {
    cat("=== Inferring donor/recipient identity ===\n")
    cat(sprintf("Samples: %d | structural = %s\n",
                length(all_samples), paste(structural_values, collapse = ", ")))
  }

  mapping <- list()
  assignments <- list()

  for (s in all_samples) {
    vpath <- file.path(vireo_dir, paste0(s, vireo_suffix), "donor_ids.tsv")

    row <- data.frame(sample = s, n_matched = 0L, n_donor0 = 0L, n_donor1 = 0L,
                      struct_frac_donor0 = NA_real_, struct_frac_donor1 = NA_real_,
                      margin = NA_real_, donor0 = NA_character_,
                      donor1 = NA_character_, call = "skipped",
                      note = NA_character_, stringsAsFactors = FALSE)

    if (!file.exists(vpath)) {
      row$note <- "donor_ids.tsv not found"
      mapping[[s]] <- row
      if (verbose) cat(sprintf("  %-8s SKIPPED - %s\n", s, row$note))
      next
    }

    vireo <- utils::read.delim(vpath, stringsAsFactors = FALSE)
    if (!all(c("cell", "donor_id") %in% colnames(vireo))) {
      row$note <- "unexpected donor_ids.tsv format"
      mapping[[s]] <- row
      if (verbose) cat(sprintf("  %-8s SKIPPED - %s\n", s, row$note))
      next
    }

    smeta <- meta[meta[[sample_col]] == s, , drop = FALSE]

    # Positional lookup keyed on barcode; merge() would reorder the rows.
    idx <- match(vireo$cell, smeta[[barcode_col]])
    vireo$.compartment <- smeta[[structural_col]][idx]
    vireo$.matched <- !is.na(idx)

    singlets <- vireo[vireo$donor_id %in% c("donor0", "donor1") & vireo$.matched, ,
                      drop = FALSE]

    present <- unique(vireo$donor_id[vireo$donor_id %in%
                                       grep("^donor[0-9]+$", unique(vireo$donor_id),
                                            value = TRUE)])
    if (length(present) > 2) {
      row$note <- sprintf("%d genotype clusters found (expected 2)", length(present))
      mapping[[s]] <- row
      if (verbose) cat(sprintf("  %-8s SKIPPED - %s\n", s, row$note))
      next
    }

    row$n_matched <- nrow(singlets)
    row$n_donor0  <- sum(singlets$donor_id == "donor0")
    row$n_donor1  <- sum(singlets$donor_id == "donor1")

    if (nrow(singlets) < min_cells) {
      row$note <- sprintf("only %d matched cells (min_cells = %d)",
                          nrow(singlets), min_cells)
      mapping[[s]] <- row
      if (verbose) cat(sprintf("  %-8s SKIPPED - %s\n", s, row$note))
      next
    }

    sfrac <- sapply(c("donor0", "donor1"), function(d) {
      sub <- singlets[singlets$donor_id == d, , drop = FALSE]
      if (nrow(sub) == 0) return(NA_real_)
      mean(sub$.compartment %in% structural_values)
    })

    row$struct_frac_donor0 <- unname(sfrac["donor0"])
    row$struct_frac_donor1 <- unname(sfrac["donor1"])

    if (any(is.na(sfrac))) {
      row$note <- "one genotype cluster has no matched cells"
      mapping[[s]] <- row
      if (verbose) cat(sprintf("  %-8s SKIPPED - %s\n", s, row$note))
      next
    }

    row$margin <- abs(sfrac["donor0"] - sfrac["donor1"])

    if (row$margin < min_margin) {
      # Refuse to guess: an inverted call silently swaps donor and recipient
      # for every downstream comparison.
      row$call <- "ambiguous"
      row$note <- sprintf("margin %.3f < min_margin %.3f", row$margin, min_margin)
      mapping[[s]] <- row
      if (verbose) cat(sprintf("  %-8s AMBIGUOUS - %s\n", s, row$note))
      next
    }

    donor_is_0 <- sfrac["donor0"] > sfrac["donor1"]
    row$donor0 <- if (donor_is_0) "Donor" else "Recipient"
    row$donor1 <- if (donor_is_0) "Recipient" else "Donor"
    row$call   <- "ok"
    mapping[[s]] <- row

    lab <- c(donor0 = row$donor0, donor1 = row$donor1)
    keep <- vireo$donor_id %in% c("donor0", "donor1")
    assignments[[s]] <- data.frame(
      sample      = s,
      barcode     = vireo$cell,
      vireo_donor = vireo$donor_id,
      donor_type  = ifelse(keep, unname(lab[vireo$donor_id]),
                           ifelse(vireo$donor_id == "doublet", "Doublet", "Unassigned")),
      stringsAsFactors = FALSE
    )

    if (verbose) {
      cat(sprintf("  %-8s donor0=%-9s (struct %.3f, n=%d) | donor1=%-9s (struct %.3f, n=%d) | margin %.3f\n",
                  s, row$donor0, sfrac["donor0"], row$n_donor0,
                  row$donor1, sfrac["donor1"], row$n_donor1, row$margin))
    }
  }

  mapping_df <- do.call(rbind, mapping)
  rownames(mapping_df) <- NULL

  assign_df <- if (length(assignments) > 0) {
    do.call(rbind, assignments)
  } else {
    data.frame(sample = character(0), barcode = character(0),
               vireo_donor = character(0), donor_type = character(0),
               stringsAsFactors = FALSE)
  }
  rownames(assign_df) <- NULL

  donor_types <- lapply(which(mapping_df$call == "ok"), function(i) {
    c(donor0 = mapping_df$donor0[i], donor1 = mapping_df$donor1[i])
  })
  names(donor_types) <- mapping_df$sample[mapping_df$call == "ok"]

  if (verbose) {
    n_ok <- sum(mapping_df$call == "ok")
    n_amb <- sum(mapping_df$call == "ambiguous")
    n_skip <- sum(mapping_df$call == "skipped")
    cat(sprintf("\nCalled %d/%d samples (%d ambiguous, %d skipped)\n",
                n_ok, nrow(mapping_df), n_amb, n_skip))
    if (n_ok > 0) {
      cat(sprintf("Margin: min %.3f, median %.3f\n",
                  min(mapping_df$margin[mapping_df$call == "ok"]),
                  stats::median(mapping_df$margin[mapping_df$call == "ok"])))
      flipped <- mapping_df$sample[mapping_df$call == "ok" & mapping_df$donor0 == "Recipient"]
      cat(sprintf("donor0 = Recipient in %d sample(s)%s\n", length(flipped),
                  if (length(flipped)) paste0(": ", paste(flipped, collapse = ", ")) else ""))
    }
    if (n_amb + n_skip > 0) {
      cat("Review these before use:\n")
      print(mapping_df[mapping_df$call != "ok",
                       c("sample", "call", "margin", "note")], row.names = FALSE)
    }
    cat("\nVireo donor labels are arbitrary per run - re-infer after any re-run.\n")
  }

  result <- list(mapping = mapping_df,
                 donor_types = donor_types,
                 assignments = assign_df)

  if (add_metadata) {
    result$object <- .attach_donor_type(object, assign_df, sample_col, barcode_col)
  }

  return(result)
}


# Pull a cell-metadata data frame out of whatever container was passed in.
.extract_cell_metadata <- function(object) {
  if (is.data.frame(object)) return(object)

  if (inherits(object, "Seurat")) {
    return(object@meta.data)
  }

  if (inherits(object, "SingleCellExperiment")) {
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      stop("SingleCellExperiment package required to read this object.")
    }
    return(as.data.frame(SummarizedExperiment::colData(object)))
  }

  stop("object must be a Seurat object, SingleCellExperiment, or data frame.")
}


# Write vireo_donor / donor_type back onto the object, matched on
# sample + barcode so a merged object with repeated barcodes stays correct.
.attach_donor_type <- function(object, assign_df, sample_col, barcode_col) {

  meta <- .extract_cell_metadata(object)
  key_meta   <- paste(as.character(meta[[sample_col]]),
                      as.character(meta[[barcode_col]]), sep = "|")
  key_assign <- paste(assign_df$sample, assign_df$barcode, sep = "|")

  idx <- match(key_meta, key_assign)
  vireo_donor <- assign_df$vireo_donor[idx]
  donor_type  <- assign_df$donor_type[idx]

  if (is.data.frame(object) && !inherits(object, c("Seurat", "SingleCellExperiment"))) {
    object$vireo_donor <- vireo_donor
    object$donor_type  <- donor_type
    return(object)
  }

  if (inherits(object, "Seurat")) {
    object@meta.data$vireo_donor <- vireo_donor
    object@meta.data$donor_type  <- donor_type
    return(object)
  }

  if (inherits(object, "SingleCellExperiment")) {
    SummarizedExperiment::colData(object)$vireo_donor <- vireo_donor
    SummarizedExperiment::colData(object)$donor_type  <- donor_type
    return(object)
  }

  object
}
