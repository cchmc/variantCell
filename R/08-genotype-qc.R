# ---------------------------------------------------------------------------
# Genotype concordance QC.
#
# Presence/absence tests of germline variants are genotype tests. They are
# informative exactly when the two groups being compared are different genomes,
# and they are structurally empty when the groups share a genome, because
# genotype is fixed within an individual. Nothing in a cell-identity label tells
# the software which situation it is in, so this module measures it from the
# data before any test runs.
#
# Calibration on the 32-sample lung transplant cohort: pseudobulk genotype
# concordance runs ~0.99 between compartments of the same individual and ~0.63
# between unrelated individuals, with donor/recipient discordance ~38%.
# ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.call_pseudobulk_genotypes <- function(AD, DP, min_dp, hom_ref, hom_alt) {
  FR <- ifelse(DP >= min_dp, AD / DP, NA_real_)
  ifelse(is.na(FR), NA_integer_,
         ifelse(FR < hom_ref, 0L, ifelse(FR > hom_alt, 2L, 1L)))
}

#' @keywords internal
#' @noRd
#' Count distinct genomes among a set of sub-pseudobulks by single-linkage
#' clustering on the concordance matrix. Derived from the data, so no patient
#' column is required and mislabelled samples cannot inflate the count.
.count_distinct_genomes <- function(conc, members, same_genome_threshold) {
  if (length(members) == 0) return(NA_integer_)
  if (length(members) == 1) return(1L)
  m <- conc[members, members, drop = FALSE]
  parent <- seq_along(members)
  find <- function(i) { while (parent[i] != i) i <- parent[i]; i }
  for (i in seq_along(members)) for (j in seq_along(members)) {
    if (i < j && !is.na(m[i, j]) && m[i, j] >= same_genome_threshold) {
      ri <- find(i); rj <- find(j)
      if (ri != rj) parent[rj] <- ri
    }
  }
  length(unique(vapply(seq_along(members), find, integer(1))))
}

#' @keywords internal
#' @noRd
#' Smallest attainable two-sided Fisher p for a perfectly separating variant
#' between n1 and n2 independent individuals. This is a permutation-count
#' ceiling: it depends only on group sizes, so no amount of sequencing depth,
#' cell number or statistical refinement can beat it.
.min_attainable_p <- function(n1, n2) {
  if (is.na(n1) || is.na(n2) || n1 < 1 || n2 < 1) return(NA_real_)
  min(1, 2 / choose(n1 + n2, n1))
}

#' @keywords internal
#' @noRd
.pairwise_concordance <- function(G, min_shared_sites) {
  n <- ncol(G)
  conc <- matrix(NA_real_, n, n, dimnames = list(colnames(G), colnames(G)))
  nsh  <- matrix(NA_integer_, n, n, dimnames = list(colnames(G), colnames(G)))
  for (i in seq_len(n)) for (j in i:n) {
    both <- !is.na(G[, i]) & !is.na(G[, j])
    if (sum(both) >= min_shared_sites) {
      conc[i, j] <- conc[j, i] <- mean(G[both, i] == G[both, j])
      nsh[i, j]  <- nsh[j, i]  <- sum(both)
    }
  }
  list(concordance = conc, n_shared = nsh)
}

#' @title Check whether two cell groups are genetically distinct
#' @name checkGenotypeConcordance
#'
#' @description
#' Measures pseudobulk genotype concordance between two sets of cells, and
#' within each set, to establish whether a presence/absence contrast between
#' them is a genotype comparison at all.
#'
#' Germline presence/absence is a genotype test. It is informative when the two
#' groups are different genomes (donor versus recipient within a sample, or two
#' individuals), structurally empty when they share a genome, and confounded
#' when either group pools several genomes -- in that last case hits reflect
#' which individuals landed in which group, not the grouping variable.
#'
#' @param snp_database The \code{snp_database} element of a variantCell project.
#' @param cells1,cells2 Cell identifiers (matching \code{cell_metadata$cell_id})
#'   or integer column indices defining the two groups.
#' @param split_by Optional vector, one entry per cell in \code{cell_metadata},
#'   used to split each group into sub-pseudobulks so within-group genome
#'   heterogeneity can be measured. Defaults to \code{sample_id} when present.
#' @param n_sites Number of best-covered autosomal sites to use. Default 30000.
#' @param min_dp Minimum pseudobulk depth for a site to be called. Default 10.
#' @param min_cells Minimum cells for a sub-pseudobulk to be used. Default 50.
#' @param min_shared_sites Minimum jointly callable sites for a comparison.
#'   Default 500.
#' @param hom_ref,hom_alt Alt-fraction bounds separating homozygous reference,
#'   heterozygous and homozygous alternative. Defaults 0.15 and 0.85.
#' @param same_genome_threshold Concordance at or above which two pseudobulks
#'   are the same individual. Default 0.90.
#' @param different_genome_threshold Concordance at or below which they are
#'   different individuals. Default 0.80.
#' @param verbose Print a summary. Default TRUE.
#'
#' @return A list with \code{concordance} (between the two groups),
#'   \code{n_shared_sites}, \code{within1}/\code{within2} (minimum concordance
#'   among each group's own sub-pseudobulks, NA if not splittable),
#'   \code{verdict}, and \code{message}. \code{verdict} is one of
#'   \code{"different_genomes"}, \code{"same_genome"},
#'   \code{"heterogeneous_groups"}, or \code{"indeterminate"}.
#'
#' @details
#' Only autosomal sites are used, so the result is unaffected by sex. Sites are
#' ranked by total depth across the whole database and the best-covered
#' \code{n_sites} retained, which keeps the aggregation cheap and the calls
#' well supported.
#'
#' A \code{"heterogeneous_groups"} verdict means at least one group contains
#' more than one individual. Presence/absence in that setting compares allele
#' frequencies between two arbitrary sets of people; at small n essentially
#' every hit is chance or ancestry, and increasing power makes it worse rather
#' than better.
#'
#' @export
checkGenotypeConcordance <- function(snp_database,
                                     cells1,
                                     cells2,
                                     split_by = NULL,
                                     n_sites = 30000,
                                     min_dp = 10,
                                     min_cells = 50,
                                     min_shared_sites = 500,
                                     hom_ref = 0.15,
                                     hom_alt = 0.85,
                                     same_genome_threshold = 0.90,
                                     different_genome_threshold = 0.80,
                                     verbose = TRUE) {

  if (is.null(snp_database$ad_matrix) || is.null(snp_database$dp_matrix)) {
    stop("snp_database must contain ad_matrix and dp_matrix")
  }
  cm <- snp_database$cell_metadata
  info <- snp_database$snp_info

  resolve <- function(x, nm) {
    if (is.character(x)) {
      idx <- match(x, cm$cell_id)
      if (anyNA(idx)) stop(sprintf("%d cells in '%s' not found in cell_metadata",
                                   sum(is.na(idx)), nm))
      idx
    } else {
      as.integer(x)
    }
  }
  i1 <- resolve(cells1, "cells1")
  i2 <- resolve(cells2, "cells2")

  if (length(intersect(i1, i2)) > 0) {
    stop("cells1 and cells2 overlap; groups must be disjoint")
  }

  # Autosomal sites only, best covered first. Sex chromosomes are excluded so a
  # sex difference between the groups cannot masquerade as genotype distance.
  dp_tot <- Matrix::rowSums(snp_database$dp_matrix)
  keep <- which(info$CHROM %in% as.character(1:22) & dp_tot > 0)
  if (length(keep) == 0) stop("no autosomal sites with coverage in snp_database")
  keep <- sort(keep[order(dp_tot[keep], decreasing = TRUE)][seq_len(min(n_sites, length(keep)))])

  # Sub-pseudobulks: each group split by split_by, so within-group genome
  # heterogeneity is measurable at the same time as the between-group distance.
  if (is.null(split_by) && "sample_id" %in% names(cm)) split_by <- cm$sample_id
  sel <- c(i1, i2)
  grp <- rep(c("group1", "group2"), c(length(i1), length(i2)))
  sub <- if (is.null(split_by)) grp else paste(grp, split_by[sel], sep = "|")

  tab <- table(sub)
  usable <- names(tab)[tab >= min_cells]
  if (length(usable) == 0) {
    stop(sprintf("no sub-group reaches min_cells = %d (largest has %d)",
                 min_cells, max(tab)))
  }
  drop_n <- sum(tab[!names(tab) %in% usable])
  ok <- sub %in% usable
  sel <- sel[ok]; grp <- grp[ok]; sub <- sub[ok]

  if (!all(c("group1", "group2") %in% unique(grp))) {
    stop(sprintf("one group has no sub-pseudobulk with at least %d cells", min_cells))
  }

  # Aggregate every sub-pseudobulk plus the two whole groups in one pass.
  levs <- c(unique(sub), "group1", "group2")
  n_sel <- length(sel)
  ii <- c(seq_len(n_sel), seq_len(n_sel))
  jj <- c(match(sub, levs), match(grp, levs))
  M <- Matrix::sparseMatrix(i = ii, j = jj, x = 1,
                            dims = c(n_sel, length(levs)),
                            dimnames = list(NULL, levs))
  AD <- as.matrix(snp_database$ad_matrix[keep, sel, drop = FALSE] %*% M)
  DP <- as.matrix(snp_database$dp_matrix[keep, sel, drop = FALSE] %*% M)
  G <- .call_pseudobulk_genotypes(AD, DP, min_dp, hom_ref, hom_alt)

  pw <- .pairwise_concordance(G, min_shared_sites)
  between <- pw$concordance["group1", "group2"]
  n_shared <- pw$n_shared["group1", "group2"]

  # Minimum concordance among a group's own sub-pseudobulks: if the group holds
  # one individual this stays high, if it pools individuals it drops.
  within_min <- function(tag) {
    s <- setdiff(grep(paste0("^", tag, "\\|"), levs, value = TRUE), c("group1", "group2"))
    if (length(s) < 2) return(NA_real_)
    m <- pw$concordance[s, s, drop = FALSE]
    diag(m) <- NA
    if (all(is.na(m))) NA_real_ else min(m, na.rm = TRUE)
  }
  w1 <- within_min("group1")
  w2 <- within_min("group2")

  # Number of independent individuals per group, and the resulting ceiling on
  # attainable evidence. A group that pools individuals is a legitimate genetic
  # association design; what it cannot do is reach genome-wide significance at
  # small n, so report the ceiling rather than discouraging the analysis.
  subs <- function(tag) setdiff(grep(paste0("^", tag, "\\|"), levs, value = TRUE),
                                c("group1", "group2"))
  n_gen1 <- .count_distinct_genomes(pw$concordance, subs("group1"), same_genome_threshold)
  n_gen2 <- .count_distinct_genomes(pw$concordance, subs("group2"), same_genome_threshold)
  if (is.na(n_gen1)) n_gen1 <- 1L
  if (is.na(n_gen2)) n_gen2 <- 1L

  n_tested <- if (!is.null(info)) nrow(info) else NA_integer_
  min_p <- .min_attainable_p(n_gen1, n_gen2)
  bonf <- if (!is.na(n_tested)) 0.05 / n_tested else NA_real_
  significance_reachable <- if (is.na(min_p) || is.na(bonf)) NA else min_p <= bonf

  het <- (!is.na(w1) && w1 < same_genome_threshold) ||
         (!is.na(w2) && w2 < same_genome_threshold)

  verdict <- if (is.na(between)) {
    "indeterminate"
  } else if (between >= same_genome_threshold) {
    "same_genome"
  } else if (het) {
    "heterogeneous_groups"
  } else if (between <= different_genome_threshold) {
    "different_genomes"
  } else {
    "indeterminate"
  }

  msg <- switch(
    verdict,
    same_genome = paste0(
      "The two groups are the same genome (concordance ", round(between, 3),
      "). Germline presence/absence is structurally empty here: genotype is ",
      "fixed within an individual, so any hit reflects transcript abundance ",
      "or coverage, not genotype. Test gene expression directly instead."),
    heterogeneous_groups = paste0(
      "Groups pool multiple genomes (", n_gen1, " vs ", n_gen2,
      " distinct individuals). This is a genetic association design, not a ",
      "within-library genotype contrast, and it is valid as such -- but the ",
      "evidence available is capped by group size, not by sequencing. The ",
      "smallest attainable two-sided p for a perfectly separating variant is ",
      signif(min_p, 3),
      if (!is.na(bonf)) paste0(", against a Bonferroni threshold of ",
                               signif(bonf, 3), " over ", n_tested, " sites") else "",
      ". ",
      if (isFALSE(significance_reachable))
        paste0("No variant can reach genome-wide significance at this n, so ",
               "rank hits as exploratory candidates rather than thresholding ",
               "them, and let prior biological support drive follow-up. ")
      else "",
      "Note also that unrelated individuals differ in ancestry, which shifts ",
      "allele frequencies genome-wide independently of the grouping variable; ",
      "this test cannot separate the two."),
    different_genomes = paste0(
      "The two groups are distinct genomes (concordance ", round(between, 3),
      ") and each is internally consistent. Presence/absence is a valid ",
      "genotype contrast here."),
    indeterminate = paste0(
      "Genotype relationship could not be resolved",
      if (!is.na(between)) paste0(" (concordance ", round(between, 3),
                                  ", between the same- and different-genome thresholds)")
      else " (too few jointly callable sites)",
      ". Related individuals fall in this band."))

  if (verbose) {
    cat(sprintf("\n=== Genotype concordance ===\n"))
    cat(sprintf("  sites used            : %d autosomal (min total depth %d)\n",
                length(keep), min(dp_tot[keep])))
    cat(sprintf("  cells                 : %d vs %d%s\n", length(i1), length(i2),
                if (drop_n > 0) sprintf("  (%d dropped, sub-group below min_cells)", drop_n) else ""))
    cat(sprintf("  between groups        : %s over %s jointly callable sites\n",
                ifelse(is.na(between), "NA", round(between, 3)),
                ifelse(is.na(n_shared), "NA", n_shared)))
    cat(sprintf("  within group1 (min)   : %s\n", ifelse(is.na(w1), "NA (not splittable)", round(w1, 3))))
    cat(sprintf("  within group2 (min)   : %s\n", ifelse(is.na(w2), "NA (not splittable)", round(w2, 3))))
    cat(sprintf("  distinct genomes      : %d vs %d\n", n_gen1, n_gen2))
    if (n_gen1 > 1 || n_gen2 > 1) {
      cat(sprintf("  best attainable p     : %s (perfect separation at this n)\n", signif(min_p, 3)))
      if (!is.na(bonf)) {
        cat(sprintf("  Bonferroni over %d : %s%s\n", n_tested, signif(bonf, 3),
                    if (isFALSE(significance_reachable))
                      sprintf("  -- unreachable, short by %.0fx", min_p / bonf) else ""))
      }
    }
    cat(sprintf("  verdict               : %s\n", verdict))
    cat(sprintf("  %s\n", msg))
  }

  list(concordance = between,
       n_shared_sites = n_shared,
       within1 = w1,
       within2 = w2,
       n_cells1 = length(i1),
       n_cells2 = length(i2),
       n_sites = length(keep),
       n_genomes1 = n_gen1,
       n_genomes2 = n_gen2,
       min_attainable_p = min_p,
       bonferroni_threshold = bonf,
       significance_reachable = significance_reachable,
       verdict = verdict,
       message = msg,
       pairwise = pw$concordance)
}


#' @title Guard a presence/absence contrast on genotype grounds
#' @name genotype_contrast_guard
#'
#' @description
#' Internal guard run by \code{findSNPsByGroup()}. Resolves the two idents to
#' cell sets using the column \code{aggregateByGroup()} grouped on, falling back
#' to the current project identity, measures genotype concordance, and warns
#' when the contrast is not a genotype contrast. Fewer sites are used
#' than the default for \code{checkGenotypeConcordance()}, since separating
#' ~0.99 from ~0.61 needs far less evidence than the calibration itself.
#'
#' @keywords internal
variantCell$set("private", "genotype_contrast_guard", function(ident.1,
                                                               ident.2 = NULL,
                                                               group_col = NULL,
                                                               donor_type = NULL,
                                                               n_sites = 10000,
                                                               verbose = TRUE) {
  cm <- self$snp_database$cell_metadata

  # The contrast is defined by the column aggregateByGroup() grouped on, which
  # is frequently NOT the project identity - aggregating by "condition" while
  # the identity is "cell_type" is the normal case. Keying on the identity
  # instead meant which(cm$cell_type == "Rejection") came back empty, the
  # concordance call errored, the error was swallowed, and the guard silently
  # did nothing. Prefer the aggregation column and fall back to the identity.
  ident_col <- if (!is.null(group_col) && group_col %in% names(cm)) {
    group_col
  } else {
    self$current_project_ident
  }
  if (is.null(cm) || is.null(ident_col) || !ident_col %in% names(cm)) {
    warning("Genotype concordance check could not run: no usable grouping ",
            "column. Presence/absence results are unguarded - see ",
            "checkGenotypeConcordance().", call. = FALSE)
    return(invisible(NULL))
  }

  # Match the cell set the test itself will use.
  keep <- rep(TRUE, nrow(cm))
  if (!is.null(donor_type) && "donor_type" %in% names(cm)) {
    keep <- !is.na(cm$donor_type) & cm$donor_type %in% donor_type
  }

  g1 <- which(keep & cm[[ident_col]] == ident.1)
  g2 <- if (is.null(ident.2)) which(keep & cm[[ident_col]] != ident.1) else
        which(keep & cm[[ident_col]] == ident.2)

  if (length(g1) == 0 || length(g2) == 0) {
    warning("Genotype concordance check could not run: '",
            if (length(g1) == 0) ident.1 else ident.2,
            "' matched no cells in column '", ident_col,
            "'. Presence/absence results are unguarded.", call. = FALSE)
    return(invisible(NULL))
  }

  res <- tryCatch(
    checkGenotypeConcordance(self$snp_database, g1, g2,
                             n_sites = n_sites, verbose = verbose),
    error = function(e) {
      if (verbose) {
        cat(sprintf("\nGenotype concordance check skipped: %s\n", conditionMessage(e)))
      }
      NULL
    })
  if (is.null(res)) return(invisible(NULL))

  if (res$verdict == "same_genome") {
    warning("Presence/absence contrast is within a single genome (concordance ",
            round(res$concordance, 3), "). Germline genotype is fixed within an ",
            "individual, so this test is structurally empty and any hit ",
            "reflects transcript abundance or coverage. Test gene expression ",
            "directly, or compare different genomes.", call. = FALSE)
  } else if (res$verdict == "heterogeneous_groups") {
    warning("Contrast spans ", res$n_genomes1, " vs ", res$n_genomes2,
            " distinct individuals, so this is a genetic association design ",
            "rather than a within-library genotype contrast. Best attainable ",
            "p at this n is ", signif(res$min_attainable_p, 3),
            if (isFALSE(res$significance_reachable))
              paste0(" against a Bonferroni threshold of ",
                     signif(res$bonferroni_threshold, 3),
                     "; no variant can reach significance, so rank hits as ",
                     "exploratory candidates rather than thresholding them")
            else "",
            ". Ancestry differences between unrelated individuals also shift ",
            "allele frequencies genome-wide and are not separable here.",
            call. = FALSE)
  } else if (res$verdict == "indeterminate") {
    warning("Could not establish whether the two groups are distinct genomes. ",
            "Interpret presence/absence results with care.", call. = FALSE)
  }

  invisible(res)
})
